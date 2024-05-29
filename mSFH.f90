SUBROUTINE CalcSFH(b_or_q, tstart, tfinal, taustar)
  ! --- this subroutine is called in the subroutine 'lum'
  use global_var; use SFHrelated
  implicit none
  INTEGER, INTENT(IN) :: b_or_q
  DOUBLE PRECISION, INTENT(IN) :: tstart, tfinal, taustar
  INTEGER :: i, k, k0
  INTEGER :: i_ts, i_tf, i_tsp1, i_tfm1
  DOUBLE PRECISION :: ab, abinv, telaps
  DOUBLE PRECISION :: taueff, taueffinv
  DOUBLE PRECISION :: zc_SFH, zc_SFH_pre
  DOUBLE PRECISION :: tmp1, tmp2

  mSFH(:,:) = 0.d0

  i_ts = bin_tSFH(tstart); i_tf = bin_tSFH(tfinal)
!!$         print '(2(I2, G10.2), 5G10.2)', i_ts, tstart, i_tf, tfinal,
!!$     $        Zc0, ssp%p(b_or_q), tfinal-tstart, taustar,
!!$     $        ssp%p(b_or_q)*(tfinal-tstart)/taustar
  IF(i_ts >= 1 .and. Mc0 > 0.d0) THEN ! i_tf >= 1 ??
     ab = ssp%alp(b_or_q) + beta; taueff = taustar / ab
     abinv = 1.d0 / ab; taueffinv = 1.d0 / taueff; tmp1 = ssp%p(b_or_q) / taustar
     IF(i_ts == i_tf) THEN
        ! the whole timestep of 'halo%tlife' is covered with a timestep
        !  of 'tSFH'
        telaps = tfinal - tstart
        zc_SFH_pre = Zc0; zc_SFH = CalZcSFH(zc_SFH_pre,telaps,tmp1)
        k0 = bin_metal_SFH(zc_SFH_pre)
        ! SFHbase%Z(k0) <= zc_SFH_pre < SFHbase%Z(k0+1)
        k  = bin_metal_SFH(zc_SFH)
        ! SFHbase%Z(k)  <= zc_SFH     < SFHbase%Z(k+1)
!!$        print '(A, 2(I1, 1X, G10.3))', '1:', k0, zc_SFH_pre, k, zc_SFH
        call CalcSFHInEachMetalBin(i_ts, k0, k, zc_SFH_pre, abinv)
     ELSE
        i_tsp1 = i_ts + 1; telaps = SFHbase%t(i_tsp1) - tstart
        zc_SFH_pre = Zc0; zc_SFH = CalZcSFH(zc_SFH_pre,telaps,tmp1)
        k0 = bin_metal_SFH(zc_SFH_pre); k = bin_metal_SFH(zc_SFH)
!!$        print '(A, 2(I1, 1X, G10.3))', '2:', k0, zc_SFH_pre, k, zc_SFH
        call CalcSFHInEachMetalBin(i_ts, k0, k, zc_SFH_pre, abinv)

        telaps = SFHbase%tstep; i_tfm1 = i_tf - 1
        DO i = i_tsp1, i_tfm1
           zc_SFH_pre = zc_SFH ! metallicity at t = SFHbase%t(i)-tstart
           zc_SFH     = CalZcSFH(zc_SFH_pre, telaps, tmp1)
                        ! metallicity at t = SFHbase%t(i)-tstart+telaps
           k0 = bin_metal_SFH(zc_SFH_pre); k = bin_metal_SFH(zc_SFH)
!!$           print '(A, 2(I1, 1X, G10.3))', '3:',k0,zc_SFH_pre,k,zc_SFH
           tmp2 = exp(-(SFHbase%t(i) - tstart) * taueffinv) * abinv
                  ! SFR at t = SFHbase%t(i)-tstart
           call CalcSFHInEachMetalBin(i, k0, k, zc_SFH_pre, tmp2)
        ENDDO

        telaps = tfinal - SFHbase%t(i_tf)
        IF(telaps > 0.d0) THEN
           zc_SFH_pre = zc_SFH; zc_SFH = CalZcSFH(zc_SFH_pre,telaps,tmp1)
           k0 = bin_metal_SFH(zc_SFH_pre); k = bin_metal_SFH(zc_SFH)
!!$           print '(A, 2(I1, 1X, G10.3))', '4:',k0,zc_SFH_pre,k,zc_SFH
           tmp2 = exp(-(SFHbase%t(i_tf) - tstart) * taueffinv) * abinv
           call CalcSFHInEachMetalBin(i_tf, k0, k, zc_SFH_pre, tmp2)
        ENDIF
     ENDIF
  ENDIF
  call CheckmSFH

CONTAINS
  INTEGER FUNCTION bin_tSFH(time) RESULT(bin)
    DOUBLE PRECISION, INTENT(IN) :: time
    ! tSFH(bin) <= time < tSFH(bin+1)

    bin = int((time - SFHbase%tl) / SFHbase%tstep) + 1
  END FUNCTION bin_tSFH

  INTEGER FUNCTION bin_metal_SFH(metal) RESULT(bin)
    ! --- the function to provide the bin of metallicity for SFH
    !     this function return the integer 'i' so that
    !      SFHbase%Z(i) <= metal < SFHbase%Z(i+1)
    DOUBLE PRECISION, INTENT(IN) :: metal ! metal mass fraction 
                                          !  (not in the unit of [Zsun])

    IF(metal < SFHbase%Z(2)) THEN
       bin = 1
    ELSEIF(metal > SFHbase%Z(paramSFH%Nzp1)) THEN
       bin = paramSFH%Nzp1
    ELSE
       bin = 3
       DO WHILE(bin <= paramSFH%Nzp1 .and. metal >= SFHbase%Z(bin))
          bin = bin + 1
       ENDDO
       bin = bin - 1
    ENDIF

    IF(bin <= paramSFH%Nz) THEN
       IF(metal < SFHbase%Z(bin) .or. metal >= SFHbase%Z(bin+1)) THEN
          print '(G9.2, 1X, I5, 1X, G9.2)', metal, bin, SFHbase%Z(bin)
       ENDIF
    ELSEIF(bin == paramSFH%Nzp1 .and. metal < SFHbase%Z(bin)) THEN
       print '(G9.2, 1X, I5, 1X, G9.2)', metal, bin, SFHbase%Z(bin)
    ENDIF

    IF(bin < 1 .or. bin > paramSFH%Nzp1) THEN
       print '(G9.2, 1X, I5)', metal, bin
    ENDIF
  END FUNCTION bin_metal_SFH

  DOUBLE PRECISION FUNCTION CalZcSFH(zc_pre, tel, tmp) RESULT(zc)
    DOUBLE PRECISION, INTENT(IN) :: zc_pre, tel
    DOUBLE PRECISION, INTENT(IN) :: tmp ! = ssp%p(b_or_q) / taustar

    zc = min(1.d0, zc_pre + tmp * tel)
  END FUNCTION CalZcSFH

  SUBROUTINE CalcSFHInEachMetalBin(itime, iZini, iZfin, Zini, SFR0)
    INTEGER, INTENT(IN) :: itime, iZini, iZfin
    DOUBLE PRECISION, INTENT(IN) :: Zini, SFR0
    INTEGER :: iZ, iZinip1, iZfinm1
    DOUBLE PRECISION :: tmp1, tmp2, tmp3

    IF(SFR0 > 0.d0 .and. itime >= 1 .and. itime <= paramSFH%Ntp1) THEN
       IF(iZfin == iZini) THEN
          ! time-evolved metallicity results in the same metallicity bin
          tmp1 = -taueffinv * telaps
          mSFH(itime, iZfin) = SFR0 * (1.d0 - exp(tmp1))
       ELSE
          ! time-evolved metallicity cross the metallicity bin
          iZinip1 = iZini + 1; iZfinm1 = iZfin - 1
          tmp1 = taustar / ssp%p(b_or_q)
          tmp2 = (SFHbase%Z(iZinip1) - Zini) * tmp1
                 ! requisite time to reach the metallicity of
                 !  SFHbase%Z(iZini+1) from Zini
          tmp2 = tmp2 * taueffinv
          mSFH(itime, iZini) = SFR0 * (1.d0 - exp(-tmp2))
          DO iZ = iZinip1, iZfinm1
             tmp2 = (SFHbase%Z(iZ) - Zini) * tmp1
             ! requisite time to reach the metallicity of
             !  SFHbase%Z(iZ) from Zini
             tmp3 = (SFHbase%Z(iZ+1) - SFHbase%Z(iZ)) * tmp1
             ! requisite time to reach the metallicity of
             !  SFHbase%Z(iZ+1) from SFHbase%Z(iZ)
             mSFH(itime, iZ) = SFR0 * exp(-tmp2 * taueffinv) &
                                 * (1.d0 - exp(-tmp3 * taueffinv))
          ENDDO
          tmp2 = (SFHbase%Z(iZfin) - Zini) * tmp1
                  ! requisite time to reach the metallicity of
                  !  SFHbase%Z(iZfin) from Zini
          tmp3 = telaps - tmp2
                  ! elapsed time to the end of SF from the time at
                  !  which the metallicity reaches to SFHbase%Z(iZfin)
                  !  from Zc0
          mSFH(itime, iZfin) = SFR0 * exp(-tmp2 * taueffinv) &
                                    * (1.d0 - exp(-tmp3 * taueffinv))
       ENDIF
    ENDIF
  END SUBROUTINE CalcSFHInEachMetalBin

  SUBROUTINE CheckmSFH
    INTEGER :: i, j

    DO i = 1, paramSFH%Ntp1
       DO j = 1, paramSFH%Nzp1
          IF(ISNAN(mSFH(i, j)) == .true.) THEN
             mSFH(i, j) = 0.d0; print *, 'NaN'
          ENDIF
       ENDDO
    ENDDO
  END SUBROUTINE CheckmSFH
END SUBROUTINE CalcSFH
