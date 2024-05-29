! ========================================================================
SUBROUTINE SetAndOpenOutputFilesForSN(fbase, file_o, run_redshift)
  use SNrelated; implicit none
  CHARACTER(LEN=*), INTENT(IN) :: fbase, file_o
  INTEGER, INTENT(IN) :: run_redshift ! 1:run at a redshift, 2:run at all redshift
  INTEGER :: ifile, ii, j, ier
  CHARACTER :: ci*1, fbaseSN*500
  CHARACTER*50 :: cerr = '# (SN) SetAndOpenOutputFilesForSFH: fail to'

  allocate(i_SN(paramSN%nfile),    stat=ier); call CheckIerr(ier, &
       trim(cerr)//' allocate: i_SN')
  allocate(fnameSN(paramSN%nfile), stat=ier); call CheckIerr(ier, &
       trim(cerr)//' allocate: fnameSN')

  print '(A)', '# (SN) Output File Names:'
  fbaseSN = trim(fbase)//trim(file_o)//'_'
  DO ifile = 1, paramSN%nfile
     ii = ifile + paramSN%base; write(ci, '(I1)') ifile
     call CheckIONum(ii, 'SetAndOpenOutputFilesForSN')
     fnameSN(ifile) = trim(fbaseSN)//'SN'//ci//'.dat'; i_SN(ifile) = ii

     IF(run_redshift == 1) THEN ! run at a redshift
        open(ii, file=trim(fnameSN(ifile)), status='unknown', iostat=ier)
        call CheckIerr(ier, trim(cerr)//' open')
     ELSEIF(run_redshift == 2) THEN ! run at all redshift
        j = index(fnameSN(ifile), '.dat') - 1
        fnameSN(ifile) = fnameSN(ifile)(1:j)//'_***.dat'
     ENDIF

     print '(A)', '  --- '//trim(fnameSN(ifile))
  ENDDO
END SUBROUTINE SetAndOpenOutputFilesForSN
! ========================================================================
SUBROUTINE AllocateSNRelated(run_redshift)
  use global_var; use SNrelated; implicit none
  INTEGER, INTENT(IN) :: run_redshift ! 1:run at a redshift, 2:run at all redshift
  INTEGER :: i, j, k, i_dtd, itime, ndata, ier, i_file
  CHARACTER*100 :: fread(18) = &
       (/CHARACTER*100:: 'cha_gr_0.0004', 'cha_gr_0.004', 'cha_gr_0.008', 'cha_gr_0.02', & !1-4
       'cha_gr_0.05', 'delta0.5Gyra0.001', 'delta1Gyra0.001', 'delta3Gyra0.001', & !5-8
       'power0.05Gyr-1a0.001', 'power0.1Gyr-1a0.001', 'power0.5Gyr-1a0.001', & !9-11
       'power1Gyr-1a0.001', 'power3Gyr-1a0.001', 'gauss0.1Gyr0.3dexa0.001', & ! 12-14
       'gauss0.5Gyr0.3dexa0.001', 'gauss1Gyr0.3dexa0.001', 'gauss3Gyr0.3dexa0.001', & ! 15-17
       'flata0.001'/) ! 18
  CHARACTER*100 :: cerr = '# (SN) AllocateSNRelated: fail to'
  CHARACTER*500 :: filename
  DOUBLE PRECISION, ALLOCATABLE :: time(:), DTDorig(:), IDTD(:)
  DOUBLE PRECISION, PARAMETER :: EPS = 1.d-20
  DOUBLE PRECISION :: x, tmp, time_th
  DOUBLE PRECISION :: LinearInterp ! function


  call SetBaseQuantitiesForSN
  call AllocateDelayTimeDistribution
  allocate(IDTD(paramSN%Ntp1), stat=ier); call CheckIerr(ier, &
       trim(cerr)//' allocate IDTD')
  call CheckIONum(paramSN%base, 'AllocateSNRelated')
  filename = 'TablesDTD/DTD.dat'
  open(paramSN%base, file = trim(filename), status='unknown', &
       iostat=ier); call CheckIerr(ier, trim(cerr)//' open file = '//trim(filename))
  i_file = paramSN%base + paramSN%nfile + 1
  call CheckIONum(i_file, 'AllocateSNRelated')
  DO i_dtd = 1, paramSN%N_DTD
     IDTD(:) = 0.d0
     write(paramSN%base, '(A, I2, A)') '# index ', i_dtd-1, &
          ': '//trim(fread(i_dtd))//'.dat'
     ndata = 0
     filename = repeat(' ', 500)
     filename = 'TablesDTD/'//trim(fread(i_dtd))//'.dat'
     open(i_file, file = trim(filename), status='old', iostat=ier); call &
       CheckIerr(ier, trim(cerr)//' open file = '//trim(filename))
     DO ndata = 1, 100000
        read(i_file, *, end=100) x, x
     ENDDO
100  close(i_file); ndata = ndata - 1

     allocate(time(ndata),    stat=ier); call CheckIerr(ier, &
          trim(cerr)//' allocate time')
     allocate(DTDorig(ndata), stat=ier); call CheckIerr(ier, &
          trim(cerr)//' allocate DTDorig')
     open(i_file, file = trim(filename), status='old', iostat=ier); call &
       CheckIerr(ier, trim(cerr)//' open file = '//trim(filename))
     DO j = 1, ndata
        read(i_file, *) time(j), DTDorig(j)
     ENDDO
     time(:) = time(:)  / param%th_yr ! [yr] --> [hubble time]
     close(i_file)

     i = 1; time_th = SubstAdeqTimeThreshold(i_dtd)
     DO itime = 1, paramSN%Nt
        IF(SNbase%t_lk(itime) - time_th < 1.d-10) THEN
           IDTD(itime) = 0.d0
        ELSE
           k = 1
           DO WHILE(time(k) < SNbase%t_lk(itime) .and. k+1 < ndata)
              k = k + 1
           ENDDO
           ! time(k-1) < SNbase%t_lk(itime) <= time(k)

           IF(k > ndata) THEN
              IDTD(itime) = DTDorig(ndata)
           ELSEIF(k == 1) THEN
              IDTD(itime) = DTDorig(1)
           ELSE
              IF(SNbase%t_lk(itime+1) - time_th < 1.d-10) THEN
                 x   = time_th
                 tmp = LinearInterp(x, time(k), time(k+1), DTDorig(k), DTDorig(k+1), &
                                    1.d0/(time(k)-time(k+1)))
                 IDTD(itime) = LinearInterp(SNbase%t_lk(itime), time(k), x, &
                                             DTDorig(k), tmp, 1.d0/(time(k)-time_th))
              ELSE
                 IDTD(itime) = LinearInterp(SNbase%t_lk(itime), time(k), time(k-1), &
                                             DTDorig(k), DTDorig(k-1), &
                                             1.d0/(time(k)-time(k-1)))
              ENDIF
           ENDIF
        ENDIF
     ENDDO
     IDTD(paramSN%Ntp1) = 0.d0
     deallocate(time,    stat=ier); call CheckIerr(ier, &
          trim(cerr)//' deallocate time')
     deallocate(DTDorig, stat=ier); call CheckIerr(ier, &
          trim(cerr)//' deallocate DTDorig')

     DO itime = 1, paramSN%Nt
        DTD(i_dtd)%r(itime) = (IDTD(itime) - IDTD(itime+1)) / SNbase%tstep_yr
                               ! mean SN rate [yr^-1] in the timestep
        IF(DTD(i_dtd)%r(itime) <= EPS &
             .or. SNbase%t_lk(itime) - time_th < 1.d-10) THEN
           DTD(i_dtd)%flag(itime) = 0; DTD(i_dtd)%r(itime) = 0.d0
        ELSE
           DTD(i_dtd)%flag(itime) = 1
        ENDIF
        write(paramSN%base, '(2I6, 4(X, E12.5), X, I1)') i_dtd, itime, &
             SNbase%t_univ(itime)*param%th, SNbase%t_lk(itime)*param%th, &
             IDTD(itime), DTD(i_dtd)%r(itime), DTD(i_dtd)%flag(itime)

        DTD(i_dtd)%r(itime) = DTD(i_dtd)%r(itime) * param%munit
     ENDDO
     write(paramSN%base, *); write(paramSN%base, *)
  ENDDO
  deallocate(IDTD, stat=ier); call CheckIerr(ier, trim(cerr)//' deallocate IDTD')
  close(paramSN%base, iostat=ier); call CheckIerr(ier, trim(cerr)//' close file = '// &
       trim(filename))
CONTAINS
  SUBROUTINE SetBaseQuantitiesForSN
    use SFHrelated
    INTEGER :: i, ier
    DOUBLE PRECISION :: Cube

    constSN%hinv = 1.d0 / param%h; constSN%hCube = Cube(param%h)
    constSN%corr = 5.d0 * log10(param%h)

    paramSN%Nt = paramSFH%Nt; paramSN%Nz = paramSFH%Nz
    paramSN%Ntp1 = paramSN%Nt + 1; paramSN%Nzp1 = paramSN%Nz + 1
    paramSN%thMstar = 1.d+6 ! [Msun]

    allocate(SNbase%t_univ(paramSN%Ntp1), SNbase%t_lk(paramSN%Ntp1), stat=ier)
    allocate(SNbase%Z(paramSN%Nzp1), stat=ier)
    IF(ier /= 0) THEN
       write(*, *) 'SetBaseQuantitiesForSN: fail to allocate stat=',&
            ier; stop
    ENDIF
    SNbase%Z(:) = SFHbase%Z(:); SNbase%t_univ(:) = SFHbase%t(:)
    SNbase%t_lk(1) = SFHbase%tu - SFHbase%tl
    DO i = 2, paramSN%Nt
       SNbase%t_lk(i) = SFHbase%tu - SFHbase%t(i)
    ENDDO
    SNbase%t_lk(paramSN%Ntp1) = 0.d0

    SNbase%tl = SFHbase%tl; SNbase%tu = SFHbase%tu
    SNbase%tstep_yr = SFHbase%tstep_yr
  END SUBROUTINE SetBaseQuantitiesForSN

  SUBROUTINE AllocateDelayTimeDistribution
    INTEGER :: i, j, ier
    CHARACTER*100 :: cerr = '# (SN) AllocateDelayTimeDistribution: '//&
         'fail to allocate:'
    CHARACTER*100 :: name(18) = (/CHARACTER*100::  &
         'Greggio & Renzini83 for Z/Zsun=1/50', &
         'Greggio & Renzini83 for Z/Zsun=1/5', &
         'Greggio & Renzini83 for Z/Zsun=1/2.5', &
         'Greggio & Renzini83 for Z/Zsun=1', &
         'Greggio & Renzini83 for Z/Zsun=2.5', &
         'delta function at 0.5 Gyr', &
         'delta function at 1 Gyr', &
         'delta function at 3 Gyr', &
         'power-law (t^-1) from 0.05 Gyr', &
         'power-law (t^-1) from 0.1 Gyr', &
         'power-law (t^-1) from 0.5 Gyr', &
         'power-law (t^-1) from 1 Gyr', &
         'power-law (t^-1) from 3 Gyr', &
         'Gaussian (sigma=0.3 dex, <t>=0.1 Gyr)', &
         'Gaussian (sigma=0.3 dex, <t>=0.5 Gyr)', &
         'Gaussian (sigma=0.3 dex, <t>=1 Gyr)', &
         'Gaussian (sigma=0.3 dex, <t>=3 Gyr)', &
         'flat (1Myr-15Gyr)' &
         /)

    allocate(R_SN(paramSN%N_DTD), stat=ier); call CheckIerr(ier, &
         trim(cerr)//' R_SN')
    allocate(CSNR(paramSN%N_DTD), stat=ier); call CheckIerr(ier, &
         trim(cerr)//' CSNR')
    CSNR = 0.d0
    allocate(DTD(paramSN%N_DTD),  stat=ier); call CheckIerr(ier, &
         trim(cerr)//' DTD')
    DO i = 1, paramSN%N_DTD
       DTD(i)%name = trim(name(i))
       allocate(DTD(i)%flag(paramSN%Ntp1), stat=ier); call CheckIerr(ier, &
            trim(cerr)//' DTD%flag')
       allocate(DTD(i)%r(paramSN%Ntp1),    stat=ier); call CheckIerr(ier, &
            trim(cerr)//' DTD%r')
       DTD(i)%flag(:) = 0; DTD(i)%r(:) = 0.d0
    ENDDO
  END SUBROUTINE AllocateDelayTimeDistribution

  DOUBLE PRECISION FUNCTION SubstAdeqTimeThreshold(i_dtd) RESULT(time_th)
    INTEGER, INTENT(IN) :: i_dtd

    time_th = -1.d0
    IF(i_dtd == 6 .or. i_dtd == 11) THEN
       time_th = 0.5d+9
    ELSEIF(i_dtd == 7 .or. i_dtd == 12) THEN
       time_th = 1.d+9
    ELSEIF(i_dtd == 8 .or. i_dtd == 13) THEN
       time_th = 3.d+9
    ELSEIF(i_dtd == 9) THEN
       time_th = 0.05d+9
    ELSEIF(i_dtd == 10) THEN
       time_th = 0.1d+9
    ENDIF
    time_th = time_th / param%th_yr
  END FUNCTION SubstAdeqTimeThreshold
END SUBROUTINE AllocateSNRelated
! ========================================================================
SUBROUTINE WriteCaptionsForSN(ilog)
 use global_var;  use SNrelated; implicit none
  INTEGER, INTENT(IN) :: ilog
  INTEGER :: i, ii, ier, i0
  CHARACTER*100, ALLOCATABLE :: ref(:)
  CHARACTER*100 :: cerr = '# (SN) WriteCaptionsForSN: fail to'
  DOUBLE PRECISION :: redshift

  redshift = param%zsp1 - 1.d0
  DO i = 1, paramSN%nfile
     write(i_SN(i), '(A, G10.2, A, F5.2)') '# Supernova Related '//&
          'Information of galaxies with Mstar >= ', paramSN%thMstar, &
          ' [Msun] at z=', redshift
     write(i_SN(i), '(A, F7.4, A)') &
          '# The age of universe at this redshift = ', const%toutGyr, &
          ' [Gyr]'
     write(i_SN(i), '(4(A, F6.3), A)') &
          '# The cosmological parameters used are h =', param%h, &
          ', Omega_M =', param%OMEGA0, ', Omega_b =', param%OMEGA, &
          ', and Omega_L =', 1.d0-param%OMEGA0, ' (flat univ.)'
     write(i_SN(i), '(A)') '# '//trim(ssp%filename(1))//&
          ' ssp file is used for starburst'
     write(i_SN(i), '(A)') '# '//trim(ssp%filename(2))//&
          ' ssp file is used for quiescent'
     write(i_SN(i), '(A, F11.4, A, I5)') '# Time step for SFH = ', &
          SNbase%tstep_yr * 1.d-6, ' [Myr]: Number of time-bin = ', &
          paramSN%Nt
  ENDDO


  ! --- write captions for i_SN(1)
  paramSN%iband = 2 ! exclude B and Br
  i0 = 11
  i = (param%nwave - paramSN%iband) + i0
  write(i_SN(1), '(4(A, I2), A)') '# (1)ID(1:central, others:satellite) '//&
       '(2)morphology(1:E, 2:S0, 3:S) (3)1:quiescent,2:starburst '//&
       '(4)weight[h_70^3/Mpc^3] (5)Mstar[Msun] (6)reff[kpc/h_70] '//&
       '(7)dust opacity in V-band (8)mass-weighted age[yr] '//&
       '(9)mean SFR during 10 Myr [Msun/yr] (10)Mhost[Msun] '//&
       '(11)Mprog[Msun] (12-', i, ')Abs.Mag-5logh_70[ABmag] (', &
       i+1, '-', i+paramSN%N_DTD, ')R_SN[yr^-1] (', i+paramSN%N_DTD+1, &
       ')Mstar(SFH)[Msun]'
  write(i_SN(1), '(A, $)') '# bandname:'
  DO i = 1+paramSN%iband, param%nwave
     ii = (i - paramSN%iband) + i0
     write(i_SN(1), '(A, I2, A, $)') ' (', ii, &
          ')'//trim(ssp%bandname(param%iwave(i)))
  ENDDO
  write(i_SN(1), *)
  write(i_SN(1), '(A, $)') '# DTD type:'
  DO i = 1, paramSN%N_DTD
     ii = i + (param%nwave - paramSN%iband) + i0
     write(i_SN(1), '(A, I2, A, $)') ' (', ii, ')'//trim(DTD(i)%name)
  ENDDO
  write(i_SN(1), *)


  ! --- write captions for i_SN(2)
  write(i_SN(2), '(A, I2, A)') '# (1)z (2)t_univ[Gyr] (3-', &
       paramSN%N_DTD+2, ')CSNR[h_70^3/Mpc^3/yr]'
  write(i_SN(2), '(A, $)') '#'
  DO i = 1, paramSN%N_DTD
     ii = i + 2
     write(i_SN(2), '(A, I2, A, $)') ' (', ii, ')'//trim(DTD(i)%name)
  ENDDO
  write(i_SN(2), *)

  ! --- write information as log into logfile
  write(ilog, *); write(ilog, *)
  write(ilog, '(A)') '# +-------------------------------------------+'
  write(ilog, '(A)') '# |          SN related parameters            |'
  write(ilog, '(A)') '# +-------------------------------------------+'
  write(ilog, '(A, I3)') '# paramSN%Nt: ', paramSN%Nt
  write(ilog, '(A, I3)') '# paramSN%Nz: ', paramSN%Nz
  write(ilog, '(A, I3)') '# paramSN%N_DTD: ', paramSN%N_DTD
  write(ilog, '(A)') '# filenames used in the SN related calculation:'
  allocate(ref(paramSN%nfile), stat=ier); call CheckIerr(ier, &
       trim(cerr)//' allocate: ref')
  ref(1) = ' = SN rate + SED of each galaxy'
  ref(2) = ' = cosmic SN rate density'
  DO i = 1, paramSN%nfile
     write(ilog, '(A, I2, A)') '# --- ', i, ':'//trim(fnameSN(i))//trim(ref(i))
  ENDDO
  deallocate(ref, stat=ier); call CheckIerr(ier, trim(cerr)//' deallocate: ref')
END SUBROUTINE WriteCaptionsForSN
! ========================================================================
SUBROUTINE DeallocateDelayTimeDistribution
  use SNrelated
  INTEGER :: i, ier
  CHARACTER*100 :: cerr = '# (SN) DeallocateDelayTimeDistribution: '//&
       'fail to deallocate:'

  DO i = 1, paramSN%N_DTD
     deallocate(DTD(i)%flag, stat=ier); call CheckIerr(ier, &
          trim(cerr)//' DTD%flag')
     deallocate(DTD(i)%r,    stat=ier); call CheckIerr(ier, &
          trim(cerr)//' DTD%r')
  ENDDO
  deallocate(DTD, stat=ier); call CheckIerr(ier, trim(cerr)//' DTD')
END SUBROUTINE DeallocateDelayTimeDistribution
! ========================================================================
SUBROUTINE CalAndWriteSNrate(id, mord, flag_burst, inv_V, Mstar, r_eff, taud_V, AgeM, &
     mSFR, Mhost, magd)
  use global_var; use SFHrelated; use SNrelated; implicit none
  INTEGER, INTENT(IN) :: id, mord, flag_burst ! flag_burst=0:quiescent, 1:starburst
  DOUBLE PRECISION, INTENT(IN) :: inv_V ! weight of galaxy [h^3/Mpc^3]
  DOUBLE PRECISION, INTENT(IN) :: Mstar ! stellar mass of galaxy [Msun]
  DOUBLE PRECISION, INTENT(IN) :: r_eff ! effective radius of galaxy [kpc/h]
  DOUBLE PRECISION, INTENT(IN) :: taud_V ! V-band opacity of galaxy
  DOUBLE PRECISION, INTENT(IN) :: AgeM ! mass-weighted age [yr]
  DOUBLE PRECISION, INTENT(IN) :: mSFR ! mean SFR during the past 10 Myr [Msun/yr]
  DOUBLE PRECISION, INTENT(IN) :: magd(param%nwave) ! dust-extincted magnitudes
  DOUBLE PRECISIOn, INTENT(IN) :: Mhost ! host halo mass [Msun]
  DOUBLE PRECISION, PARAMETER :: EPS = 1.d-30
  INTEGER :: i, i_dtd, itime, iZ, iSFH, sftype ! sftype=1:quiescent, 2:starburst
  DOUBLE PRECISION :: tmp, Mstar_SFH

  sftype = 1 ! quiescent
  IF(flag_burst == 1) sftype = 2 ! starburst

  ! --- write physical and observable quantities of the galaxy into output file
  write(i_SN(1), '(I8, 2I2, 8E12.4, $)') id, mord, sftype, &
       inv_V * constSN%hCube, Mstar, r_eff * constSN%hinv, taud_V, &
       AgeM, mSFR, Mhost, gal(id)%Mhalo * param%munit
  DO i = 1+paramSN%iband, param%nwave
     write(i_SN(1), '(F8.3, $)') magd(i) + constSN%corr
  ENDDO

  Mstar_SFH = 0.d0
  DO itime = 1, paramSN%Nt
     DO iZ = 1, paramSN%Nzp1
        IF(ISNAN(gal(id)%SFH(itime, iZ)) .eqv. .false.) &
             Mstar_SFH = Mstar_SFH + gal(id)%SFH(itime, iZ)
     ENDDO
  ENDDO
  DO i_dtd = 1, paramSN%N_DTD
     R_SN(i_dtd) = 0.d0
     DO itime = 1, paramSN%Nt
        IF(DTD(i_dtd)%flag(itime) == 1) THEN
           tmp = 0.d0
           IF(i_dtd == 1) THEN
              DO iZ = 1, 2
                 IF(ISNAN(gal(id)%SFH(itime, iZ)) .eqv. .false.) &
                      tmp = tmp + gal(id)%SFH(itime, iZ)
              ENDDO
           ELSEIF(i_dtd == 2) THEN
              IF(ISNAN(gal(id)%SFH(itime, 3)) .eqv. .false.) &
                   tmp = gal(id)%SFH(itime, 3)
           ELSEIF(i_dtd == 3) THEN
              IF(ISNAN(gal(id)%SFH(itime, 4)) .eqv. .false.) &
                   tmp = gal(id)%SFH(itime, 4)
           ELSEIF(i_dtd == 4) THEN
              IF(ISNAN(gal(id)%SFH(itime, 5)) .eqv. .false.) &
                   tmp = gal(id)%SFH(itime, 5)
           ELSEIF(i_dtd == 5) THEN
              IF(ISNAN(gal(id)%SFH(itime, 6)) .eqv. .false.) &
                   tmp = gal(id)%SFH(itime, 6)
           ELSE
              DO iZ = 1, paramSN%Nzp1
                 IF(ISNAN(gal(id)%SFH(itime, iZ)) .eqv. .false.) THEN
                    tmp = tmp + gal(id)%SFH(itime, iZ)
                 ENDIF
              ENDDO
           ENDIF
           R_SN(i_dtd) = R_SN(i_dtd) + tmp * DTD(i_dtd)%r(itime)
        ENDIF
     ENDDO
     IF(R_SN(i_dtd) < EPS) R_SN(i_dtd) = 0.d0
     write(i_SN(1), '(E12.4, $)') R_SN(i_dtd)
  ENDDO
  write(i_SN(1), '(E12.4)') Mstar_SFH * param%munit

  ! --- adding the contribution of the galaxy into cosmic SN rate density
  CSNR(:) = CSNR(:) + R_SN(:) * inv_V
END SUBROUTINE CalAndWriteSNrate
! ========================================================================
SUBROUTINE WriteCosmicSNRateDensity
  use global_var; use SNrelated; implicit none
  INTEGER :: i, ier
  CHARACTER*100 :: cerr = '# (SN) WriteCosmicSNRateDensity: fail to deallocate'

  write(i_SN(2), '(F6.2, F10.6, $)') &
       param%zsp1-1.d0, SNbase%tu*param%th_yr*1.d-9
  DO i = 1, paramSN%N_DTD
     write(i_SN(2), '(G11.4, $)') CSNR(i)*constSN%hCube
  ENDDO
  write(i_SN(2), *)

  deallocate(R_SN, stat=ier); call CheckIerr(ier, trim(cerr)//' R_SN')
  deallocate(CSNR, stat=ier); call CheckIerr(ier, trim(cerr)//' CSNR')
END SUBROUTINE WriteCosmicSNRateDensity
! ========================================================================
SUBROUTINE CloseFilesForSN
  use SNrelated; implicit none
  INTEGER :: i, ii, ier
  CHARACTER*100 :: cerr = '# (SN) CloseFilesForSN: fail to'

  deallocate(i_SN)
  DO i = 1, paramSN%nfile
     ii = i + paramSN%base
     close(ii, iostat=ier); call CheckIerr(ier, trim(cerr)//' close file = '//&
          trim(fnameSN(i)))
  ENDDO
  deallocate(fnameSN, stat=ier); call CheckIerr(ier, trim(cerr)//' deallocate fnameSN')
END SUBROUTINE CloseFilesForSN
! ========================================================================
