!!$========================================================================
SUBROUTINE SetAndOpenOutputFilesForSFH(fbase, file_o, run_redshift)
  use SFHrelated
  implicit none
  CHARACTER(LEN=*), INTENT(IN) :: fbase, file_o
  INTEGER, INTENT(IN) :: run_redshift ! 1:run at a redshift, 2:run at all redshift
  INTEGER :: ifile, ii, j, ier
  CHARACTER :: ci*1, fbaseSFH*500
  CHARACTER*50 :: cerr = '# (SFH) SetAndOpenOutputFilesForSFH: fail to '


  allocate(i_SFH(paramSFH%nfile),    stat=ier); call CheckIerr(ier, &
       trim(cerr)//' allocate: i_SFH')
  allocate(fnameSFH(paramSFH%nfile), stat=ier); call CheckIerr(ier, &
       trim(cerr)//' allocate: fnameSFH')

  print '(A)', '# (SFH) Output File Names:'
  fbaseSFH = trim(fbase)//trim(file_o)//'_'
  DO ifile = 1, paramSFH%nfile
     ii = ifile + paramSFH%base; write(ci, '(I1)') ifile
     call CheckIONum(ii, 'SetAndOpenOutputFilesForSFH')
     fnameSFH(ifile) = trim(fbaseSFH)//'SFH'//ci//'.dat'
     i_SFH(ifile)    = ii

     IF(run_redshift == 1) THEN ! run at a redshift
        open(ii, file=trim(fnameSFH(ifile)), status='unknown', iostat=ier)
        call CheckIerr(ier, trim(cerr)//' open')
     ELSEIF(run_redshift == 2) THEN ! run at all redshift
        j = index(fnameSFH(ifile), '.dat') - 1
        fnameSFH(ifile) = fnameSFH(ifile)(1:j)//'_***.dat'
     ENDIF

     print '(A)', '  --- '//trim(fnameSFH(ifile))
  ENDDO
END SUBROUTINE SetAndOpenOutputFilesForSFH
!!$========================================================================
SUBROUTINE WriteCaptionsForSFH(ilog)
  use global_var; use SFHrelated; implicit none
  INTEGER, INTENT(IN) :: ilog
  INTEGER :: i, j, k, ier, iband, i1, itime
  CHARACTER*100, ALLOCATABLE :: ref(:)
  CHARACTER*50 :: cerr = '# (SFH) WriteCaptionsForSFH: fail to '
  CHARACTER(LEN=20) :: type_mag(2) = (/CHARACTER(LEN=20):: '[Vega mag]', '[AB mag]'/)
  DOUBLE PRECISION :: redshift


  redshift = param%zsp1 - 1.d0
  DO i = 1, paramSFH%nfile
     write(i_SFH(i), '(2(A, F7.4))') '# SFH calculation in the redshift '// &
          'range of z= ', SFHbase%zp1l-1.d0, '--', SFHbase%zp1u-1.d0
     write(i_SFH(i), '(A, I5, A)') '# Time step is linearly separated with '// &
          't_univ(z)/Nt (Nt =', paramSFH%Nt, ')'
     write(i_SFH(i), '(A, F9.2, A)') '# Time Step for SFH = ', &
          SFHbase%tstep_yr * 1.d-6, ' [Myr]'
     write(i_SFH(i), '(A, F9.4, A)') '# t_hubble(z) = ', const%toutGyr, ' [Gyr]'
  ENDDO

  ! --- *SFH1.dat: cosmic SFH
  write(i_SFH(1), '(A, F5.2)') '# Cosmic SFH at z = ', redshift
  i = paramSFH%N_AV - 2; j = paramSFH%N_AV + 2; k = paramSFH%N_AV - 1
  write(i_SFH(1), '(3(A, I2))') '# (1)Hubble time[yr] (2)look back time '// &
       'from z [yr] (3+i)CSFH[Msun/yr/Mpc^3] with stepAV*i <= Av '// &
       '< stepAV*(i+1) for i = 0-', i, ' (', j, ')CSFH[Msun/'// &
       'yr/Mpc^3] with Av >= stepAV*', k
  write(i_SFH(1), '(A, F5.2, A)') '# Step of Av is ', paramSFH%stepAV, &
       ' [mag]'

  ! --- *SFH{2,3,4}.dat: physical quantities and SFH of each galaxy
  write(i_SFH(2), '(A)') '# physical quantities of each galaxy'
  write(i_SFH(2), '(A, F6.2, A, $)') '# (1)ID (2)ID of central galaxy '//&
       '(3,4)mor,mor_d[1:E,2:S0,3:S] (5)flag_c (6)flag_burst '//&
       '(7)iforest (8)ID(host) (9)mpi(gal) (10)mpi(host) (11)f_des '//&
       '(12,13,14)Mstarb,Mstard,Mcool[Msun] (15,16)Mhalo^{host,prog}'//&
       '[Msun] (17,18)Vcirc^{host,prog}[km/s] (19)SFR[Msun/yr] '//&
       '(20)mean SFR during the past ', const%t0yr * 1.d-6, ' Myr'//&
       '[Msun/yr] (21)Mbh[Msun] (22,23,24)mass-weighted mean stellar '//&
       'metallicity(bulge,disk,all) '
  DO iband = 2, param%nwave
     i1 = 25 + 2 * (iband - 2)
     IF(i1 < 99) THEN
        write(i_SFH(2), '(2(A, I2), A, $)') '(', i1, ',', i1+1, &
             ')'//trim(ssp%bandname(param%iwave(iband)))//'-band mag (w/o,'//&
             'w/ dust)'//trim(type_mag(param%type_mag))//' '
     ELSEIF(i1 == 99) THEN
        write(i_SFH(2), '(A, I2, A, I3, A, $)') '(', i1, ',', i1+1, &
             ')'//trim(ssp%bandname(param%iwave(iband)))//'-band mag (w/o,'//&
             'w/ dust)'//trim(type_mag(param%type_mag))//' '
     ELSE
        write(i_SFH(2), '(2(A, I3), A, $)') '(', i1, ',', i1+1, &
             ')'//trim(ssp%bandname(param%iwave(iband)))//'-band mag (w/o,'//&
             'w/ dust)'//trim(type_mag(param%type_mag))//' '
     ENDIF
  ENDDO
  write(i_SFH(2), *)
  write(i_SFH(3), '(A)') '# Timeslices for SFH'
  write(i_SFH(3), '(A)') '# (1)t_univ[yr] (2)t_lookback[yr]'
  DO itime = 1, paramSFH%Nt
     write(i_SFH(3), '(2G13.5)') SFHbase%t_yr(itime), SFHbase%tlk_yr(itime)
  ENDDO
  write(i_SFH(4), '(A)') '# SFH of each galaxy'
  write(i_SFH(4), '(A)') '# (1)SFH(t)[Msun/yr]'
!!$  write(i_SFH(4), '(A, $)') '# (1)t_univ[yr] (2)t_lookback[yr] (3-'
!!$  IF(paramSFH%Nzp1+2 < 10) THEN
!!$     write(i_SFH(4), '(I1, A)') paramSFH%Nzp1+2, ')SFR(Z)[Msun/yr]'
!!$  ELSE
!!$     write(i_SFH(4), '(I2, A)') paramSFH%Nzp1+2, ')SFR(Z)[Msun/yr]'
!!$  ENDIF
!!$  DO i = 1, paramSFH%Nzp1
!!$     IF(i+2 < 10) THEN
!!$        write(i_SFH(4), '(A, I1, A, $)') '# (', i+2, ')'
!!$     ELSE
!!$        write(i_SFH(4), '(A, I2, A, $)') '# (', i+2, ')'
!!$     ENDIF
!!$     IF(i == paramSFH%Nzp1) THEN
!!$        write(i_SFH(4), '(A, G9.2)') 'Z > ', SFHbase%Z(i)
!!$     ELSE
!!$        write(i_SFH(4), '(2(A, G9.2))') &
!!$               'Z = ', SFHbase%Z(i), ' --', SFHbase%Z(i+1)
!!$     ENDIF
!!$  ENDDO

  ! --- *SFH{5,6}.dat: information for debug
  DO i = 5, 6
     write(i_SFH(i), '(A)') '# Information of the galaxies '// &
          'with abs[Mstar(Mitaka)-Mstar(SFH)]/Mstar(Mitaka) '// &
          '> 1.0 and Mstar(Mitaka) >= 10^5 Msun'
  ENDDO
  write(i_SFH(5), '(A)') '# physical quantities of each galaxy'
  write(i_SFH(6), '(A)') '# SFH of each galaxy'

  write(ilog, *); write(ilog, *)
  write(ilog, '(A)') '# +--------------------------------------------+'
  write(ilog, '(A)') '# |          SFH related parameters            |'
  write(ilog, '(A)') '# +--------------------------------------------+'
  write(ilog, '(A, I3)') '# paramSFH%Nt:   ', paramSFH%Nt
  write(ilog, '(A, I3)') '# paramSFH%Nz:   ', paramSFH%Nz
  write(ilog, '(A, I3)') '# paramSFH%N_AV: ', paramSFH%N_AV
  write(ilog, '(A)') '# filenames used in the SFH related calculation:'
  allocate(ref(paramSFH%nfile), stat=ier)
  call CheckIerr(ier, trim(cerr)//' allocate: ref')
  ref(1) = ' = CSFRD'
  ref(2) = ' = Physical Quantities of Each Galaxy'
  ref(3) = ' = Timeslices for SFH'
  ref(4) = ' = SFH of each galaxy'
  ref(5) = ' = for check in the calc of SFRD: Physical Quantities'
  ref(6) = ' = for check in the calc of SFRD: SFH'
  DO i = 1, paramSFH%nfile
     write(ilog, '(A, I2, A)') '# --- ', i, ':'//trim(fnameSFH(i))//trim(ref(i))
  ENDDO
  deallocate(ref, stat=ier); call CheckIerr(ier, trim(cerr)//' deallocate: ref')
END SUBROUTINE WriteCaptionsForSFH
!!$========================================================================
SUBROUTINE AllocateSFHRelated(run_redshift)
  use global_var; use SFHrelated
  implicit none
  INTEGER, INTENT(IN) :: run_redshift
  INTEGER :: i, ier, Ntp1, Nzp1
  CHARACTER*50 :: cerr = '# (SFH) AllocateSFHRelated: fail to allocate:'
  DOUBLE PRECISION, PARAMETER :: tlk_max = 10.d+9 ! the maximum look-back time in yr
  DOUBLE PRECISION :: tmp
  DOUBLE PRECISION :: z2t ! function

  paramSFH%loop = 0

  SFHbase%zp1l = param%zsp1;    SFHbase%tu = z2t(SFHbase%zp1l)
  SFHbase%zp1u = mrgp%zp1ar(1); SFHbase%tl = z2t(SFHbase%zp1u)
  IF(paramSFH%div_type == 1) THEN ! equal number of time-bin among different redshift
     ! linearly separated timestep
     SFHbase%tstep    = (SFHbase%tu - SFHbase%tl) / dble(paramSFH%Nt)
     SFHbase%tstep_yr = SFHbase%tstep * param%th_yr
  ELSEIF(paramSFH%div_type == 2) THEN ! equally spaced timestep among different redshift
     ! linearly separated timestep
     SFHbase%tstep_yr = 100.d+6 ! [yr] 10, 50, 100, 500 Myr
     SFHbase%tstep    = SFHbase%tstep_yr / param%th_yr

     IF(SFHbase%tu - SFHbase%tl < tlk_max / param%th_yr) THEN
        paramSFH%Nt = 1 + int((SFHbase%tu - SFHbase%tl) / SFHbase%tstep)
     ELSE ! lookback time is set to be tlk_max
        paramSFH%Nt = 1 + int(tlk_max / SFHbase%tstep_yr)
     ENDIF
     SFHbase%tl   = SFHbase%tu - dble(paramSFH%Nt - 1) * SFHbase%tstep
     SFHbase%zp1u = SFHbase%zp1u * &
          (SFHbase%tl / z2t(SFHbase%zp1u))**(-2.d0/3.d0)
     print '(5(A, G12.4), A)', '  zu =', SFHbase%zp1u-1.d0, &
          ', tl(given) =', SFHbase%tl * param%th_yr * 1.d-9, ' Gyr'// &
          ', tl(zu) = ', z2t(SFHbase%zp1u) * param%th_yr * 1.d-9, ' Gyr'// &
          ', tl(z=', mrgp%zp1ar(1)-1.d0, ') = ', &
          z2t(mrgp%zp1ar(1)) * param%th_yr * 1.d-9, ' Gyr'
  ENDIF
  print '(A, F10.4, A, I5)', '# (SFH) Time Step for SFH = ', &
       SFHbase%tstep_yr * 1.d-6, ' [Myr] with paramSFH%Nt =', paramSFH%Nt

  paramSFH%SFRunit = param%munit / SFHbase%tstep_yr
                     ! conversion factor for SFR from the nuGC unit to Msun/yr

  paramSFH%Ntp1 = paramSFH%Nt + 1; Ntp1 = paramSFH%Ntp1
  allocate(SFHbase%t(Ntp1), stat=ier); call CheckIerr(ier, &
       trim(cerr)//' SFHbase%t')
  SFHbase%t(1) = SFHbase%tl
  DO i = 2, paramSFH%Nt
     SFHbase%t(i) = SFHbase%tstep * dble(i-1) + SFHbase%tl
  ENDDO
  SFHbase%t(Ntp1) = SFHbase%tu

  tmp = 0.5d0 * SFHbase%tstep_yr
  allocate(SFHbase%t_yr(Ntp1), stat=ier); call CheckIerr(ier, &
       trim(cerr)//' SFHbase%t_yr')
  SFHbase%t_yr(:) = SFHbase%t(:) * param%th_yr + tmp ! time [yr]
  allocate(SFHbase%tlk_yr(Ntp1), stat=ier); call CheckIerr(ier, &
       trim(cerr)//' SFHbase%tlk_yr')
  SFHbase%tlk_yr(:) = const%tout * param%th_yr - SFHbase%t_yr(:) ! look-back time [yr]

  paramSFH%Nzp1 = paramSFH%Nz + 1; Nzp1 = paramSFH%Nzp1
  allocate(SFHbase%Z(Nzp1), stat=ier); call CheckIerr(ier, &
       trim(cerr)//' SFHbase%Z')
  SFHbase%Z(1) = 0.d0
  ! --- metallicity bin for BC03: paramSFH%Nz = 6
  SFHbase%Z(2) = 1.d-4; SFHbase%Z(3) = 4.d-4; SFHbase%Z(4) = 4.d-3
  SFHbase%Z(5) = 8.d-3; SFHbase%Z(6) = 2.d-2; SFHbase%Z(7) = 5.d-2
  ! correction for the definition of metallicity in the Mitaka model
  tmp = const%Zsun / 2.d-2 ! 2.d-2: solar metallicity in BC03
  DO i = 2, Nzp1
     SFHbase%Z(i) = SFHbase%Z(i) * tmp
  ENDDO

  allocate(mSFH(Ntp1, Nzp1), stat=ier); call CheckIerr(ier, trim(cerr)//' mSFH')
  IF(run_redshift == 1) THEN
     allocate(CSFH(Ntp1, Nzp1, paramSFH%N_AV), stat=ier); call CheckIerr(ier, &
          trim(cerr)//' CSFH')
  ENDIF
END SUBROUTINE AllocateSFHRelated
!!$========================================================================
SUBROUTINE DeallocateSFHRelated(run_redshift)
  use global_var; use SFHrelated
  implicit none
  INTEGER, INTENT(IN) :: run_redshift
  INTEGER :: ier
  CHARACTER*50 :: cerr = '# (SFH) DeallocateSFHRelated: fail to deallocate:'

  deallocate(SFHbase%t,      stat=ier); call CheckIerr(ier, &
       trim(cerr)//' SFHbase%t')
  deallocate(SFHbase%t_yr,   stat=ier); call CheckIerr(ier, &
       trim(cerr)//' SFHbase%t_yr')
  deallocate(SFHbase%tlk_yr, stat=ier); call CheckIerr(ier, &
       trim(cerr)//' SFHbase%tlk_yr')
  deallocate(SFHbase%Z,      stat=ier); call CheckIerr(ier, &
       trim(cerr)//' SFHbase%Z')

  deallocate(mSFH, stat=ier); call CheckIerr(ier, trim(cerr)//' mSFH')
  IF(run_redshift == 1) THEN
     deallocate(CSFH, stat=ier); call CheckIerr(ier, trim(cerr)//' CSFH')
  ENDIF
END SUBROUTINE DeallocateSFHRelated
!!$========================================================================
SUBROUTINE AllocSFH(Ngal)
  use global_var; use SFHrelated
  implicit none
  INTEGER, INTENT(IN) :: Ngal
  INTEGER :: igal, ier, Ntp1, Nzp1
  CHARACTER*50 :: cerr = '# (SFH) AllocSFH: fail to allocate:'

  allocate(reff(Ngal),     stat=ier); call CheckIerr(ier, trim(cerr)//' reff')
  allocate(tauV_SFH(Ngal), stat=ier); call CheckIerr(ier, trim(cerr)//' tauV_SFH')
  reff(:) = 0.d0; tauV_SFH(:) = 0.d0

  Ntp1 = paramSFH%Ntp1; Nzp1 = paramSFH%Nzp1
  DO igal = 1, Ngal
     allocate(gal(igal)%SFH(Ntp1, Nzp1),      stat=ier); call CheckIerr(ier, &
          trim(cerr)//' gal(igal)%SFH');      gal(igal)%SFH(:,:) = 0.d0
     allocate(gal_prev(igal)%SFH(Ntp1, Nzp1), stat=ier); call CheckIerr(ier, &
          trim(cerr)//' gal_prev(igal)%SFH'); gal_prev(igal)%SFH(:,:) = 0.d0
     allocate(gal_next(igal)%SFH(Ntp1, Nzp1), stat=ier); call CheckIerr(ier, &
          trim(cerr)//' gal_next(igal)%SFH'); gal_next(igal)%SFH(:,:) = 0.d0

!!$     allocate(gal(igal)%SFH0(Ntp1, Nzp1),      stat=ier); call &
!!$          CheckIerr(ier, trim(cerr)//' gal(igal)%SFH')
!!$     allocate(gal_prev(igal)%SFH0(Ntp1, Nzp1), stat=ier); call &
!!$          CheckIerr(ier, trim(cerr)//' gal_prev(igal)%SFH')
!!$     allocate(gal_next(igal)%SFH0(Ntp1, Nzp1), stat=ier); call &
!!$          CheckIerr(ier, trim(cerr)//' gal_next(igal)%SFH')
!!$     gal(igal)%SFH0(:,:) = 0.d0; gal_prev(igal)%SFH0(:,:) = 0.d0; &
!!$          gal_next(igal)%SFH0(:,:) = 0.d0
  ENDDO
END SUBROUTINE AllocSFH
!!$========================================================================
SUBROUTINE DeallocSFH(Ngal)
  use global_var; use SFHrelated
  implicit none
  INTEGER, INTENT(IN) :: Ngal
  INTEGER :: igal, ier
  CHARACTER*50 :: cerr = '# (SFH) DeallocSFH: fail to deallocate:'

  deallocate(reff,     stat=ier); call CheckIerr(ier, trim(cerr)//' reff')
  deallocate(tauV_SFH, stat=ier); call CheckIerr(ier, trim(cerr)//' tauV_SFH')

  DO igal = 1, Ngal
     deallocate(gal(igal)%SFH,      stat=ier); call CheckIerr(ier, &
          trim(cerr)//' gal(igal)%SFH')
     deallocate(gal_prev(igal)%SFH, stat=ier); call CheckIerr(ier, &
          trim(cerr)//' gal_prev(igal)%SFH')
     deallocate(gal_next(igal)%SFH, stat=ier); call CheckIerr(ier, &
          trim(cerr)//' gal_next(igal)%SFH')
  ENDDO
END SUBROUTINE DeallocSFH
!!$========================================================================
SUBROUTINE ReplacePhysQuantsForSFH(id_new, id_old)
  use global_var; use SFHrelated
  implicit none
  INTEGER, INTENT(IN) :: id_new, id_old

  gal(id_new)%SFH(:,:) = gal(id_old)%SFH(:,:)
END SUBROUTINE ReplacePhysQuantsForSFH
!!$========================================================================
SUBROUTINE EraseGalaxyForSFH(id)
  use global_var; use SFHrelated
  implicit none
  INTEGER, INTENT(IN) :: id

  gal(id)%SFH(:,:) = 0.d0
!!$  gal(id)%SFH0(:,:) = 0.d0
END SUBROUTINE EraseGalaxyForSFH
!!$========================================================================
SUBROUTINE CalSFHrelatedInCom(sftype, id)
  use global_var; use SFHrelated
  implicit none
  INTEGER, INTENT(IN) :: sftype, id

  gal(id)%SFH(:,:) = gal(id)%SFH(:,:) + mSFH(:,:) * Mc0
END SUBROUTINE CalSFHrelatedInCom
!!$========================================================================
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


  i_ts = bin_tSFH(tstart); i_tf = bin_tSFH(tfinal)
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
        call CalcSFHInEachMetalBin(i_ts, k0, k, zc_SFH_pre, abinv)
     ELSE
        i_tsp1 = i_ts + 1; telaps = SFHbase%t(i_tsp1) - tstart
        zc_SFH_pre = Zc0; zc_SFH = CalZcSFH(zc_SFH_pre,telaps,tmp1)
        k0 = bin_metal_SFH(zc_SFH_pre); k = bin_metal_SFH(zc_SFH)
        call CalcSFHInEachMetalBin(i_ts, k0, k, zc_SFH_pre, abinv)

        telaps = SFHbase%tstep; i_tfm1 = i_tf - 1
        DO i = i_tsp1, i_tfm1
           zc_SFH_pre = zc_SFH ! metallicity at t = SFHbase%t(i)-tstart
           zc_SFH     = CalZcSFH(zc_SFH_pre, telaps, tmp1)
                        ! metallicity at t = SFHbase%t(i)-tstart+telaps
           k0 = bin_metal_SFH(zc_SFH_pre); k = bin_metal_SFH(zc_SFH)
           tmp2 = exp(-(SFHbase%t(i) - tstart) * taueffinv) * abinv
                  ! SFR at t = SFHbase%t(i)-tstart
           call CalcSFHInEachMetalBin(i, k0, k, zc_SFH_pre, tmp2)
        ENDDO

        telaps = tfinal - SFHbase%t(i_tf)
        IF(telaps > 0.d0) THEN
           zc_SFH_pre = zc_SFH; zc_SFH = CalZcSFH(zc_SFH_pre,telaps,tmp1)
           k0 = bin_metal_SFH(zc_SFH_pre); k = bin_metal_SFH(zc_SFH)
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
          IF(ISNAN(mSFH(i, j)) .eqv. .true.) THEN
             mSFH(i, j) = 0.d0; print *, 'NaN'
          ENDIF
       ENDDO
    ENDDO
  END SUBROUTINE CheckmSFH
END SUBROUTINE CalcSFH
!!$========================================================================
SUBROUTINE CheckSFH(number, id)
  use global_var; use SFHrelated
  INTEGER, INTENT(IN) :: number, id
  INTEGER :: itime, iZ

  DO itime = 1, paramSFH%Ntp1
     DO iZ = 1, paramSFH%Nzp1
        IF(ISNAN(gal(id)%SFH(itime,iZ)) .eqv. .true.) THEN
           IF(number == 7) THEN
              print '(I1, A, I8, 2(A, I2), 2A, 3G10.3)', number, &
                   ': gal(', id, ')%SFH(', itime, ',', iZ, ')=NaN, ', &
                   '(Mc0,mSFH,Mcool)=', Mc0, mSFH(itime,iZ), gal(id)%Mcoold
           ELSE
              print '(I1, A, I8, 2(A, I2), A)', number, &
                   ': gal(', id, ')%SFH(', itime, ',', iZ, ')=NaN'
           ENDIF
        ENDIF
     ENDDO
  ENDDO
END SUBROUTINE CheckSFH
!!$========================================================================
SUBROUTINE SubstReffTauV(id, mor, rb, rd, x)
  use SFHrelated
  implicit none
  INTEGER, INTENT(IN) :: id, mor
  DOUBLE PRECISION, INTENT(IN) :: rd, rb
  DOUBLE PRECISION, INTENT(IN) :: x ! = param%tauV0 * (MZc(id)/const%Zsun)
  DOUBLE PRECISION :: Square ! function

  IF(mor == 3) THEN ! S
     reff(id) = rd
  ELSE ! E or S0
     reff(id) = rb
  ENDIF

  IF(reff(id) > 0.d0 .and. x > 0.d0) THEN
     tauV_SFH(id) = x / Square(reff(id))
  ELSE
     tauV_SFH(id) = 0.d0
  ENDIF
END SUBROUTINE SubstReffTauV
!!$========================================================================
SUBROUTINE WriteSFH(igal, nz, iforest, inv_V, rb0, rd0, mor, mor_d, Mssat)
  use global_var; use SFHrelated
  implicit none
  INTEGER, INTENT(IN) :: igal, nz, iforest, mor, mor_d
  DOUBLE PRECISION, INTENT(IN) :: inv_V, rb0, rd0, Mssat
  INTEGER :: icall = 0, itime, iZ, iband
  DOUBLE PRECISION :: SFRtot, Zmassb, Zmassd, SFRtot0


  ! print '(A, I10)', '# calling WriteSFH: icall = ', icall
  ! --- writing physical quantities of each galaxy
  write(i_SFH(2), '(2I10, 4I2, I8, 4I16, 10G13.5, $)') &
       igal, gal(igal)%id_cgal, mor, mor_d, gal(igal)%flag_c, & ! 1,2,3,4,5
       gal(igal)%flag_burst, iforest, gal(igal)%IDhost, & ! 6,7,8
       gal(igal)%mpi, mrgt%mpi(gal(igal)%IDhost), & ! 9,10
       mrgt%f_des(gal(igal)%IDhost), & ! 11
       gal(igal)%Mstarb * param%munit, & ! 12
       gal(igal)%Mstard * param%munit, & ! 13
       gal(igal)%Mcoold  * param%munit, & ! 14
       mrgt%mhalo(gal(igal)%IDhost) * param%munit, & ! 15
       gal(igal)%Mhalo * param%munit,            & ! 16
       gal(gal(igal)%id_cgal)%Vcent, gal(igal)%Vcent, & ! 17,18
       gal(igal)%SFR, gal(igal)%mSFR, & ! 19,20
       gal(igal)%Mbh * param%munit ! 21
  Zmassb = 0.d0; Zmassd = 0.d0
  IF(gal(igal)%Mstarb > 0.d0) Zmassb = gal(igal)%Zmassb / gal(igal)%Mstarb
  IF(gal(igal)%Mstard > 0.d0) Zmassd = gal(igal)%Zmassd / gal(igal)%Mstard
  write(i_SFH(2), '(3G13.5, $)') &
       Zmassb, Zmassd, (gal(igal)%Zmassb + gal(igal)%Zmassd) / Mssat
  DO iband = 2, param%nwave
     write(i_SFH(2), '(2F11.4, $)') mag(iband),  mag_d(iband)
  ENDDO
  write(i_SFH(2), *)

  ! --- writing SFH of each galaxy
!!$  SFRtot = 0.d0
!!$  DO itime = 1, paramSFH%Nt-1
!!$     DO iZ = 1, paramSFH%Nzp1
!!$        SFRtot = SFRtot + gal(igal)%SFH(itime, iZ)
!!$     ENDDO
!!$  ENDDO
!!$  write(i_SFH(4), '(2(A, I10), A, I16, A, G12.4, A, G12.4)') &
!!$       '# index ', icall, ': igal = ', igal, ', mpi = ', gal(igal)%mpi, &
!!$       ', Mstar = ', Mssat * param%munit, ' [Msun]: ', &
!!$       SFRtot * param%munit * ssp%alp(1)
  write(i_SFH(4), '(2(A, I10), A, I16, A, G12.4, A)') &
       '# index ', icall, ': igal = ', igal, ', mpi = ', gal(igal)%mpi, &
       ', Mstar = ', Mssat * param%munit, ' [Msun]'
  DO itime = 1, paramSFH%Nt-1
!!$     write(i_SFH(4), '(2G13.5, $)') &
!!$       SFHbase%t_yr(itime), SFHbase%tlk_yr(itime)
!!$     DO iZ = 1, paramSFH%Nzp1
!!$        IF(gal(igal)%SFH(itime, iZ) * paramSFH%SFRunit < 1.d-30) &
!!$             gal(igal)%SFH(itime, iZ) = 0.d0
!!$        write(i_SFH(4), '(G13.5, $)') &
!!$          gal(igal)%SFH(itime, iZ) * paramSFH%SFRunit
!!$     ENDDO
!!$     write(i_SFH(4), *)
     SFRtot = 0.d0
!!$     SFRtot0 = 0.d0
     DO iZ = 1, paramSFH%Nzp1
        SFRtot = SFRtot + gal(igal)%SFH(itime, iZ)
!!$        SFRtot0 = SFRtot0 + gal(igal)%SFH0(itime, iZ)
     ENDDO
     IF(SFRtot * paramSFH%SFRunit < 1.d-30) SFRtot = 0.d0
!!$     IF(SFRtot0 * paramSFH%SFRunit < 1.d-30) SFRtot0 = 0.d0
!!$     write(i_SFH(4), '(2G13.5)') SFRtot * paramSFH%SFRunit, &
!!$          (SFRtot-SFRtot0) * paramSFH%SFRunit
     write(i_SFH(4), '(G13.5)') SFRtot * paramSFH%SFRunit
  ENDDO
  write(i_SFH(4), *); write(i_SFH(4), *)
  icall = icall + 1
END SUBROUTINE WriteSFH
!!$========================================================================
SUBROUTINE CalcCSFH(igal, flag_burst, weight, Mhalo)
  use global_var; use SFHrelated
  implicit none
  INTEGER, INTENT(IN) :: igal, flag_burst
  DOUBLE PRECISION, INTENT(IN) :: weight, Mhalo
  INTEGER :: itime, iZ, bin, b_or_q ! 1:starburst, 2:quiescent
  DOUBLE PRECISION :: Mstar, Av


  Mstar = gal(igal)%Mstarb + gal(igal)%Mstard; b_or_q = 2
  IF(flag_burst == 1) b_or_q = 1
  IF(Mstar > 0.d0) THEN
     Mstar = Mstar * param%munit; Av = CalAv(tauV_SFH(igal))
     bin   = int(aint(Av/paramSFH%stepAV)) + 1
     ! paramSFH%stepAV*a <= Av < paramSFH%stepAV*(a+1) --> bin = a+1
     ! the center of bin=a+1: Av = paramSFH%stepAV*(a+0.5)
     IF(bin > paramSFH%N_AV) bin = paramSFH%N_AV

     DO itime = 1, paramSFH%Nt
        DO iZ = 1, paramSFH%Nzp1
           call FinalCheckSFH(itime, iZ, igal, gal(igal)%SFH(itime,iZ)*weight)
           CSFH(itime,iZ, bin) = CSFH(itime,iZ, bin) &
                + gal(igal)%SFH(itime,iZ) * weight
           call CheckCSFH(itime, iZ, igal, bin)
        ENDDO
     ENDDO

     call CheckConsistency
  ENDIF
CONTAINS
  DOUBLE PRECISION FUNCTION CalAv(tauV) RESULT(Av)
    DOUBLE PRECISION, INTENT(IN) :: tauV
    DOUBLE PRECISION :: f_slab ! function

    IF(tauV > 0.d0) THEN
       Av = -2.5d0 * log10(f_slab(tauV))
    ELSE
       Av = 0.d0
    ENDIF
  END FUNCTION CalAv

  SUBROUTINE FinalCheckSFH(itime, iZ, igal, x)
    INTEGER, INTENT(IN) :: itime, iZ, igal
    DOUBLE PRECISION, INTENT(IN) :: x ! = gal(igal)%SFH(itime,iZ)*weight

    IF(ISNAN(x) .eqv. .true.) THEN
       gal(igal)%SFH(itime,iZ) = 0.d0
       print '(A, I8, 2(A, I2), A)', 'NaN: gal(', igal, ')%SFH(',itime,',',iZ,')'
    ENDIF
  END SUBROUTINE FinalCheckSFH

  SUBROUTINE CheckCSFH(itime, iZ, igal, bin)
    INTEGER, INTENT(IN) :: itime, iZ, igal, bin

    IF(ISNAN(CSFH(itime,iZ,bin)) .eqv. .true.) THEN
       print '(2(A, I2), A, I5, A, I8, 2(A, I2), A, G10.3)', &
            'CSFH(',itime,',',iZ,',',bin,')=NaN: gal(', igal, &
            ')%SFH(',itime,',',iZ, ')*weight=', &
            gal(igal)%SFH(itime,iZ)*weight
    ENDIF
  END SUBROUTINE CheckCSFH

  SUBROUTINE CheckConsistency
    DOUBLE PRECISION :: tmp

    tmp = 0.d0
    DO itime = 1, paramSFH%Nt-1
       DO iZ = 1, paramSFH%Nzp1
          tmp = tmp + gal(igal)%SFH(itime, iZ)
       ENDDO
    ENDDO
    tmp = (ssp%alp(b_or_q) * tmp * param%munit - Mstar) / Mstar
    IF(abs(tmp) > 1.d0 .and. Mstar >= 1.d+5) THEN
       write(i_SFH(4), '(2I8, 7G10.2)') paramSFH%loop, igal, &
            Mhalo*param%munit, Mstar, tmp, Av, reff(igal), &
            tauV_SFH(igal), weight

       write(i_SFH(5), '(A, I8)') '# index ', paramSFH%loop
       DO itime = 1, paramSFH%Nt-1
          write(i_SFH(5), '(2G10.2, $)') &
               SFHbase%t_yr(itime), SFHbase%tlk_yr(itime)
          DO iZ = 1, paramSFH%Nzp1
             IF(iZ == paramSFH%Nzp1) THEN
                write(i_SFH(5), '(G10.2)') &
                     gal(igal)%SFH(itime, iZ) * paramSFH%SFRunit
             ELSE
                write(i_SFH(5), '(G10.2, $)') &
                     gal(igal)%SFH(itime, iZ) * paramSFH%SFRunit
             ENDIF
          ENDDO
       ENDDO
       write(i_SFH(5), *); write(i_SFH(5), *)

       paramSFH%loop = paramSFH%loop + 1
    ENDIF
  END SUBROUTINE CheckConsistency
END SUBROUTINE CalcCSFH
!!$========================================================================
SUBROUTINE WriteCSFH
  use global_var; use SFHrelated
  implicit none
  INTEGER :: indx, iZ, itime, iAV
  DOUBLE PRECISION :: tmp1, tmp2
  DOUBLE PRECISION :: Cube

  print '(A)', '# (SFH) Writing the CSFH data into output file'

  indx = 0
  tmp1 = Cube(param%h) * paramSFH%SFRunit
  DO iZ = 1, paramSFH%Nzp1
     write(i_SFH(1), '(A, I2, $)') '# index ', indx
     IF(iZ == paramSFH%Nzp1) THEN
        write(i_SFH(1), '(A, G10.2)') ': Z > ', SFHbase%Z(iZ)
     ELSE
        write(i_SFH(1), '(2(A, G10.2))') &
             ': Z = ', SFHbase%Z(iZ), '--', SFHbase%Z(iZ+1)
     ENDIF

     DO itime = 1, paramSFH%Nt-1
        write(i_SFH(1), '(2G12.4, $)') SFHbase%t_yr(itime), SFHbase%tlk_yr(itime)
        DO iAV = 1, paramSFH%N_AV
           tmp2 = CSFH(itime, iZ, iAV) * tmp1
           IF(ISNAN(tmp2) .eqv. .true.) THEN
              print '(6(A, I3), 2(A, G10.3))', 'CSFH(',itime,',',iZ,',',iAV, &
                   ')*tmp1=NaN: CSFH(',itime,',',iZ,',',iAV,')=', &
                   CSFH(itime,iZ,iAV), ', tmp1=', tmp1
              write(i_SFH(1), '(G14.6, $)') 0.d0
           ELSEIF(tmp2 < 1.d-90) THEN
              write(i_SFH(1), '(G14.6, $)') 0.d0
           ELSE
              write(i_SFH(1), '(G14.6, $)') tmp2
           ENDIF
        ENDDO
        write(i_SFH(1), *)
     ENDDO
     write(i_SFH(1), *); write(i_SFH(1), *)
     indx = indx + 1
  ENDDO
END SUBROUTINE WriteCSFH
!!$========================================================================
SUBROUTINE CloseFilesForSFH
  use SFHrelated; implicit none
  INTEGER :: i, ii, ier
  CHARACTER(LEN=50) :: cerr = '# (SFH) CloseFilesForSFH: fail to'

  deallocate(i_SFH, stat=ier); call CheckIerr(ier, trim(cerr)//' deallocate: i_SFH')
  DO i = 1, paramSFH%nfile
     ii = i + paramSFH%base
     close(ii, iostat=ier); call CheckIerr(ier, trim(cerr)//' close file:'//&
          trim(fnameSFH(i)))
  ENDDO
  deallocate(fnameSFH, stat=ier); call CheckIerr(ier, trim(cerr)//' deallocate: fnameSFH')
END SUBROUTINE CloseFilesForSFH
!!$========================================================================
