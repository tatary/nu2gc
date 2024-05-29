PROGRAM main_nugc
  use global_var; use LAErelated; use SFHrelated; use SNrelated
  use CLCOOLrelated 
  implicit none

  include 'mpif.h' ! for parallel

!!$  CHARACTER(LEN=50) :: fbase = '/home/makiya/nugc/data/', fbase_nbody = '/work/'
!!$  CHARACTER(LEN=50) :: fbase = 'data/', fbase_nbody = '/work/'
  CHARACTER(LEN=50) :: fbase = 'data/', fbase_nbody

  INTEGER :: i, j, k, ier
  INTEGER :: i_ssp, iforest, me, i_file, iout, iout_q, igal, iband, ihost
  INTEGER :: nop ! No. of parameters for MCMC
  DOUBLE PRECISION :: x, inv_V
  CHARACTER(LEN=50) :: cerr_a = '# main_nugc: failed to allocate: '
  CHARACTER(LEN=50) :: cerr_d = '# main_nugc: failed to deallocate: '

  DOUBLE PRECISION, PARAMETER :: EPS = 1.d-12
  CHARACTER(LEN=100) :: c_mftype(11) = (/ CHARACTER(LEN=100) ::&
         'bulge stellar mass', 'disk stellar mass', 'cold gas mass', &
         'SMBH mass', 'hot gas mass', 'bulge mass (stellar + SMBH)', &
         'disk mass (stellar + cold gas)', &
         'galaxy mass (stellar + SMBH + cold gas)', &
         'stellar mass (bulge + disk)', 'hydrogen gas mass', 'halo mass'/)
  TYPE DistributionFunction2
     INTEGER :: Nbin, N1, N2
     DOUBLE PRECISION :: step, invstep, base
     DOUBLE PRECISION, ALLOCATABLE :: bin(:)
     DOUBLE PRECISION, ALLOCATABLE :: x(:,:,:), n(:,:,:), xn(:,:,:), xxn(:,:,:)
  END type DistributionFunction2
  TYPE(DistributionFunction2), ALLOCATABLE :: ml(:)
                                              ! M_HI/L_B as a function of M_B
  TYPE(DistributionFunction2), ALLOCATABLE :: MbhMbulge
                                              ! M_BH as a function of Mbulge
  TYPE(DistributionFunction2), ALLOCATABLE :: DiskScale
                                              ! R_d -- V_d relation of Spirals
  TYPE(DistributionFunction2), ALLOCATABLE :: SphScale
                                              ! fundamental plane for spheroids 

  INTEGER, PARAMETER :: ionum = 30
  CHARACTER(LEN=100) :: file_out(0:NFILE)
  CHARACTER(LEN=50)  :: file_o, file_i
  DOUBLE PRECISION, ALLOCATABLE :: lumtmp(:), lumtmp_d(:)

  ! CHARACTER(LEN=50) :: filessp(2)
  DOUBLE PRECISION :: rb0, rd0, zg, zt
  DOUBLE PRECISION :: tots, totc, toth, totd, totlum

  DOUBLE PRECISION :: MZctmp, temp

  CHARACTER(LEN=20), PARAMETER :: base_ssp(N_SSP) &
       = (/CHARACTER(LEN=20)::'_KA97', '_PEGASE', '_BC03'/)
  CHARACTER(LEN=200) :: header
  DOUBLE PRECISION :: Vdomi, rdomi, Vmino, rmino
  DOUBLE PRECISION :: Zmassb, Zmassd

  INTEGER :: idum_org
  CHARACTER(LEN=20) :: c_sftype(2) = (/'# starburst --- ', '# quiescent --- '/)

  DOUBLE PRECISION ::fobs 
  ! --- functions
  DOUBLE PRECISION :: Cube

  ! --- surface brightness related, added by Makiya 2014/07/09
  DOUBLE PRECISION :: muSB

  ! for CSFR
  CHARACTER(LEN=100) :: CSFR_out

  call MPI_INIT(istat)
  call MPI_Comm_size(MPI_COMM_WORLD, nnode, istat)
  call MPI_Comm_rank(MPI_COMM_WORLD, inode, istat)
  NgalMax = int(NgalMax/nnode); NgalMax2 = int(NgalMax2/nnode)


  ! --- Initialize random functions. added  by Makiya
  param%iy = 0; param%iv(:) = 0; param%iset = 0

  IF(inode == 0) THEN
     print *
     print '(A)', '#####################################################'
     print '(A)', '#        Execution of the nuGC model starts         #'
     print '(A)', '#####################################################'
  ENDIF

  tots = 0.d0; totc = 0.d0; toth = 0.d0; totd = 0.d0

!!$       === Set Input File Name & Read the Input Parameters ===
  call getarg(1, file_i)
  file_o = file_i(2:len_trim(file_i))
  call ReadParameters
  fbase_nbody = '../'
  IF(param%file_nbody(1:4) == '1120') fbase_nbody = '/work2/'
  call SubstFilterNumb
  idum_org = param%idum ! the original value of idum is memorized

!!$       === Identify the N-body Run and Calculate the N-body Related Quantities ===
  i_file = 1
  open(i_file, file = trim(fbase_nbody)//trim(param%file_nbody)//'/param', &
       status = 'old', iostat=ier)
  call CheckIerr(ier, 'main_nugc: fail to open file: '//trim(fbase_nbody)//&
       trim(param%file_nbody)//'/param')
  read(i_file, *) param%n_forest
  read(i_file, *) x ! simulation box size [Mpc/h]
  read(i_file, *) param%OMEGA0
  read(i_file, *) param%h
  read(i_file, *) param%OMEGA
  close(i_file)

  IF(param%n_forest < nnode) THEN
      print *, "number of node should be less than number of forest"
      stop
  ENDIF

  inv_V = 1.d0 / Cube(x) ! [h^3/Mpc^3]
  param%OMEGA_L   = 1.d0 - param%OMEGA0 ! Omega_L
  param%OMEGA_rat = param%OMEGA_L / param%OMEGA0 ! Omega_L / Omega_M
  param%bar_rat   = param%OMEGA / param%OMEGA0 ! = Omega_b / Omega_M

  param%th     = 1.d0 / (0.1024d0 * param%h) ! hubble time [Gyr] (not in Gyr/h)
  param%th_yr  = param%th * 1.d+9 ! hubble time [yr]
  param%tau0st = param%tau0st / param%th ! [Gyr] --> [hubble time]
  call SetConstants
  IF(param%SB) call ReadPetSBcor ! reading the relevant data for calculation of LFobs
                                 !   w/ a given limited surface brightness

  ! --- concentration parameter, added by Makiya (2015/10/07)
  call ReadCST

  call CreateOutputFileNames(param%run_all_ssp, param%run_type)

!!$      Reading cooling data (revised by Shirakata 18/Aug/17)
  IF(param%CoolFN == 1) THEN  ! Sutherland & Dopita 93
     call ReadCoolingData
  ELSE IF(param%CoolFN == 2) THEN ! Okamoto-cloudy
     call ReadCLCoolingData
  ENDIF


!!$       === LAE Related Preparation ===
  IF(param%LAE) THEN
     call SetAndOpenOutputFilesForLAE(fbase, file_o, param%run_type, nnode, inode)
     call AllocSchaererArrays
     call ReadSchaererSSP
     call CorrectionSchaererSSP(inode)
          ! correction to the different cutoff mass used in the S03 SSPs
     call InitializeArraysForLAE
  ENDIF

!!$       === SFH Related Preparation ===
  IF(param%SFH) call SetAndOpenOutputFilesForSFH(fbase, file_o, param%run_type)

!!$       === SN Related Preparation ===
  IF(param%SN) call SetAndOpenOutputFilesForSN(fbase, file_o, param%run_type)

  IF(param%run_type == 1) THEN ! run at a redshift
     IF(inode == 0) THEN
        print '(A)', '# Calculation for a redshift is executed'
        print '(A, $)', '# -------------------------------------------------'
        print '(A)', '--------------------------'
     ENDIF
     IF(param%run_all_ssp) THEN
!!$             open(,file = 'list_sspfiles_AB', status = 'old')
!!$             read()
        DO i_ssp = 1, N_SSP_TOT
!!$                 read() ssp%filename(1)
!!$                 ssp%filename(2) = ssp%filename(1)
           call CalcGalaxiesAtARedshift(param%zsp1)
        ENDDO
!!$             close()
     ELSE
        call CalcGalaxiesAtARedshift(param%zsp1)
     ENDIF
  ELSEIF(param%run_type == 2) THEN ! run at all redshift
     IF(inode == 0) THEN
        print '(A)', '# Calculation for all redshift is executed'
        print '(A, $)', '# -------------------------------------------------'
        print '(A)', '--------------------------'
     ENDIF
     call CalcGalaxiesAtAllRedshift
  ELSEIF(param%run_type == 3) THEN ! MCMC
     nop = 0
     call ReadMCMCParameters(nop)
     IF(inode == 0) THEN
        print '(A)', '# MCMC Calculation is executed'
        print '(A, $)', '# -------------------------------------------------'
        print '(A)', '--------------------------'
     ENDIF
     call CalcMCMC(nop)
  ENDIF

  IF(param%SB) THEN
     deallocate(SB%Pet2TotRatio, stat=ier); call CheckIerr(ier, &
          trim(cerr_d)//' SB%Pet2TotRatio')
  ENDIF
  deallocate(param%iwave, stat=ier); call CheckIerr(ier, &
       trim(cerr_d)//' param%iwave')
  deallocate(mag,      stat=ier); call CheckIerr(ier, trim(cerr_d)//' mag')
  deallocate(mag_d,    stat=ier); call CheckIerr(ier, trim(cerr_d)//' mag_d')
  deallocate(xi,       stat=ier); call CheckIerr(ier, trim(cerr_d)//' xi')
  deallocate(lumtmp,   stat=ier); call CheckIerr(ier, trim(cerr_d)//' lumtmp')
  deallocate(lumtmp_d, stat=ier); call CheckIerr(ier, trim(cerr_d)//' lumtmp_d')
  deallocate(mag_q,    stat=ier); call CheckIerr(ier, trim(cerr_d)//' mag_q')

  IF(inode == 0) THEN
     print *
     print '(A)', '#####################################################'
     print '(A)', '# Execution of the nuGC model is done successfully  #'
     print '(A)', '#####################################################'
     print *
  ENDIF

  call MPI_FINALIZE(istat)
!!$============================================================================
CONTAINS
!!$============================================================================
  SUBROUTINE CalcGalaxiesAtARedshift(zsp1)
    DOUBLE PRECISION, INTENT(INOUT) :: zsp1 ! param%zsp1
    INTEGER :: endhalo, end_step, nz, nz_end, i_file, ier
    INTEGER :: i, bin, i25, i50, i75
    INTEGER :: mor, mor_d, b_or_d
    DOUBLE PRECISION :: Mssat ! = Mstar^bulge + Mstar^disk
!!$    DOUBLE PRECISION :: bt, bt_d ! B/T ratio in B-band w/o and w/ dust
    DOUBLE PRECISION :: galmass, Mdisk, Mbulge, M_H ! for mass functions
    DOUBLE PRECISION :: MassLumi ! for M_HI/L_B distribution
!!$    DOUBLE PRECISION :: log10ML  ! for M_HI/L_B distribution
    DOUBLE PRECISION, ALLOCATABLE :: arr(:) ! for M_HI/L_B distribution
    ! --- functions
    INTEGER :: DetMorType
    DOUBLE PRECISION :: Square, AbsDiff, z2t, CalBulgeReff, CalDiskReff
    DOUBLE PRECISION :: ObsFracAGN_UV

    CHARACTER(LEN=50) :: cerr_o = '# CalcGalaxiesAtARedshift: fail to open file '
    CHARACTER(LEN=50) :: cerr_a = '# CalcGalaxiesAtARedshift: fail to allocate '
    CHARACTER(LEN=50) :: cerr_d = '# CalcGalaxiesAtARedshift: fail to deallocate '

    DO i_file = 0, NFILE
       i = i_file + ionum
       call CheckIONum(i, 'CalcGalaxiesAtARedshift')
       open(i, file=file_out(i_file), status='unknown', iostat=ier)
       call CheckIerr(ier, trim(cerr_o)//' '//trim(file_out(i_file)))
    ENDDO
    nz_end = 1; nz = 1

    !-- for parallel run, added by Makiya (2015/09/27)
    ! galaxy catalog
    i = ionum + NFILE + inode + 1
    call CheckIONum(i, 'CalcGalaxiesAtARedshift')
    open(i, file=file_catalog, status='unknown', iostat=ier)
    call CheckIerr(ier, trim(cerr_o)//' '//trim(file_catalog))
    ! quasar catalog
    i = ionum + NFILE + nnode + inode + 1
    call CheckIONum(i, 'CalcGalaxiesAtARedshift')
    open(i, file=file_catalog_q, status='unknown', iostat=ier)
    call CheckIerr(ier, trim(cerr_o)//' '//trim(file_catalog_q))

    IF(param%traceIDs == 1) THEN
       ! target catalog
       i = ionum + NFILE + nnode + nnode + inode + 1
       call CheckIONum(i, 'CalcGalaxiesAtARedshift')
       open(i, file=file_catalog_t, status='unknown', iostat=ier)
       call CheckIerr(ier, trim(cerr_o)//' '//trim(file_catalog_t))
    ENDIF

    N1 = param%nwave; N2 = 4
    N1N2Nbin = N1 * N2 * NbinLF
    allocate(lf_n_all(N1, N2, NbinLF), stat=ier); call CheckIerr(ier, &
         trim(cerr_a)//' lf_n_all'); lf_n_all(:,:,:) = 0.d0
    allocate(mf_n_all(NtypeMF, N2, NbinMF), stat=ier); call CheckIerr(ier, &
         trim(cerr_a)//' mf_n_all'); mf_n_all(:,:,:) = 0.d0
    allocate(lf_q_all(param%nwaveq, 1, NbinLF), stat=ier); call CheckIerr(ier, &
         trim(cerr_a)//' lf_q_all'); lf_q_all(:,:,:) = 0.d0
    IF(param%Mbh) THEN
       allocate(MbhMbulge_bin(50), stat=ier); call CheckIerr(ier, &
            trim(cerr_a)//' MbhMbulge_bin'); MbhMbulge_bin(:) = 0.d0
       allocate(MbhMbulge_n_all(1,   N2, 50), stat=ier); call CheckIerr(ier, &
            trim(cerr_a)//' MbhMbulge_n_all');   MbhMbulge_n_all(:,:,:)   = 0.d0
       allocate(MbhMbulge_xn_all(1,  N2, 50), stat=ier); call CheckIerr(ier, &
            trim(cerr_a)//' MbhMbulge_xn_all');  MbhMbulge_xn_all(:,:,:)  = 0.d0
       allocate(MbhMbulge_xxn_all(1, N2, 50), stat=ier); call CheckIerr(ier, &
            trim(cerr_a)//' MbhMbulge_xxn_all'); MbhMbulge_xxn_all(:,:,:) = 0.d0
    ENDIF
    IF(param%ML) THEN
       allocate(ML_bin(50), stat=ier); call CheckIerr(ier, &
            trim(cerr_a)//' ML_bin'); ML_bin(:) = 0.d0
       allocate(ML_n_all(2, 4, 50), stat=ier); call CheckIerr(ier, &
            trim(cerr_a)//' ML_n_all'); ML_n_all(:,:,:) = 0.d0
       IF(inode == 0) THEN
          allocate(ML_x_all(NgalMax*nnode, 4, 50), stat=ier)
       ELSE
          allocate(ML_x_all(NgalMax, 4, 50), stat=ier)
       ENDIF
       call CheckIerr(ier, trim(cerr_a)//' ML_x_all'); ML_x_all(:,:,:) = 0.d0
       allocate(ML_med(3, 4, 50), stat=ier); call CheckIerr(ier, &
            trim(cerr_a)//' ML_med'); ML_med(:,:,:) = 0.d0
    ENDIF


    ! write the captions of each data file
    IF(inode == 0) THEN
       call WriteCaptions(ionum, nz, nz_end)
       call WriteInputParameters(1) ! write the input parameters
                                    !  into the log file of file_out(0)
    ENDIF

    ! obtain the parameters related to merger tree
    call ReadMRGB(0,nz); call AllocateGalaxy

    ! --- CSFR, by Makiya
    allocate(CSFR(2,     mrgp%num_step), stat=ier); call CheckIerr(ier, &
         trim(cerr_a)//' CSFR');     CSFR(:,:)     = 0.d0
    allocate(CSFR_all(2, mrgp%num_step), stat=ier); call CheckIerr(ier, &
         trim(cerr_a)//' CSFR_all'); CSFR_all(:,:) = 0.d0
    CSFR(1, :)    = mrgp%zp1ar(:) ! redshift bin
    CSFR_all(1,:) = mrgp%zp1ar(:) ! redshift bin

    ! --- SMF at all z, by Makiya
    allocate(SMF_z(mrgp%num_step,     NbinMF), stat=ier); call CheckIerr(ier, &
         trim(cerr_a)//' SMF_z');     SMF_z(:,:)     = 0.d0
    allocate(SMF_z_all(mrgp%num_step, NbinMF), stat=ier); call CheckIerr(ier, &
         trim(cerr_a)//' SMF_z_all'); SMF_z_all(:,:) = 0.d0

    ! --- the calculation is stop if the requested redshift is larger than
    !       the redshift of the 2nd N-body slice
    IF(zsp1 > mrgp%zp1ar(2)) THEN
       print *
       print '(A, $)', '# ---------------------------------------------------'
       print '(A)', '---------------------------------------'
       print '(2(A, F8.5), A)', '# requested redshift (=', zsp1-1.d0, &
            ') is larger than the redshift of the 2nd N-body slice (=',&
            mrgp%zp1ar(2)-1.d0, ')'
       print '(A, F8.5, A)', '# Rerun with different redshift, which is '//&
            'smaller than ', mrgp%zp1ar(2)-1.d0, ', or change N-body data'
       print '(A, $)', '# ---------------------------------------------------'
       print '(A)', '---------------------------------------'
       print *
       stop
    ENDIF

    ! search the closest redshift slice in N-body run to param%zsp1
    ! --- zsp1: the closest redshift slice in N-bory run
    ! --- param%loopz: the correspondent integer to zsp1
    !                  = mrgp%zp1ar(mrgp%num_step-(param%loopz-1))
    !                  = mrgp%zp1ar(param%izout)
    param%loopz = 1
    DO WHILE(mrgp%zp1ar(mrgp%num_step - (param%loopz - 1)) < zsp1)
       param%loopz = param%loopz + 1
    ENDDO
    i = mrgp%num_step - (param%loopz - 1)
    IF(inode == 0) THEN
       print '(2(A,I3), A, I3, 2(A,F8.5), $)', &
            '# param%loopz = ', param%loopz, ', i = ', i, &
            ': mrgp%zp1ar(', i,   ') = ', mrgp%zp1ar(i), ', zsp1 = ', zsp1
       IF(i > 1) print '(A, I3, A, F8.5, $)', &
            ', mrgp%zp1ar(', i-1, ') = ', mrgp%zp1ar(i-1)
       print *
    ENDIF
    IF(param%loopz > 1 .and. i < mrgp%num_step) THEN
       IF(mrgp%zp1ar(i) - zsp1 > zsp1 - mrgp%zp1ar(i+1)) THEN
          param%loopz = param%loopz - 1
          IF(inode == 0) print '(2(A, I3))', '# param%loopz is modified: ', &
               param%loopz+1, ' --> ', param%loopz
       ENDIF
    ENDIF
    param%izout = mrgp%num_step - (param%loopz - 1)
    zsp1 = mrgp%zp1ar(param%izout)
    ! IF(AbsDiff(zsp1_orig, zsp1) > 1.d-20) THEN
    IF(AbsDiff(param%zsp1_input, zsp1) > 1.d-20) THEN
       ! modified by MARK (2016/06/28)
       IF(inode == 0) THEN
          print '(2(A,F8.5), A,I3, $)', &
               ! '# param%zsp1 is modified: ', zsp1_orig, ' --> ', zsp1, &
               '# param%zsp1 is modified: ', param%zsp1_input, ' --> ', zsp1, &
                                             ! modified by MARK (2016/06/28)
               ' (param%izout = ', param%izout
          IF(param%izout < mrgp%num_step) THEN
             print '(A, I3, A, F8.5, A)', &
                  ', mrgp%zp1ar(', param%izout+1, ') = ', &
                  mrgp%zp1ar(param%izout+1), ')'
          ELSE
             print '(A)', ')'
          ENDIF
       ENDIF
       param%zsp1 = zsp1 ! added by MARK (2016/Jun/24)
       ! recalculate the age of the universe at output redshift
       const%tout = z2t(param%zsp1); const%toutGyr = const%tout * param%th
    ENDIF

    call AllocateDistributionFunctions(nz_end)

    IF(param%LAE) THEN ! allocate LAE related arrays
       call AllocLAE(mrgt%num_galarray, param%nwave)
       call AllocateDistributionFunctions_LAE(nz_end, param%nwave)
    ENDIF
    IF(param%SFH) THEN ! allocate SFH related arrays
       call AllocateSFHRelated(param%run_type)
       call AllocSFH(mrgt%num_galarray)
    ENDIF

    IF(param%SN) call AllocateSNRelated(param%run_type)

    call ReadSSPFile(param%loopz, nz_end, nz)
    IF(inode == 0) call WriteInputParameters(2) ! write the input parameters
                                                !  into the log file of file_out(0)
    call optdepth(zsp1, ionum, 0, nz)
    IF(param%LAE) THEN
       IF(inode == 0) call WriteCaptionsForLAE(ionum, param%run_type)
       call optdepthForLAE(ionum)
    ENDIF
    IF(param%SFH) call WriteCaptionsForSFH(ionum)
    IF(param%SN) call WriteCaptionsForSN(ionum)

    nforest_this = param%n_forest / nnode
    IF(mod(param%n_forest, nnode) /= 0) THEN
        nforest_this = nforest_this+1
    ENDIF

    iforest_start = inode * nforest_this

    IF(param%LAE) call DeallocLAE
    IF(param%SFH) call DeallocSFH(mrgt%num_galarray)
    call DeallocateGalaxy; call DeallocMRGB

    DO iforest = iforest_start, min(iforest_start+nforest_this-1, param%n_forest-1)

       ! initialize ran1 and gasdev, for check (Makiya, 2015/09/27)
       param%idum = -131; param%iy = 0; param%iv(:) = 0; param%iset = 0

       call ReadMRGB(iforest,nz); call AllocateGalaxy
       IF(param%LAE) call AllocLAE(mrgt%num_galarray, param%nwave)
       IF(param%SFH) call AllocSFH(mrgt%num_galarray)

       ! calculate star formation in galaxies
       !  and return 'endhalo' and 'end_step' as the total number of
       !  galaxies at param%zsp1 and the corresponding number of param%zsp1
       !  for mrgp
       call star(iforest, endhalo, end_step, ionum, nz)
       gal(1:endhalo) = gal_next(1:endhalo)
!!$              call CheckNumGtot(endhalo, end_step)
       IF(param%LAE) allgal(1:endhalo) = allgal_next(1:endhalo)

       include "out_norm_l9a.f90"

       IF(param%LAE) call DeallocLAE
       IF(param%SFH) call DeallocSFH(mrgt%num_galarray)
       call DeallocateGalaxy; call DeallocMRGB
    ENDDO

    call flush(ionum+NFILE+inode+1); call flush(ionum+NFILE+nnode+inode+1)
    IF(param%traceIDs == 1) call flush(ionum+NFILE+nnode+nnode+inode+1)
    call MPI_BARRIER(MPI_COMM_WORLD, ier)
    call GatherResults(nz)

    IF(param%LAE) THEN
       call flush(paramLAE%base + paramLAE%nfile + inode + 1)
       call flush(paramLAE%base + paramLAE%nfile + nnode + inode + 1)
       call MPI_BARRIER(MPI_COMM_WORLD, ier)
       call GatherResultsLAE(param%run_type, inode)
    ENDIF

    ! --- CSFR
    IF(inode == 0) THEN
       CSFR_out = trim(fbase)//trim(file_o)//'_CSFR.dat'
       open(1, file = trim(trim(fbase)//trim(file_o)//'_CSFR.dat'), &
            status='replace', iostat=ier)
       call CheckIerr(ier, 'main_nugc: fail to open file: '//trim(CSFR_out))
       write(1, '(A)') '# (1)z+1 (2)CSFR [h^3 Mpc^-3 Msun yr^-1]'
       temp = inv_V * param%munit / param%th_yr
       DO j = 1, mrgp%num_step
           write(1, '(2F11.4)') CSFR(1,j), CSFR(2,j) * temp
       ENDDO
       close(1)
    ENDIF

    IF(inode == 0) THEN
       include "out_norm_l9b.f90"
    ENDIF


    call MPI_BARRIER(MPI_COMM_WORLD, ier)
    DO i_file = 0, NFILE
       close(ionum+i_file)
    ENDDO

    ! delete temp. catalog files
    i = ionum + NFILE + inode + 1
    close(i,       status='delete')
    close(i+nnode, status='delete')
    IF(param%traceIDs == 1) close(i+nnode+nnode, status='delete')

    call DeallocSSP
    call DeallocateDistributionFunctions(nz_end)

    IF(param%LAE) THEN
       ! for parallel
       i = paramLAE%base + paramLAE%nfile + inode + 1
       close(i,       status='delete')
       close(i+nnode, status='delete')

       call DeallocateDistributionFunctions_LAE(nz_end)
       call DeallocSchaererArrays
    ENDIF

    ! added by MARK (2017/Mar/18)
    deallocate(lf_n_all, stat=ier); call CheckIerr(ier, trim(cerr_d)//' lf_n_all')
    deallocate(mf_n_all, stat=ier); call CheckIerr(ier, trim(cerr_d)//' mf_n_all')
    deallocate(lf_q_all, stat=ier); call CheckIerr(ier, trim(cerr_d)//' lf_q_all')
    IF(param%Mbh) THEN
       deallocate(MbhMbulge_bin,     stat=ier); call CheckIerr(ier, &
            trim(cerr_d)//' MbhMbulge_bin')
       deallocate(MbhMbulge_n_all,   stat=ier); call CheckIerr(ier, &
            trim(cerr_d)//' MbhMbulge_n_all')
       deallocate(MbhMbulge_xn_all,  stat=ier); call CheckIerr(ier, &
            trim(cerr_d)//' MbhMbulge_xn_all')
       deallocate(MbhMbulge_xxn_all, stat=ier); call CheckIerr(ier, &
            trim(cerr_d)//' MbhMbulge_xxn_all')
    ENDIF
    IF(param%ML) THEN
       deallocate(ML_bin,   stat=ier); call CheckIerr(ier, &
            trim(cerr_d)//' ML_bin')
       deallocate(ML_n_all, stat=ier); call CheckIerr(ier, &
            trim(cerr_d)//' ML_n_all')
       deallocate(ML_x_all, stat=ier); call CheckIerr(ier, &
            trim(cerr_d)//' ML_x_all')
       deallocate(ML_med,   stat=ier); call CheckIerr(ier, &
            trim(cerr_d)//' ML_med')
    ENDIF
    ! added by Makiya (2015/Nov/18)
    deallocate(CSFR,     stat=ier); call CheckIerr(ier, trim(cerr_d)//' CSFR')
    deallocate(CSFR_all, stat=ier); call CheckIerr(ier, trim(cerr_d)//' CSFR_all')
    deallocate(SMF_z,    stat=ier); call CheckIerr(ier, trim(cerr_d)//' SMF_z')
    deallocate(SMF_z_all,stat=ier); call CheckIerr(ier, trim(cerr_d)//' SMF_z_all')
    IF(param%traceIDs == 1) call DealloctraceIDs
  END SUBROUTINE CalcGalaxiesAtARedshift
!!$============================================================================
  SUBROUTINE CalcGalaxiesAtAllRedshift
    INTEGER, PARAMETER :: NzSliceMax = 200
                          ! maximum # of the redshift slices in N-body run
    INTEGER :: endhalo, end_step, nz, nz_end, i_file, num, ier
    INTEGER :: i, bin, i25, i50, i75
    INTEGER :: mor, mor_d, b_or_d
    CHARACTER(LEN=3)  :: ci
    DOUBLE PRECISION :: zsp1_max = 20.d0 ! maximum redshift to output the results
    DOUBLE PRECISION :: Mssat ! = Mstar^bulge + Mstar^disk
!!$    DOUBLE PRECISION :: bt, bt_d ! B/T ratio in B-band w/o and w/ dust
    DOUBLE PRECISION :: galmass, Mdisk, Mbulge, M_H ! for mass functions
    DOUBLE PRECISION :: MassLumi ! for M_HI/L_B distribution
!!$    DOUBLE PRECISION :: log10ML  ! for M_HI/L_B distribution
    DOUBLE PRECISION, ALLOCATABLE :: arr(:) ! for M_HI/L_B distribution
    DOUBLE PRECISION, ALLOCATABLE :: tout(:), toutGyr(:)
    ! --- functions
    INTEGER :: DetMorType
    DOUBLE PRECISION :: z2t, Square, CalBulgeReff, CalDiskReff
    DOUBLE PRECISION :: ObsFracAGN_UV
    
    CHARACTER(LEN=50) :: cerr_o = '# CalcGalaxiesAtAllRedshift: fail to open file '
    CHARACTER(LEN=50) :: cerr_a = '# CalcGalaxiesAtAllRedshift: fail to allocate '
    CHARACTER(LEN=50) :: cerr_d = '# CalcGalaxiesAtAllRedshift: fail to deallocate '

    i_file = 0
    i = i_file + ionum
    call CheckIONum(i, 'CalcGalaxiesAtAllRedshift')
    open(i, file=file_out(i_file), status='unknown', iostat=ier)
    call CheckIerr(ier, trim(cerr_o)//' '//trim(file_out(i_file)))

    ! obtain the parameters related to merger tree
    call ReadMRGB(0,0); call AllocateGalaxy

    N1 = param%nwave; N2 = 4
    N1N2Nbin = N1 * N2 * NbinLF
    allocate(lf_n_all(N1, N2, NbinLF), stat=ier); call CheckIerr(ier, &
         trim(cerr_a)//' lf_n_all'); lf_n_all(:,:,:) = 0.d0
    allocate(mf_n_all(NtypeMF, N2, NbinMF), stat=ier); call CheckIerr(ier, &
         trim(cerr_a)//' mf_n_all'); mf_n_all(:,:,:) = 0.d0
    allocate(lf_q_all(param%nwaveq, 1, NbinLF), stat=ier); call CheckIerr(ier, &
         trim(cerr_a)//' lf_q_all'); lf_q_all(:,:,:) = 0.d0
    IF(param%Mbh) THEN
       allocate(MbhMbulge_bin(50), stat=ier); call CheckIerr(ier, &
            trim(cerr_a)//' MbhMbulge_bin'); MbhMbulge_bin(:) = 0.d0
       allocate(MbhMbulge_n_all(1,   N2, 50), stat=ier); call CheckIerr(ier, &
            trim(cerr_a)//' MbhMbulge_n_all');   MbhMbulge_n_all(:,:,:)   = 0.d0
       allocate(MbhMbulge_xn_all(1,  N2, 50), stat=ier); call CheckIerr(ier, &
            trim(cerr_a)//' MbhMbulge_xn_all');  MbhMbulge_xn_all(:,:,:)  = 0.d0
       allocate(MbhMbulge_xxn_all(1, N2, 50), stat=ier); call CheckIerr(ier, &
            trim(cerr_a)//' MbhMbulge_xxn_all'); MbhMbulge_xxn_all(:,:,:) = 0.d0
    ENDIF
    IF(param%ML) THEN
       allocate(ML_bin(50), stat=ier); call CheckIerr(ier, &
            trim(cerr_a)//' ML_bin'); ML_bin(:) = 0.d0
       allocate(ML_n_all(2, 4, 50), stat=ier); call CheckIerr(ier, &
            trim(cerr_a)//' ML_n_all'); ML_n_all(:,:,:) = 0.d0
       IF(inode == 0) THEN
          allocate(ML_x_all(NgalMax*nnode, 4, 50), stat=ier)
       ELSE
          allocate(ML_x_all(NgalMax, 4, 50), stat=ier)
       ENDIF
       call CheckIerr(ier, trim(cerr_a)//' ML_x_all'); ML_x_all(:,:,:) = 0.d0
       allocate(ML_med(3, 4, 50), stat=ier); call CheckIerr(ier, &
            trim(cerr_a)//' ML_med'); ML_med(:,:,:) = 0.d0
    ENDIF

    ! --- CSFR
    allocate(CSFR(2,     mrgp%num_step), stat=ier); call CheckIerr(ier, &
       trim(cerr_a)//' CSFR');     CSFR(:,:)     = 0.d0
    allocate(CSFR_all(2, mrgp%num_step), stat=ier); call CheckIerr(ier, &
       trim(cerr_a)//' CSFR_all'); CSFR_all(:,:) = 0.d0
    CSFR(1, :)    = mrgp%zp1ar(:) ! redshift bin
    CSFR_all(1,:) = mrgp%zp1ar(:) ! redshift bin

    ! --- SMF at all z
    allocate(SMF_z(mrgp%num_step,     NbinMF), stat=ier); call CheckIerr(ier, &
       trim(cerr_a)//' SMF_z');     SMF_z(:,:)     = 0.d0
    allocate(SMF_z_all(mrgp%num_step, NbinMF), stat=ier); call CheckIerr(ier, &
       trim(cerr_a)//' SMF_z_all'); SMF_z_all(:,:) = 0.d0

    nz_end = mrgp%num_step - 1
    IF(zsp1_max < mrgp%zp1ar(mrgp%num_step)) THEN
       print '(A)', '# Choose adequate zsp1_max!!'
       print '(A, I2, A, F9.4)', '# zsp1_max must be larger than '//&
            'mrgp%zp1ar(mrgp%num_step=', mrgp%num_step, ') = ', &
            mrgp%zp1ar(mrgp%num_step)
       stop
    ELSEIF(zsp1_max < mrgp%zp1ar(1)) THEN
       DO WHILE(mrgp%zp1ar(mrgp%num_step-(nz_end-1)) >= zsp1_max)
          nz_end = nz_end - 1
       ENDDO
    ENDIF
    print '(A, I3, 2(A, F9.4))', 'nz_end = ', nz_end, &
         ': mrgp%zp1ar(nz_end) = ', mrgp%zp1ar(nz_end), &
         ', mrgp%zp1ar(nz_end-1) = ', mrgp%zp1ar(nz_end-1)
    call AllocateDistributionFunctions(nz_end)
    allocate(tout(nz_end),    stat=ier); call CheckIerr(ier, trim(cerr_a)//' tout')
    allocate(toutGyr(nz_end), stat=ier); call CheckIerr(ier, trim(cerr_a)//' toutGyr')

    IF(param%LAE) THEN ! allocate LAE related arrays
       call AllocLAE(mrgt%num_galarray, param%nwave)
       call AllocateDistributionFunctions_LAE(nz_end, param%nwave)
    ENDIF

    nforest_this = param%n_forest / nnode
    IF(mod(param%n_forest, nnode) /= 0) THEN
        nforest_this = nforest_this+1
    ENDIF
    iforest_start = inode * nforest_this

    IF(param%LAE) call DeallocLAE
    IF(param%SFH) call DeallocSFH(mrgt%num_galarray)
    call DeallocateGalaxy; call DeallocMRGB

    ID_iforest: DO iforest = iforest_start, min(iforest_start+nforest_this-1, param%n_forest-1)

       ID_nz: DO nz = 1, nz_end
          ! initialize ran1 and gasdev, for check (Makiya, 2015/09/27)
          param%idum = -131; param%iy = 0; param%iv(:) = 0; param%iset = 0

          call ReadMRGB(iforest,nz); call AllocateGalaxy
          param%loopz = nz
          param%izout = mrgp%num_step - (param%loopz - 1)
          param%zsp1  = mrgp%zp1ar(param%izout)


          tout(nz) = z2t(param%zsp1); toutGyr(nz) = tout(nz) * param%th
          const%tout = tout(nz); const%toutGyr = toutGyr(nz)
          IF(param%SFH) THEN
             call AllocateSFHRelated(param%run_type)
             call AllocSFH(mrgt%num_galarray)
          ENDIF
          IF(param%SN) call AllocateSNRelated(param%run_type)
          call ReadSSPFile(param%loopz, nz_end, nz)

          IF(iforest == iforest_start) THEN
             ! --- changing the last characters in the name of output files
             write(ci, '(I3.3)') param%izout-1; num = index(file_out(1), '.dat') - 3
             DO i_file = 1, NFILE
                file_out(i_file)(num:num+2) = ci
                i = ionum + NFILE * (nz-1) + i_file
                call CheckIONum(i, 'CalcGalaxiesAtAllRedshift')
                open(i, file = file_out(i_file), status = 'unknown', iostat=ier)
                call CheckIerr(ier, trim(cerr_o)//' '//trim(file_out(i_file)))
             ENDDO
             !-- for parallel run, added by Shirakata (2018/05/26)
             ! galaxy catalog
             num = index(file_catalog, '.dat') - 3
             file_catalog(num:num+2) = ci
             i = ionum + NFILE * nz_end + nnode * (nz-1) + inode + 1
             call CheckIONum(i, 'CalcGalaxiesAtAllRedshift')
             open(i, file=file_catalog, status='unknown', iostat=ier)
             call CheckIerr(ier, trim(cerr_o)//' '//trim(file_catalog))
             ! quasar catalog
             num = index(file_catalog_q, '.dat') - 3
             file_catalog_q(num:num+2) = ci
             i = ionum + NFILE * nz_end + nnode * (nz_end + nz - 1) + inode + 1
             call CheckIONum(i, 'CalcGalaxiesAtAllRedshift')
             open(i, file=file_catalog_q, status='unknown', iostat=ier)
          ENDIF

          IF(param%LAE) THEN
             DO i_file = 1, paramLAE%nfile
                i = i_file + paramLAE%base
                num = index(fnameLAE(i_file), '.dat') - 3
                fnameLAE(i_file)(num:num+2) = ci
                IF(iforest == 0) THEN
                   open(i, file = fnameLAE(i_file), status = 'new', iostat=ier)
                ELSE
                   open(i, file = fnameLAE(i_file), status = 'old', &
                        position = 'append', iostat=ier)
                ENDIF
                call CheckIerr(ier, trim(cerr_o)//' '//trim(fnameLAE(i_file)))
             ENDDO
          ENDIF
          IF(param%SFH) THEN
             DO i_file = 1, paramSFH%nfile
                i = i_file + paramSFH%base
                IF(iforest == 0) close(i)
                num = index(fnameSFH(i_file), '.dat') - 3
                fnameSFH(i_file)(num:num+2) = ci
                IF(iforest == 0) THEN
                   open(i, file = fnameSFH(i_file), status = 'new', iostat=ier)
                ELSE
                   open(i, file = fnameSFH(i_file), status = 'old', &
                        position = 'append', iostat=ier)
                ENDIF
                call CheckIerr(ier, trim(cerr_o)//' '//trim(fnameSFH(i_file)))
             ENDDO
          ENDIF
          IF(iforest == 0) THEN
             call WriteCaptions(ionum, nz, nz_end) ! write the captions of each data file
             IF(nz == 1) THEN
               call WriteInputParameters(1) ! write the input parameters
                                            !  into the log file of file_out(0)
               call WriteInputParameters(2) ! write the input parameters at each z
                                            !  into the log file of file_out(0)
             ENDIF
             IF(param%LAE) THEN
                call WriteCaptionsForLAE(ionum, param%run_type)
                call optdepthForLAE(ionum)
             ENDIF
             IF(param%SFH) call WriteCaptionsForSFH(ionum)
             IF(param%SN) call WriteCaptionsForSN(ionum)
          ENDIF
          call optdepth(param%zsp1, ionum, iforest, nz)

          ! calculate star formation in galaxies
          !  and return 'endhalo' and 'end_step' as the total number of
          !  galaxies at param%zsp1 and the corresponding number of param%zsp1
          !  for mrgp
          call star(iforest, endhalo, end_step, ionum, nz)
          gal(1:endhalo) = gal_next(1:endhalo)
          IF(param%LAE) allgal(1:endhalo) = allgal_next(1:endhalo)

          include "out_norm_l9a.f90"

          IF(param%LAE) call DeallocLAE
          IF(param%SFH) call DeallocSFH(mrgt%num_galarray)
          call DeallocateGalaxy; call DeallocMRGB

          call flush(ionum+NFILE*nz_end+nnode*(nz-1)+inode+1)
          call flush(ionum+NFILE*nz_end+nnode*(nz_end+nz-1)+inode+1)

          IF(param%LAE) THEN
             DO i_file = 1, paramLAE%nfile
                close(i_file+paramLAE%base)
             ENDDO
          ENDIF
          IF(param%SFH) THEN
             call DeallocateSFHRelated(param%run_type)
             call DeallocSFH(mrgt%num_galarray)
             DO i_file = 1, paramSFH%nfile
                close(i_file+paramSFH%base)
             ENDDO
          ENDIF
          call DeallocSSP

          print *
       ENDDO ID_nz
    ENDDO ID_iforest 

    call MPI_BARRIER(MPI_COMM_WORLD, ier)
    DO nz = 1, nz_end
       call GatherResults(nz)
    ENDDO

    IF(param%LAE) THEN
       call flush(paramLAE%base + paramLAE%nfile + inode + 1)
       call flush(paramLAE%base + paramLAE%nfile + nnode + inode + 1)
       call MPI_BARRIER(MPI_COMM_WORLD, ier)
       call GatherResultsLAE(param%run_type, inode)
    ENDIF

    deallocate(tout,    stat=ier); call CheckIerr(ier, trim(cerr_d)//' tout')
    deallocate(toutGyr, stat=ier); call CheckIerr(ier, trim(cerr_d)//' toutGyr')

    ! added by MARK (2017/Mar/18)
    deallocate(lf_n_all, stat=ier); call CheckIerr(ier, trim(cerr_d)//' lf_n_all')
    deallocate(mf_n_all, stat=ier); call CheckIerr(ier, trim(cerr_d)//' mf_n_all')
    deallocate(lf_q_all, stat=ier); call CheckIerr(ier, trim(cerr_d)//' lf_q_all')
    IF(param%Mbh) THEN
       deallocate(MbhMbulge_bin,     stat=ier); call CheckIerr(ier, &
            trim(cerr_d)//' MbhMbulge_bin')
       deallocate(MbhMbulge_n_all,   stat=ier); call CheckIerr(ier, &
            trim(cerr_d)//' MbhMbulge_n_all')
       deallocate(MbhMbulge_xn_all,  stat=ier); call CheckIerr(ier, &
            trim(cerr_d)//' MbhMbulge_xn_all')
       deallocate(MbhMbulge_xxn_all, stat=ier); call CheckIerr(ier, &
            trim(cerr_d)//' MbhMbulge_xxn_all')
    ENDIF
    IF(param%ML) THEN
       deallocate(ML_bin,   stat=ier); call CheckIerr(ier, &
            trim(cerr_d)//' ML_bin')
       deallocate(ML_n_all, stat=ier); call CheckIerr(ier, &
            trim(cerr_d)//' ML_n_all')
       deallocate(ML_x_all, stat=ier); call CheckIerr(ier, &
            trim(cerr_d)//' ML_x_all')
       deallocate(ML_med,   stat=ier); call CheckIerr(ier, &
            trim(cerr_d)//' ML_med')
    ENDIF

    ! --- CSFR
    IF(inode == 0) THEN
       CSFR_out = trim(fbase)//trim(file_o)//'_CSFR.dat'
       open(1, file = trim(trim(fbase)//trim(file_o)//'_CSFR.dat'), &
          status='replace', iostat=ier)
       call CheckIerr(ier, 'main_nugc: fail to open file: '//trim(CSFR_out))
       write(1, '(A)') '# (1)z+1 (2)CSFR [h^3 Mpc^-3 Msun yr^-1]'
       temp = inv_V * param%munit / param%th_yr
       DO j = 1, mrgp%num_step
          write(1, '(2F11.4)') CSFR(1,j), CSFR(2,j) * temp
       ENDDO
       close(1)
    ENDIF

    DO nz = 1, nz_end
       write(ci, '(I3.3)') mrgp%num_step - nz
       num = index(file_out(1), '.dat') - 3
       DO i = 1, NFILE
          file_out(i)(num:num+2) = ci
          open(ionum+NFILE*(nz-1)+i, file = file_out(i), status = 'old', &
               position = 'append', iostat=ier)
          call CheckIerr(ier, trim(cerr_o)//' '//trim(file_out(i)))
       ENDDO

       IF(inode == 0) THEN
          include "out_norm_l9b.f90"
       ENDIF

       call MPI_BARRIER(MPI_COMM_WORLD, ier)
       DO i_file = 1, NFILE
          close(ionum+NFILE*(nz-1)+i_file)
       ENDDO

       IF(nz == 1) THEN
         i_file = 0
         close(i_file+ionum)
       ENDIF

       ! delete temp. calalog files
       i = ionum + NFILE * nz_end + nnode * (nz-1) + inode + 1
       close(i,       status='delete')
       close(i+nnode*nz_end, status='delete')
    ENDDO
    call DeallocateDistributionFunctions(nz_end)
    IF(param%LAE) THEN
       ! for parallel
       i = paramLAE%base + paramLAE%nfile + inode + 1
       close(i,       status='delete')
       close(i+nnode, status='delete')

       call DeallocateDistributionFunctions_LAE(nz_end)
       call DeallocSchaererArrays
    ENDIF
    deallocate(CSFR,     stat=ier); call CheckIerr(ier, trim(cerr_d)//' CSFR')
    deallocate(CSFR_all, stat=ier); call CheckIerr(ier, trim(cerr_d)//' CSFR_all')
    deallocate(SMF_z,    stat=ier); call CheckIerr(ier, trim(cerr_d)//' SMF_z')
    deallocate(SMF_z_all,stat=ier); call CheckIerr(ier, trim(cerr_d)//' SMF_z_all')
  END SUBROUTINE CalcGalaxiesAtAllRedshift
!!$============================================================================
  SUBROUTINE CalcMCMC(nop)
    use MCMCrelated

    INTEGER, INTENT(IN) :: nop
    INTEGER :: endhalo, end_step, nz, nz_end, i_file, ier
    INTEGER :: i, bin, i25, i50, i75
    INTEGER :: mor, mor_d, b_or_d
    CHARACTER(LEN=500) :: fname
    DOUBLE PRECISION :: Mssat ! = Mstar^bulge + Mstar^disk
    DOUBLE PRECISION :: galmass, Mdisk, Mbulge, M_H ! for mass functions
    DOUBLE PRECISION :: MassLumi ! for M_HI/L_B distribution
    DOUBLE PRECISION, ALLOCATABLE :: arr(:) ! for M_HI/L_B distribution
    INTEGER :: dof ! degree of freedom & No. of free parameters
    DOUBLE PRECISION :: zsp1_max = 20.d0 ! maximum redshift to output the results
    DOUBLE PRECISION, ALLOCATABLE :: tout(:), toutGyr(:)

    ! --- functions
    INTEGER :: DetMorType
    DOUBLE PRECISION :: Square, AbsDiff, z2t, CalBulgeReff, CalDiskReff
    DOUBLE PRECISION :: ObsFracAGN_UV

    ! --- for parallel
    !INTEGER :: nforest_this, iforest_start
    !INTEGER :: N1, N2, N1N2Nbin
    INTEGER :: clock,nRand
    INTEGER, ALLOCATABLE :: seed(:)

    CHARACTER(LEN=50) :: cerr_a = '# CalcMCMC: fail to allocate '
    CHARACTER(LEN=50) :: cerr_d = '# CalcMCMC: fail to deallocate '

    call MPI_Comm_size(MPI_COMM_WORLD, nnode, istat)
    call MPI_Comm_rank(MPI_COMM_WORLD, inode, istat)

    call random_seed(size = nRand)
    allocate(seed(nRand), stat=ier)
    call system_clock(count=clock)
    seed = clock
    call random_seed(put=seed)

    !--- write input parameters ---
    IF(inode == 0) THEN
       open(ionum, file = file_out(0), status = 'unknown', iostat=ier)
       call CheckIerr(ier, 'CalcMCMC: fail to open file: '//&
            trim(file_out(0)))
       write(ionum, '(A)') '# log file: the adopted parameter values in '//&
            'the calculation are written'
       call WriteInputParameters(1) ! write the input parameters
                                    !  into the log file of file_out(0)
       call WriteInputParameters(2) ! write the input parameters
       call WriteInputParameters(3)
    ENDIF

    call ReadObsData

    call ReadMRGB(0,0); call AllocateGalaxy

    N1 = param%nwave; N2 = 4
    N1N2Nbin = N1 * N2 * NbinLF
    allocate(lf_n_all(N1, N2, NbinLF), stat=ier); call CheckIerr(ier, &
         trim(cerr_a)//' lf_n_all'); lf_n_all(:,:,:) = 0.d0
    allocate(mf_n_all(NtypeMF, N2, NbinMF), stat=ier); call CheckIerr(ier, &
         trim(cerr_a)//' mf_n_all'); mf_n_all(:,:,:) = 0.d0
    allocate(lf_q_all(param%nwaveq, 1, NbinLF), stat=ier); call CheckIerr(ier, &
         trim(cerr_a)//' lf_q_all'); lf_q_all(:,:,:) = 0.d0

    ! Mbh -- Mbulge
    allocate(MbhMbulge_bin(50), stat=ier); call CheckIerr(ier, &
         trim(cerr_a)//' MbhMbulge_bin')
    allocate(MbhMbulge_n_all(1,   N2, 50), stat=ier); call CheckIerr(ier, &
         trim(cerr_a)//' MbhMbulge_n_all')
    allocate(MbhMbulge_xn_all(1,  N2, 50), stat=ier); call CheckIerr(ier, &
         trim(cerr_a)//' MbhMbulge_xn_all')
    allocate(MbhMbulge_xxn_all(1, N2, 50), stat=ier); call CheckIerr(ier, &
         trim(cerr_a)//' MbhMbulge_xxn_all')
    ! Rd -- Vd 
    allocate(DiskScale_bin(50), stat=ier); call CheckIerr(ier, &
         trim(cerr_a)//' DiskScale_bin')
    allocate(DiskScale_n_all(1,   N2, 50), stat=ier); call CheckIerr(ier, &
         trim(cerr_a)//' DiskScale_n_all')
    allocate(DiskScale_xn_all(1,  N2, 50), stat=ier); call CheckIerr(ier, &
         trim(cerr_a)//' DiskScale_xn_all')
    allocate(DiskScale_xxn_all(1, N2, 50), stat=ier); call CheckIerr(ier, &
         trim(cerr_a)//' DiskScale_xxn_all')

    ! Fundamental plane of spheroids 
    allocate(SphScale_bin(50), stat=ier); call CheckIerr(ier, &
         trim(cerr_a)//' SphScale_bin')
    allocate(SphScale_n_all(1,   N2, 50), stat=ier); call CheckIerr(ier, &
         trim(cerr_a)//' SphScale_n_all')
    allocate(SphScale_xn_all(1,  N2, 50), stat=ier); call CheckIerr(ier, &
         trim(cerr_a)//' SphScale_xn_all')
    allocate(SphScale_xxn_all(1, N2, 50), stat=ier); call CheckIerr(ier, &
         trim(cerr_a)//' SphScale_xxn_all')

    IF(param%ML) THEN
       allocate(ML_bin(50), stat=ier); call CheckIerr(ier, &
            trim(cerr_a)//' ML_bin')
       allocate(ML_n_all(2, N2, 50), stat=ier); call CheckIerr(ier, &
            trim(cerr_a)//' ML_n_all')
       allocate(ML_x_all(NgalMax*nnode, N2, 50), stat=ier); call CheckIerr(ier, &
            trim(cerr_a)//' ML_x_all')
       allocate(ML_med(3, N2, 50), stat=ier); call CheckIerr(ier, &
            trim(cerr_a)//' ML_med')
    ENDIF

    ! --- CSFR, by Makiya
    allocate(CSFR(2,     mrgp%num_step), stat=ier); call CheckIerr(ier, &
         trim(cerr_a)//' CSFR');     CSFR(:,:)     = 0.d0
    allocate(CSFR_all(2, mrgp%num_step), stat=ier); call CheckIerr(ier, &
         trim(cerr_a)//' CSFR_all'); CSFR_all(:,:) = 0.d0
    CSFR(1,:)     = mrgp%zp1ar(:) ! redshift bin
    CSFR_all(1,:) = mrgp%zp1ar(:) ! redshift bin

    ! --- SMF at all z, by Makiya
    allocate(SMF_z(mrgp%num_step,     NbinMF), stat=ier); call CheckIerr(ier, &
         trim(cerr_a)//' SMF_z');     SMF_z(:,:)     = 0.d0
    allocate(SMF_z_all(mrgp%num_step, NbinMF), stat=ier); call CheckIerr(ier, &
         trim(cerr_a)//' SMF_z_all'); SMF_z_all(:,:) = 0.d0

    nz_end = mrgp%num_step - 1
    IF(zsp1_max < mrgp%zp1ar(mrgp%num_step)) THEN
       print '(A)', '# Choose adequate zsp1_max!!'
       print '(A, I2, A, F9.4)', '# zsp1_max must be larger than '//&
            'mrgp%zp1ar(mrgp%num_step=', mrgp%num_step, ') = ', &
            mrgp%zp1ar(mrgp%num_step)
       stop
    ELSEIF(zsp1_max < mrgp%zp1ar(1)) THEN
       DO WHILE(mrgp%zp1ar(mrgp%num_step-(nz_end-1)) >= zsp1_max)
          nz_end = nz_end - 1
       ENDDO
    ENDIF

    call AllocateDistributionFunctions(nz_end)
    allocate(tout(nz_end),    stat=ier); call CheckIerr(ier, trim(cerr_a)//' tout')
    allocate(toutGyr(nz_end), stat=ier); call CheckIerr(ier, trim(cerr_a)//' toutGyr')

    nforest_this = param%n_forest / nnode
    IF(mod(param%n_forest, nnode) /= 0) THEN
        nforest_this = nforest_this+1
    ENDIF
    iforest_start = inode * nforest_this
    
    call DeallocateGalaxy; call DeallocMRGB

    chisqtot = 0.d0; chisq_tmp = 0.d0

    IF(inode == 0) THEN
       print *, ""
       print *, "++++++++++++++++++++++++++"
       print *, "+ MCMC fitting start!!!! +"
       print *, "++++++++++++++++++++++++++"
    ENDIF

    ! call InitParams
    dof = 0
    DO iMCMC = 1, mcmc%iMCMC_max
       IF(inode == 0) print *, "iMCMC :: ", iMCMC

       call MPI_BARRIER(MPI_COMM_WORLD, ier)
       ! --- Synthesize the paranmeters of each node
       call SynthesizeParams

       !--- Initialize distribution functions ---
       ! Mbh -- Mbulge
       MbhMbulge%n(:,:,:)   = 0.d0; MbhMbulge%xn(:,:,:) = 0.d0
       MbhMbulge%xxn(:,:,:) = 0.d0

       ! Rd -- Vd 
       DiskScale%n(:,:,:)   = 0.d0; DiskScale%xn(:,:,:) = 0.d0
       DiskScale%xxn(:,:,:) = 0.d0

       ! Fundamental plane
       SphScale%n(:,:,:)   = 0.d0; SphScale%xn(:,:,:) = 0.d0
       SphScale%xxn(:,:,:) = 0.d0

       CSFR(2,:) = 0.d0
       SMF_z(:,:) = 0.d0

       DO nz= 1, nz_end
          lf(nz)%n(:,:,:)      = 0.d0; lf_d(nz)%n(:,:,:)      = 0.d0
          lf(nz)%n_brst(:,:,:) = 0.d0; lf_d(nz)%n_brst(:,:,:) = 0.d0
          mf(nz)%n(:,:,:) = 0.d0; mf(nz)%n_brst(:,:,:) = 0.d0
          IF(param%ML) THEN
             ml(nz)%n(:,:,:)   = 0.d0; ml(nz)%xn(:,:,:) = 0.d0
             ml(nz)%xxn(:,:,:) = 0.d0; ml(nz)%x(:,:,:)  = 0.d0
          ENDIF
          lf_q(nz)%n(:,:,:) = 0.d0
       ENDDO
       !-----------------------------------------

       ! if new parameters exceeds the bounds of prior, they are rejected
       !  and calculation of galaxy formation is skipped
       IF(param_bounds_key == 0) THEN
          ID_iforest:DO iforest = iforest_start, &
                        min(iforest_start+nforest_this-1, param%n_forest-1)
             ID_nz: DO nz = 1, nz_end
                IF(nz == 1 .or. nz == 8 .or. nz == 17 &
                   .or. nz == 26 .or. nz == 32) THEN
                ! z ~ 0, 0.4, 1, 2, 3 (nugc-SS sim.)
                   call ReadMRGB(iforest,nz); call AllocateGalaxy
                   ! initialize ran1 and gasdev
                   param%idum = -131; param%iy = 0; param%iv(:) = 0; param%iset = 0
                   param%loopz = nz
                   param%izout = mrgp%num_step - (param%loopz - 1)
                   param%zsp1  = mrgp%zp1ar(param%izout)

                   tout(nz) = z2t(param%zsp1); toutGyr(nz) = tout(nz) * param%th
                   const%tout = tout(nz); const%toutGyr = toutGyr(nz)

                   call ReadSSPFile(param%loopz, nz_end, nz)
                   call optdepth(param%zsp1, ionum, 0, nz)

                   ! calculate star formation in galaxies
                   !  and return 'endhalo' and 'end_step' as the total number of
                   !  galaxies at param%zsp1 and the corresponding number of param%zsp1
                   !  for mrgp
                   call star(iforest, endhalo, end_step, ionum, nz)
                   gal(1:endhalo) = gal_next(1:endhalo)

                   include "out_norm_l9a.f90"

                   call DeallocateGalaxy; call DeallocMRGB
                   call DeallocSSP
                ENDIF
             ENDDO ID_nz
          ENDDO ID_iforest

          call MPI_BARRIER(MPI_COMM_WORLD, ier)
          ! --- gather the calculation results of each node
          DO nz = 1, nz_end
             call GatherResults(nz)
          ENDDO

          IF(inode == 0) THEN
             ! --- calculate the chi^2 of this run
             call CalcChisq(inv_V, dof)
             p_tmp = exp(-0.5*(chisqtot-chisq_tmp))
             call random_number(d)
             IF(p_tmp > d .or. iMCMC == 1) THEN ! adopt new parameter set
                IF(inode == 0) print *, "Accept new parameter set"
                chisq_tmp = chisqtot

                alpst_tmp   = param%alpst
                Vst_tmp     = param%Vst
                tau0st_tmp  = param%tau0st; 
                eps_SF_tmp = param%eps_SF
                emin_tmp    = param%emin
                emax_tmp    = param%emax
                Zch_tmp    = param%Zch
                taud_th_tmp    = param%taud_th

                Vhot2_tmp   = param%Vhot(2);   Vhot1_tmp   = param%Vhot(1)
                alphot2_tmp = param%alphot(2); alphot1_tmp = param%alphot(1)
                alp_ret_tmp = param%alp_ret

                fmerge_tmp     = param%fmerge
                fmajor_tmp     = param%fmajor
                Krem_tmp       = param%Krem
                fdiss_tmp      = param%fdiss
                Mh0_tmp        = param%Mh0
                alpha_halo_tmp = param%alpha_halo

                Vcut_tmp    = param%Vcut
                kappa_tmp   = param%kappa_croton; eta_tmp = param%eta_croton
                alp_bow_tmp = param%alp_bower;    eps_bow_tmp = param%eps_bower
                Mbhseed_tmp = param%Mbhseed

                em_tmp      = param%em
                fdi_tmp     = param%fdi

                fbh1_tmp     = param%fbh(1); fbh2_tmp = param%fbh(2)
                Mbhseed_tmp  = param%Mbhseed
                tacc0_tmp    = param%tacc_0; tacc1_tmp = param%tacc_1
                alpha_ad_tmp = param%alpha_ad
                gamma_ad_tmp = param%gamma_ad
                Edd_max_tmp  = param%Edd_max

                ! output temp. results of LF, MF...
                ! DO nz = 1, nz_end
                !    call WriteTmpResults(nz)
                ! ENDDO
             ELSE ! reject new parameter set
                IF(inode == 0) print *, "Reject new parameter set"
             ENDIF
          ENDIF
       ELSE
          IF(inode == 0) print *, "Reject new parameter set (param outside bounds)"
       ENDIF ! if param_bounds_key end

       IF(inode == 0) THEN
          IF(iMCMC == 1) THEN
             dof = dof - nop
             print *, 'dof = ', dof
             IF(dof <= 0) THEN
                print *, 'dof <= 0 so dof is fixed as 1.'
                dof = 1
             ENDIF
             i_file = 98; fname = 'history_'//trim(file_o)//'.dat'
             open(i_file, file = trim(fname), status = 'replace')
             write(i_file, '(A, I3)') '# d.o.f = ', dof
             write(i_file, '(A)') '# (1)iteration (2)reduced chi^2 (3)alpst '// &
                  '(4)Vst (5)tau0st (6)eps_SF (7)emax (8)emin (9)Zch '//&
                  '(10)taud_th (11)Vhot(1) (12)Vhot(2) (13)alphot(1) '//&
                  '(14)alphot(2) (15)alp_ret (16)fmerge (17)fmajor '//&
                  '(18)Krem (19)fdiss (20)Mh0 (21)alpha_halo '//&
                  '(22)Vcut (23)kappa_croton (24)eta_croton '//&
                  '(25)alp_bower (26)eps_bower (27)em (28)fdi '//&
                  '(29)fbh(1) (30)fbh(2) (31)Mbhseed (32)tacc_0 '//&
                  '(33)tacc_1 (34)alpha_ad (35)gamma_ad (36)Edd_max'
          ENDIF

          write(i_file, '(I6, E12.4e2, F9.4, E12.4e2, 6F9.4, 2E12.4e2, 9F9.4, E12.4e2, 8F9.4, E12.4e2, 5F9.4)') &
               iMCMC, chisq_tmp/dof, alpst_tmp, Vst_tmp, tau0st_tmp*param%th, eps_SF_tmp, &
               emin_tmp, emax_tmp, Zch_tmp, taud_th_tmp, & 
               Vhot1_tmp, Vhot2_tmp, alphot1_tmp, alphot2_tmp, &
               alp_ret_tmp, fmerge_tmp, fmajor_tmp, Krem_tmp, fdiss_tmp, &
               Mh0_tmp, alpha_halo_tmp, Vcut_tmp, kappa_tmp, eta_tmp, &
               alp_bow_tmp, eps_bow_tmp, em_tmp, fdi_tmp, &
               fbh1_tmp, fbh2_tmp, Mbhseed_tmp, tacc0_tmp, &
               tacc1_tmp, alpha_ad_tmp, gamma_ad_tmp, Edd_max_tmp
          call flush(i_file)
       ENDIF

       !! Renewal parameters
       IF(inode == 0) THEN
          call RenewalParams
          print *, "# =*=*=*=*=*=*=*=*=*=*=*= #"
          ! --- If you want to calculate mean LF and MF,
          ! --- you should call following two subroutines.
          ! --- for MCMC, you must comment out them.
          ! call RenewalParams_mean
          ! call WriteMeanResults(nz, iMCMC)
       ENDIF
    ENDDO
    !--- MCMC loop end---!

    IF(inode == 0) THEN
       include "out_norm_l9b.f90"
    ENDIF

    close(ionum)

    call DeallocateDistributionFunctions(nz_end)


    deallocate(tout,    stat=ier); call CheckIerr(ier, trim(cerr_d)//' tout')
    deallocate(toutGyr, stat=ier); call CheckIerr(ier, trim(cerr_d)//' toutGyr')

    ! added by MARK (2017/Mar/18)
    deallocate(seed,     stat=ier); call CheckIerr(ier, trim(cerr_d)//' seed')
    deallocate(lf_n_all, stat=ier); call CheckIerr(ier, trim(cerr_d)//' lf_n_all')
    deallocate(mf_n_all, stat=ier); call CheckIerr(ier, trim(cerr_d)//' mf_n_all')
    deallocate(lf_q_all, stat=ier); call CheckIerr(ier, trim(cerr_d)//' lf_q_all')

    deallocate(MbhMbulge_bin,     stat=ier); call CheckIerr(ier, &
         trim(cerr_d)//' MbhMbulge_bin')
    deallocate(MbhMbulge_n_all,   stat=ier); call CheckIerr(ier, &
         trim(cerr_d)//' MbhMbulge_n_all')
    deallocate(MbhMbulge_xn_all,  stat=ier); call CheckIerr(ier, &
         trim(cerr_d)//' MbhMbulge_xn_all')
    deallocate(MbhMbulge_xxn_all, stat=ier); call CheckIerr(ier, &
         trim(cerr_d)//' MbhMbulge_xxn_all')

    deallocate(DiskScale_bin,     stat=ier); call CheckIerr(ier, &
         trim(cerr_d)//' DiskScale_bin')
    deallocate(DiskScale_n_all,   stat=ier); call CheckIerr(ier, &
         trim(cerr_d)//' DiskScale_n_all')
    deallocate(DiskScale_xn_all,  stat=ier); call CheckIerr(ier, &
         trim(cerr_d)//' DiskScale_xn_all')
    deallocate(DiskScale_xxn_all, stat=ier); call CheckIerr(ier, &
         trim(cerr_d)//' DiskScale_xxn_all')

    deallocate(SphScale_bin,     stat=ier); call CheckIerr(ier, &
         trim(cerr_d)//' SphScale_bin')
    deallocate(SphScale_n_all,   stat=ier); call CheckIerr(ier, &
         trim(cerr_d)//' SphScale_n_all')
    deallocate(SphScale_xn_all,  stat=ier); call CheckIerr(ier, &
         trim(cerr_d)//' SphScale_xn_all')
    deallocate(SphScale_xxn_all, stat=ier); call CheckIerr(ier, &
         trim(cerr_d)//' SphScale_xxn_all')

    IF(param%ML) THEN
       deallocate(ML_bin,   stat=ier); call CheckIerr(ier, &
            trim(cerr_d)//' ML_bin')
       deallocate(ML_n_all, stat=ier); call CheckIerr(ier, &
            trim(cerr_d)//' ML_n_all')
       deallocate(ML_x_all, stat=ier); call CheckIerr(ier, &
            trim(cerr_d)//' ML_x_all')
       deallocate(ML_med,   stat=ier); call CheckIerr(ier, &
            trim(cerr_d)//' ML_med')
    ENDIF
    ! added by Makiya (2015/Nov/18)
    deallocate(CSFR,     stat=ier); call CheckIerr(ier, trim(cerr_d)//' CSFR')
    deallocate(CSFR_all, stat=ier); call CheckIerr(ier, trim(cerr_d)//' CSFR_all')
    deallocate(SMF_z,    stat=ier); call CheckIerr(ier, trim(cerr_d)//' SMF_z')
    deallocate(SMF_z_all,stat=ier); call CheckIerr(ier, trim(cerr_d)//' SMF_z_all')
  END SUBROUTINE CalcMCMC
!!$============================================================================
  SUBROUTINE ReadCoolingData
    use global_var
    implicit none
    INTEGER :: i
    DOUBLE PRECISION :: T, x


    open(1, file = 'inpdata/coolfn.dat', status = 'old'); i = 1
100 read(1, *, end = 101) T, x
    lambda(i) = dble(x)
    i = i + 1
    goto 100
101 continue
    close(1)

    open(1, file = 'inpdata/coolsol.dat', status = 'old'); i = 1
102 read(1, *, end = 103) T, x
    lambdas(i) = dble(x)
    i = i + 1
    goto 102
103 continue
    close(1)
  END SUBROUTINE ReadCoolingData
!!$============================================================================
  SUBROUTINE ReadParameters
    INTEGER :: i, ier, i_file
    CHARACTER(LEN=500) :: buf
    CHARACTER(LEN=50) :: cerr_o = '# ReadParameters: fail to open file'
    CHARACTER(LEN=50) :: cerr_a = '# ReadParameters: fail to allocate'


    param%munit = 1.d+14
    param%log10munit = log10(param%munit)

    i_file = 1
    open(i_file, file='inpdata/'//trim(file_i)//'.dat', status='old', iostat=ier)
    call CheckIerr(ier, trim(cerr_o)//' inpdata/'//trim(file_i)//&
         '.dat')

    read(i_file, '(A)') buf
    i = index(buf, ' ')
    param%file_nbody = buf(1:i-1)
    read(i_file, *) param%run_type ! 1:run at a redshift
                                   ! 2:run at all redshift
                                   ! 3:MCMC
    read(i_file, *) param%traceIDs ! 0: off 1: on
    read(i_file, *) param%idum     ! seed for random nuber
    ! read(i_file, *) param%zsp1      ! output redshift
    read(i_file, *) param%zsp1_input ! output redshift
                                     ! (modified by MARK on 2016/06/28)
    param%zsp1 = param%zsp1_input   ! added by MARK (2016/06/28)
    read(i_file, *) param%alpst     ! index of SF timescale
    read(i_file, *) param%Vst       ! normalization  of SF timescale
    read(i_file, *) param%Vhot(1)   ! normalization of SN feedback for starubrst
    read(i_file, *) param%Vhot(2)   !                                  quiescent
    read(i_file, *) param%alphot(1) ! index of SN feedback for starburst
    read(i_file, *) param%alphot(2) !                          quiescent
    read(i_file, *) param%SFmodel   ! 1:CSF, 2:DSF
    read(i_file, *) param%tau0st    ! minimum SF timescale for quiescent
    read(i_file, *) param%eps_SF    ! SF efficiency in DSF model
    read(i_file, *) param%emax      ! maximum SFE in Makiya14 model
    read(i_file, *) param%emin      ! minimum SFE in Makiya14 model
    read(i_file, *) param%Zch       ! characteristic metallicity in Makiya14 model
    read(i_file, *) param%taud_th   ! Dust opacity threshold in Makiya14 model
    read(i_file, *) param%fmerge    ! parameter for the timescale of dyn. fric.
    read(i_file, *) param%fmajor    ! parameter for mass ratio of merging galaxy 
    read(i_file, *) param%Krem      ! parameter for the remained energy fraction
    read(i_file, *) param%alp_ret   ! gas mass fraction returned from reservoir

    ! --- logical parameters
    read(i_file, *) param%off_collision     ! for random collision
    read(i_file, *) param%equal_mass_merger ! for euqal-mass merger
    read(i_file, *) param%dyn_resp_bulge    ! for dyn. response for bulge
    read(i_file, *) param%dyn_resp_disk     ! for dyn. response for disk
    read(i_file, *) param%dyn_resp_halo     ! for dyn. response for halo

    read(i_file, *) param%freheat ! param. for halo major merger
    read(i_file, *) param%fb      ! param. for bulge size
    read(i_file, *) param%Vcut    ! cooling cut-off circular vel. (max.)
    read(i_file, *) param%Vlow    ! cooling cut-off circular vel. (min.)
    read(i_file, *) param%tauV0   ! V-band dust opacity
                                  !  (original: 9.d+5 from Galactic Dynamics)
    read(i_file, *) param%exttype ! dust ext. curve
                                  !  (1:MW, 2:LMC, 3:SMC, 4: Calzetti-law)
    read(i_file, *) param%fdm     ! param. for Vcent, DM contribution
    read(i_file, *) param%fdiss   ! param. for dissipation fraction
    read(i_file, *) param%em      ! param. for disk instability
    read(i_file, *) param%fdi     ! param. for disk instability
    read(i_file, *) param%alpha_tau ! param. for redshift evolution of dust opacity
    read(i_file, *) param%alpha_rad   ! param. for redshift evolution of disk size
    read(i_file, *) param%alpha_burst ! param. for dust opacity in starburst
    read(i_file, *) param%Mh0         ! param. for Bulbe velocity and size 
    read(i_file, *) param%alpha_halo  ! param. for Bulbe velocity and size
    read(i_file, *) param%f_ngal    ! param. for the number of elements of
                                    !  galaxies' arrays (added by MARK,
                                    !  2016/Jun/29)
    read(i_file, *) param%UVfb      ! UV feedback (1:Vlow, 2:Okamoto+08)
    read(i_file, *) param%zreion    ! reionization redshift
                                    !  (only used at param%UVfb=2)
    IF(param%UVfb == 2) THEN ! Okamoto+08 UV feedback mode
       param%iO08_p  = 0; param%iO08_n0 = 0; param%iO08_n1 = 0
       param%iO08_n2 = 0; param%iO08_n3 = 0; param%iO08_n4 = 0
    ENDIF
    read(i_file, *) param%CoolFN    ! Cooling function (Added by Shirakata) 
                                    ! 1: Sotherland & Dopita, 2: Okamoto cloudy

    ! --- SMBH related parameters
    read(i_file, *) param%fbh(1)     ! accreting mass fraction (default: 0.01d0)
    read(i_file, *) param%fbh(2)     ! accreting mass fraction (default: 0.01d0)
    read(i_file, *) param%Mbhseed ! mass of the seed BH (default: 0.d0)
    read(i_file, *) param%eps_agn ! radiative efficiency in B-band
                                  !  (default: 0.00551d0)
    read(i_file, *) param%tacc_0  ! AGN lifetime [10^7 Myr]
    read(i_file, *) param%tacc_1  ! AGN lifetime [10^7 Myr]
    read(i_file, *) param%acc     ! AGN light curve related
    read(i_file, *) param%alpha_ad  ! AGN light curve related
    read(i_file, *) param%gamma_ad  ! AGN light curve related
    read(i_file, *) param%tad_0     ! AGN light curve related
    read(i_file, *) param%nwaveq  ! # of bands for AGN

    ! --- AGN feedback related parameters
    read(i_file, *) param%AGNFB_key ! 1:Vcut, 2:Croton, 3:Bower
    read(i_file, *) param%kappa_croton
    read(i_file, *) param%eta_croton
    read(i_file, *) param%alp_bower
    read(i_file, *) param%eps_bower
!    IF(param%AGNFB_key == 3) THEN ! Bower+06 mode AGN feedback
       param%Lcool0 = 1.5d+41  * param%munit * const%kB / (param%mu * const%mp)
       param%Ledd0  = 1.26d+38 * param%munit
!    ENDIF
    ! --- parameters (added by MARK)
    read(i_file, *) param%ML  ! for ML calculation
    read(i_file, *) param%Mbh ! for Mbh as a function of Mbulge calculation
    read(i_file, *) param%LAE ! for LAE related calculation
    read(i_file, *) param%SFH ! for SFH related calculation
    read(i_file, *) param%SN  ! for SN related calculation
    IF((param%SFH .eqv. .false. .or. param%LAE .eqv. .false.) &
         .and. param%SN) THEN
       print '(A)', '# both param%SFH and param%LAE should be .true. '//&
            'if you want to calculate SN related quantities'; stop
    ENDIF

    ! --- surface brightness limit related parameters
    !      (added by MARK, 2015/Apr/11)
    read(i_file, *) param%SB      ! for LF calculation w/ SB limit
    read(i_file, *) param%SBlimit ! SB limit for LF calculation
    IF(param%SB .eqv. .false.) param%SBlimit = 100.d0
    read(i_file, *) param%DI
    read(i_file, *) param%ML_AGN ! 1: L = eMc^2 2: Watarai+00
    read(i_file, *) param%Edd_max 
    read(i_file, *) param%BolCor  ! 1: Marconi+04
    read(i_file, *) param%ObsFrac ! 1: Shirakata+18 2: Hopkins+07
    read(i_file, *) param%delay
    ! --- ssp related parameters
    read(i_file, *) param%run_all_ssp ! .true. : run all SSPs
                                      ! .false.: run a certain SSP
    read(i_file, *) ssp%filename(1) ! ssp filename for starburst
    read(i_file, *) ssp%filename(2) ! ssp filename for quiescent
    IF(ssp%filename(2)(1:1) /= ssp%filename(1)(1:1)) THEN
       print '(A)', '# The Name of sspfile Is Inappropriate!!!'; stop
    ELSE
       IF(ssp%filename(1)(1:1) == 'K') THEN
          param%type_ssp = 1
       ELSEIF(ssp%filename(1)(1:1) == 'P') THEN
          param%type_ssp = 2
       ELSEIF(ssp%filename(1)(1:1) == 'B') THEN
          param%type_ssp = 3
       ENDIF
    ENDIF
    ! added by MARK (2014/Feb/02)
    param%type_mag = index(trim(ssp%filename(1)), 'AB') + 1 ! 1 for Vega mag.
    IF(param%type_mag /= 1) param%type_mag = 2 ! AB mag.
    read(i_file, *) param%nwave ! # of calculated filters
    IF(param%nwave < 1 .or. param%nwave > TNWAVE_ALL) THEN
       print '(A, I3, A)', &
            '# The Number of Filters Should Be 1-', TNWAVE_ALL, '!!!'; stop
    ELSE
       param%nwp1 = param%nwave + 1; param%tnw  = 2 * param%nwave
       allocate(param%iwave(param%nwave), stat=ier); call CheckIerr(ier, &
            trim(cerr_a)//' param%iwave')
       allocate(mag(param%nwave),    stat=ier); call CheckIerr(ier, &
            trim(cerr_a)//' mag')
       allocate(mag_d(param%nwave),  stat=ier); call CheckIerr(ier, &
            trim(cerr_a)//' mag_d')
       allocate(xi(param%nwp1),      stat=ier); call CheckIerr(ier, &
            trim(cerr_a)//' xi')
       allocate(lumtmp(param%tnw),   stat=ier); call CheckIerr(ier, &
            trim(cerr_a)//' lumtmp')
       allocate(lumtmp_d(param%tnw), stat=ier); call CheckIerr(ier, &
            trim(cerr_a)//' lumtmp_d')
       IF(param%nwave /= TNWAVE_ALL) THEN
          read(i_file, '(A)') buf ! bandnames of the calculated filters
          IF(index(buf, '#') == 0) THEN
             param%line_filters = trim(buf)
          ELSE
             param%line_filters = buf(1:index(buf, '#')-1)
          ENDIF
       ENDIF
    ENDIF

    IF(param%nwaveq < 1) THEN
       print '(A, I3, A)', &
            '# The Number of QSO Filters is ', param%nwaveq, '!!!'; stop
    ELSE
       allocate(mag_q(param%nwaveq), stat=ier); call CheckIerr(ier, &
            trim(cerr_a)//' mag_q')
    ENDIF

    close(i_file)

    IF(param%traceIDs == 1) THEN
       i_file = 2; num_tar = 0
       open(i_file,file='inpdata/IDlist.dat', status='old', iostat=ier)
       call CheckIerr(ier, trim(cerr_o)//' inpdata/IDlist.dat')
       read(i_file, '()')

       DO
          read(i_file, '()', end=900)
          num_tar = num_tar + 1
       END DO
       900  close(i_file)
       allocate(targ(num_tar), stat=ier)
       call CheckIerr(ier, 'AllocateTarg: fail to allocate  targ')

       open(i_file,file='inpdata/IDlist.dat', status='old', iostat=ier)
       call CheckIerr(ier, trim(cerr_o)//' inpdata/IDlist.dat')
       read(i_file, '()')

       DO i = 1, num_tar 
          read(i_file, *) &
             targ(i)%fdes, targ(i)%iforest, targ(i)%hstart 
       END DO
       close(i_file)
    ENDIF
  END SUBROUTINE ReadParameters
!!$============================================================================
  SUBROUTINE ReadMCMCParameters(nop)
    INTEGER, INTENT(INOUT) :: nop
    INTEGER :: i, ier, i_file
    CHARACTER(LEN=500) :: buf
    CHARACTER(LEN=50) :: cerr_o = '# ReadMCMCParameters: fail to open file'
    CHARACTER(LEN=50) :: cerr_a = '# ReadMCMCParameters: fail to allocate'


    param%munit = 1.d+14
    param%log10munit = log10(param%munit)

    i_file = 1
    open(i_file, file='inpdata/'//trim(file_i)//'_MCMCparams.dat', status='old', iostat=ier)
    call CheckIerr(ier, trim(cerr_o)//' inpdata/'//trim(file_i)//&
         '_MCMCparams.dat')

    read(i_file, '()') 
    read(i_file, '()')  
    read(i_file, '()')  
    read(i_file, '()')  
    read(i_file, '()')  

    read(i_file, *) mcmc%iMCMC_max ! maximum No. of mcmc loop

    read(i_file, '()')  
    read(i_file, '()')  

    read(i_file, *) mcmc%alpst;     IF(mcmc%alpst) nop = nop + 1  
    read(i_file, *) mcmc%alpst_min  
    read(i_file, *) mcmc%alpst_max
    IF(param%alpst < mcmc%alpst_min .or. mcmc%alpst_max < param%alpst) THEN
       print *, "ReadMCMCParameters: value error at alpst", &
                 param%alpst, "[min, max] = ", mcmc%alpst_min, mcmc%alpst_max
       stop
    ENDIF
    read(i_file, *) mcmc%Vst;         IF(mcmc%Vst) nop = nop + 1  
    read(i_file, *) mcmc%Vst_min
    read(i_file, *) mcmc%Vst_max
    IF(param%Vst < mcmc%Vst_min .or. mcmc%Vst_max < param%Vst) THEN
       print *, "ReadMCMCParameters: value error at Vst", &
                 param%Vst, "[min, max] = ", mcmc%Vst_min, mcmc%Vst_max
       stop
    ENDIF
    read(i_file, *) mcmc%tau0st;   IF(mcmc%tau0st) nop = nop + 1  
    read(i_file, *) mcmc%tau0st_min 
    read(i_file, *) mcmc%tau0st_max
    IF(param%tau0st < mcmc%tau0st_min .or. mcmc%tau0st_max < param%tau0st) THEN
       print *, "ReadMCMCParameters: value error at tau0st", &
                 param%tau0st, "[min, max] = ", mcmc%tau0st_min, mcmc%tau0st_max
       stop
    ENDIF
    read(i_file, *) mcmc%eps_SF;   IF(mcmc%eps_SF) nop = nop + 1  
    read(i_file, *) mcmc%eps_SF_min
    read(i_file, *) mcmc%eps_SF_max
    IF(param%eps_SF < mcmc%eps_SF_min .or. mcmc%eps_SF_max < param%eps_SF) THEN
       print *, "ReadMCMCParameters: value error at eps_SF", &
                 param%eps_SF, "[min, max] = ", mcmc%eps_SF_min, mcmc%eps_SF_max
       stop
    ENDIF
    read(i_file, *) mcmc%emax;       IF(mcmc%emax) nop = nop + 1  
    read(i_file, *) mcmc%emax_min
    read(i_file, *) mcmc%emax_max
    IF(param%emax < mcmc%emax_min .or. mcmc%emax_max < param%emax) THEN
       print *, "ReadMCMCParameters: value error at emax", &
                 param%emax, "[min, max] = ", mcmc%emax_min, mcmc%emax_max
       stop
    ENDIF
    read(i_file, *) mcmc%emin;       IF(mcmc%emin) nop = nop + 1  
    read(i_file, *) mcmc%emin_min
    read(i_file, *) mcmc%emin_max
    IF(log10(param%emin) < mcmc%emin_min .or. mcmc%emin_max < log10(param%emin)) THEN
       print *, "ReadMCMCParameters: value error at emin", &
                 log10(param%emin), "[min, max] = ", &
                 mcmc%emin_min, mcmc%emin_max
       stop
    ENDIF
    read(i_file, *) mcmc%Zch;         IF(mcmc%Zch) nop = nop + 1  
    read(i_file, *) mcmc%Zch_min
    read(i_file, *) mcmc%Zch_max
    IF(param%Zch < mcmc%Zch_min .or. mcmc%Zch_max < param%Zch) THEN
       print *, "ReadMCMCParameters: value error at Zch", &
                 param%Zch, "[min, max] = ", mcmc%Zch_min, mcmc%Zch_max
       stop
    ENDIF
    read(i_file, *) mcmc%taud_th; IF(mcmc%taud_th) nop = nop + 1  
    read(i_file, *) mcmc%taud_th_min
    read(i_file, *) mcmc%taud_th_max
    IF(param%taud_th < mcmc%taud_th_min .or. mcmc%taud_th_max < param%taud_th) THEN
       print *, "ReadMCMCParameters: value error at taud_th", &
                 param%taud_th, "[min, max] = ", mcmc%taud_th_min, mcmc%taud_th_max
       stop
    ENDIF

    read(i_file, '()')  
    read(i_file, '()')  

    read(i_file, *) mcmc%Vhot(1);     IF(mcmc%Vhot(1)) nop = nop + 1  
    read(i_file, *) mcmc%Vhot_min(1)
    read(i_file, *) mcmc%Vhot_max(1)
    IF(param%Vhot(1) < mcmc%Vhot_min(1) .or. mcmc%Vhot_max(1) < param%Vhot(1)) THEN
       print *, "ReadMCMCParameters: value error at Vhot(1)", &
                 param%Vhot(1), "[min, max] = ", mcmc%Vhot_min(1), mcmc%Vhot_max(1)
       stop
    ENDIF
    read(i_file, *) mcmc%Vhot(2);     IF(mcmc%Vhot(2)) nop = nop + 1  
    read(i_file, *) mcmc%Vhot_min(2)
    read(i_file, *) mcmc%Vhot_max(2)
    IF(param%Vhot(2) < mcmc%Vhot_min(2) .or. mcmc%Vhot_max(2) < param%Vhot(2)) THEN
       print *, "ReadMCMCParameters: value error at Vhot(2)", &
                 param%Vhot(2), "[min, max] = ", mcmc%Vhot_min(2), mcmc%Vhot_max(2)
       stop
    ENDIF
    read(i_file, *) mcmc%alphot(1); IF(mcmc%alphot(1)) nop = nop + 1  
    read(i_file, *) mcmc%alphot_min(1)
    read(i_file, *) mcmc%alphot_max(1)
    IF(param%alphot(1) < mcmc%alphot_min(1) .or. mcmc%alphot_max(1) < param%alphot(2)) THEN
       print *, "ReadMCMCParameters: value error at alphot(1)", &
                 param%alphot(1), "[min, max] = ", mcmc%alphot_min(1), mcmc%alphot_max(1)
       stop
    ENDIF
    read(i_file, *) mcmc%alphot(2); IF(mcmc%alphot(2)) nop = nop + 1  
    read(i_file, *) mcmc%alphot_min(2)
    read(i_file, *) mcmc%alphot_max(2)
    IF(param%alphot(2) < mcmc%alphot_min(2) .or. mcmc%alphot_max(2) < param%alphot(2)) THEN
       print *, "ReadMCMCParameters: value error at alphot(2)", &
                 param%alphot(2), "[min, max] = ", mcmc%alphot_min(2), mcmc%alphot_max(2)
       stop
    ENDIF
    read(i_file, *) mcmc%alp_ret;     IF(mcmc%alp_ret) nop = nop + 1  
    read(i_file, *) mcmc%alp_ret_min
    read(i_file, *) mcmc%alp_ret_max
    IF(param%alp_ret < mcmc%alp_ret_min .or. mcmc%alp_ret_max < param%alp_ret) THEN
       print *, "ReadMCMCParameters: value error at alp_ret", &
                 param%alp_ret, "[min, max] = ", mcmc%alp_ret_min, mcmc%alp_ret_max
       stop
    ENDIF

    read(i_file, '()')  
    read(i_file, '()')  

    read(i_file, *) mcmc%fmerge;         IF(mcmc%fmerge) nop = nop + 1
    read(i_file, *) mcmc%fmerge_min
    read(i_file, *) mcmc%fmerge_max
    IF(param%fmerge < mcmc%fmerge_min .or. mcmc%fmerge_max < param%fmerge) THEN
       print *, "ReadMCMCParameters: value error at fmerge", &
                 param%fmerge, "[min, max] = ", mcmc%fmerge_min, mcmc%fmerge_max
       stop
    ENDIF
    read(i_file, *) mcmc%fmajor;         IF(mcmc%fmajor) nop = nop + 1
    read(i_file, *) mcmc%fmajor_min
    read(i_file, *) mcmc%fmajor_max
    IF(param%fmajor < mcmc%fmajor_min .or. mcmc%fmajor_max < param%fmajor) THEN
       print *, "ReadMCMCParameters: value error at fmajor", &
                 param%fmajor, "[min, max] = ", mcmc%fmajor_min, mcmc%fmajor_max
       stop
    ENDIF
    read(i_file, *) mcmc%Krem;             IF(mcmc%Krem) nop = nop + 1
    read(i_file, *) mcmc%Krem_min
    read(i_file, *) mcmc%Krem_max
    IF(param%Krem < mcmc%Krem_min .or. mcmc%Krem_max < param%Krem) THEN
       print *, "ReadMCMCParameters: value error at Krem", &
                 param%Krem, "[min, max] = ", mcmc%Krem_min, mcmc%Krem_max
       stop
    ENDIF
    read(i_file, *) mcmc%fdiss;           IF(mcmc%fdiss) nop = nop + 1
    read(i_file, *) mcmc%fdiss_min
    read(i_file, *) mcmc%fdiss_max
    IF(param%fdiss < mcmc%fdiss_min .or. mcmc%fdiss_max < param%fdiss) THEN
       print *, "ReadMCMCParameters: value error at fdiss", &
                 param%fdiss, "[min, max] = ", mcmc%fdiss_min, mcmc%fdiss_max
       stop
    ENDIF
    read(i_file, *) mcmc%Mh0;               IF(mcmc%Mh0) nop = nop + 1
    read(i_file, *) mcmc%Mh0_min
    read(i_file, *) mcmc%Mh0_max
    IF(param%Mh0 < mcmc%Mh0_min .or. mcmc%Mh0_max < param%Mh0) THEN
       print *, "ReadMCMCParameters: value error at Mh0", &
                 param%Mh0, "[min, max] = ", mcmc%Mh0_min, mcmc%Mh0_max
       stop
    ENDIF
    read(i_file, *) mcmc%alpha_halo; IF(mcmc%alpha_halo) nop = nop + 1
    read(i_file, *) mcmc%alpha_halo_min
    read(i_file, *) mcmc%alpha_halo_max
    IF(param%alpha_halo < mcmc%alpha_halo_min .or. mcmc%alpha_halo_max < param%alpha_halo) THEN
       print *, "ReadMCMCParameters: value error at alpha_halo", &
                 param%alpha_halo, "[min, max] = ", mcmc%alpha_halo_min, mcmc%alpha_halo_max
       stop
    ENDIF

    read(i_file, '()')  
    read(i_file, '()')  

    read(i_file, *) mcmc%Vcut;                   IF(mcmc%Vcut) nop = nop + 1
    read(i_file, *) mcmc%Vcut_min
    read(i_file, *) mcmc%Vcut_max
    IF(log10(param%Vcut) < mcmc%Vcut_min .or. mcmc%Vcut_max < log10(param%Vcut)) THEN
       print *, "ReadMCMCParameters: value error at Vcut", &
                 log10(param%Vcut), "[min, max] = ", mcmc%Vcut_min, mcmc%Vcut_max
       stop
    ENDIF
    read(i_file, *) mcmc%kappa_croton;   IF(mcmc%kappa_croton) nop = nop + 1
    read(i_file, *) mcmc%kappa_croton_min
    read(i_file, *) mcmc%kappa_croton_max
    IF(param%kappa_croton < mcmc%kappa_croton_min .or. mcmc%kappa_croton_max < param%kappa_croton) THEN
       print *, "ReadMCMCParameters: value error at kappa_croton", &
                 param%kappa_croton, "[min, max] = ", mcmc%kappa_croton_min, mcmc%kappa_croton_max
       stop
    ENDIF
    read(i_file, *) mcmc%eta_croton;       IF(mcmc%eta_croton) nop = nop + 1
    read(i_file, *) mcmc%eta_croton_min
    read(i_file, *) mcmc%eta_croton_max
    IF(param%eta_croton < mcmc%eta_croton_min .or. mcmc%eta_croton_max < param%eta_croton) THEN
       print *, "ReadMCMCParameters: value error at eta_croton", &
                 param%eta_croton, "[min, max] = ", mcmc%eta_croton_min, mcmc%eta_croton_max
       stop
    ENDIF
    read(i_file, *) mcmc%alp_bower;         IF(mcmc%alp_bower) nop = nop + 1
    read(i_file, *) mcmc%alp_bower_min
    read(i_file, *) mcmc%alp_bower_max
    IF(param%alp_bower < mcmc%alp_bower_min .or. mcmc%alp_bower_max < param%alp_bower) THEN
       print *, "ReadMCMCParameters: value error at alp_bower", &
                 param%alp_bower, "[min, max] = ", mcmc%alp_bower_min, mcmc%alp_bower_max
       stop
    ENDIF
    read(i_file, *) mcmc%eps_bower; IF(mcmc%eps_bower) nop = nop + 1
    read(i_file, *) mcmc%eps_bower_min
    read(i_file, *) mcmc%eps_bower_max
    IF(log10(param%eps_bower) < mcmc%eps_bower_min .or. mcmc%eps_bower_max < log10(param%eps_bower)) THEN
       print *, "ReadMCMCParameters: value error at eps_bower", &
                 log10(param%eps_bower), "[min, max] = ", &
                 mcmc%eps_bower_min, mcmc%eps_bower_max
       stop
    ENDIF

    read(i_file, '()')  
    read(i_file, '()')  

    read(i_file, *) mcmc%em;   IF(mcmc%em) nop = nop + 1
    read(i_file, *) mcmc%em_min
    read(i_file, *) mcmc%em_max
    IF(param%em < mcmc%em_min .or. mcmc%em_max < param%em) THEN
       print *, "ReadMCMCParameters: value error at em", &
                 param%em, "[min, max] = ", mcmc%em_min, mcmc%em_max
       stop
    ENDIF
    read(i_file, *) mcmc%fdi; IF(mcmc%fdi) nop = nop + 1
    read(i_file, *) mcmc%fdi_min
    read(i_file, *) mcmc%fdi_max
    IF(param%fdi < mcmc%fdi_min .or. mcmc%fdi_max < param%fdi) THEN
       print *, "ReadMCMCParameters: value error at fdi", &
                 param%fdi, "[min, max] = ", mcmc%fdi_min, mcmc%fdi_max
       stop
    ENDIF

    read(i_file, '()')  
    read(i_file, '()')  

    read(i_file, *) mcmc%fbh(1);     IF(mcmc%fbh(1)) nop = nop + 1
    read(i_file, *) mcmc%fbh_min(1)
    read(i_file, *) mcmc%fbh_max(1)
    IF(param%fbh(1) < mcmc%fbh_min(1) .or. mcmc%fbh_max(1) < param%fbh(1)) THEN
       print *, "ReadMCMCParameters: value error at fbh(1)", &
                 param%fbh(1), "[min, max] = ", mcmc%fbh_min(1), mcmc%fbh_max(1)
       stop
    ENDIF
    read(i_file, *) mcmc%fbh(2);     IF(mcmc%fbh(2)) nop = nop + 1
    read(i_file, *) mcmc%fbh_min(2)
    read(i_file, *) mcmc%fbh_max(2)
    IF(param%fbh(2) < mcmc%fbh_min(2) .or. mcmc%fbh_max(2) < param%fbh(2)) THEN
       print *, "ReadMCMCParameters: value error at fbh(2)", &
                 param%fbh(2), "[min, max] = ", mcmc%fbh_min(2), mcmc%fbh_max(2)
       stop
    ENDIF
    read(i_file, *) mcmc%Mbhseed;   IF(mcmc%Mbhseed) nop = nop + 1
    read(i_file, *) mcmc%Mbhseed_min
    read(i_file, *) mcmc%Mbhseed_max
    IF(log10(param%Mbhseed) < mcmc%Mbhseed_min .or. mcmc%Mbhseed_max < log10(param%Mbhseed)) THEN
       print *, "ReadMCMCParameters: value error at Mbhseed", &
                 param%Mbhseed, "[min, max] = ", mcmc%Mbhseed_min, mcmc%Mbhseed_max
       stop
    ENDIF
    read(i_file, *) mcmc%tacc_0;     IF(mcmc%tacc_0) nop = nop + 1
    read(i_file, *) mcmc%tacc_0_min
    read(i_file, *) mcmc%tacc_0_max
    IF(param%tacc_0 < mcmc%tacc_0_min .or. mcmc%tacc_0_max < param%tacc_0) THEN
       print *, "ReadMCMCParameters: value error at tacc_0", &
                 param%tacc_0, "[min, max] = ", mcmc%tacc_0_min, mcmc%tacc_0_max
       stop
    ENDIF
    read(i_file, *) mcmc%tacc_1;     IF(mcmc%tacc_1) nop = nop + 1
    read(i_file, *) mcmc%tacc_1_min
    read(i_file, *) mcmc%tacc_1_max
    IF(param%tacc_1 < mcmc%tacc_1_min .or. mcmc%tacc_1_max < param%tacc_1) THEN
       print *, "ReadMCMCParameters: value error at tacc_1", &
                 param%tacc_1, "[min, max] = ", mcmc%tacc_1_min, mcmc%tacc_1_max
       stop
    ENDIF
    read(i_file, *) mcmc%alpha_ad; IF(mcmc%alpha_ad) nop = nop + 1
    read(i_file, *) mcmc%alpha_ad_min
    read(i_file, *) mcmc%alpha_ad_max
    IF(param%alpha_ad < mcmc%alpha_ad_min .or. mcmc%alpha_ad_max < param%alpha_ad) THEN
       print *, "ReadMCMCParameters: value error at alpha_ad", &
                 param%alpha_ad, "[min, max] = ", mcmc%alpha_ad_min, mcmc%alpha_ad_max
       stop
    ENDIF
    read(i_file, *) mcmc%gamma_ad; IF(mcmc%gamma_ad) nop = nop + 1
    read(i_file, *) mcmc%gamma_ad_min
    read(i_file, *) mcmc%gamma_ad_max
    IF(param%gamma_ad < mcmc%gamma_ad_min .or. mcmc%gamma_ad_max < param%gamma_ad) THEN
       print *, "ReadMCMCParameters: value error at gamma_ad", &
                 param%gamma_ad, "[min, max] = ", mcmc%gamma_ad_min, mcmc%gamma_ad_max
       stop
    ENDIF
    read(i_file, *) mcmc%Edd_max;   IF(mcmc%Edd_max) nop = nop + 1
    read(i_file, *) mcmc%Edd_max_min
    read(i_file, *) mcmc%Edd_max_max
    IF(param%Edd_max < mcmc%Edd_max_min .or. mcmc%Edd_max_max < param%Edd_max) THEN
       print *, "ReadMCMCParameters: value error at Edd_max", &
                 param%Edd_max, "[min, max] = ", mcmc%Edd_max_min, mcmc%Edd_max_max
       stop
    ENDIF

    close(i_file)
  END SUBROUTINE ReadMCMCParameters
!!$============================================================================
  SUBROUTINE SubstFilterNumb
    INTEGER :: i, j, ispace, icheck
    CHARACTER(LEN=500) :: charac
    CHARACTER(LEN=5)   :: c_rest = '_rest'
    CHARACTER(LEN=20)  :: bandname
    CHARACTER(LEN=20), PARAMETER :: bandname_base(NWAVE_ALL) = (/CHARACTER(LEN=20)::&
         'B', 'U', 'V', 'Rc', 'RJ', 'Ic', 'IJ', 'z', 'J', 'H', 'K', 'Kp', 'L', & ! 1--13
         'Suprime_B', 'Suprime_g', 'Suprime_V', 'Suprime_r', 'Suprime_R', &
           'Suprime_ip', 'Suprime_I', 'Suprime_zp', 'Suprime_zp_red', & ! 14--22
         'obs5um', & ! 23
         '2MASS_J', '2MASS_H', '2MASS_Ks', & ! 24--26
         'ACS_F775W', 'ACS_F850LP', & ! 27--28
         'CISCO_z', 'CISCO_J', 'CISCO_H', 'CISCO_K', 'CISCO_Kp', & ! 29--33
         'GALEX_FUV', 'GALEX_NUV', & ! 34--35
         'GOODS_RB', 'GOODS_RG', 'GOODS_RS', & ! 36--38
         'HST_F300w', 'HST_F450w', 'HST_F555w', 'HST_F606w', 'HST_F702w', &
           'HST_F814w', & ! 39--44
         'SDSS_up', 'SDSS_gp', 'SDSS_rp', 'SDSS_ip', 'SDSS_zp', & ! 45--49
         'NLyC', 'L1216', 'L1400', 'L1500', 'L1600', 'L1700', 'L2800', &
           'L4861', 'L6563', & ! 50--58
         'IRACch1', 'IRACch2', 'IRACch3', 'IRACch4', & ! 59--62
         'WFCAM_z', 'WFCAM_Y', 'WFCAM_J', 'WFCAM_H', 'WFCAM_K', & ! 63--67
         'VIRCAM_Y', 'VIRCAM_J', 'VIRCAM_H', 'VIRCAM_K', & ! 68--71
         'CIBER_I', 'CIBER_H', & ! 72--73
         'NEWFIRM_J1', 'NEWFIRM_J2', 'NEWFIRM_J3', 'NEWFIRM_H1', 'NEWFIRM_H2', &
           'NEWFIRM_Ks', & ! 74--79
         'AKARI_N2', 'AKARI_N3', 'AKARI_N4', 'AKARI_S7', 'AKARI_S9W', &
           'AKARI_S11', & ! 80--85
         'WISH_0', 'WISH_1', 'WISH_2', 'WISH_3', 'WISH_4', 'WISH_5', & ! 86--91
         'ACS_F435W', 'ACS_F475W', 'ACS_F606W', 'ACS_F814W', 'NICMOS_F160W', &
           'WFPC2_F300W', 'WFPC2_F450W', 'WFPC2_F555W', 'WFC3_F125W', &
           'WFC3_F140W', 'WFC3_F160W', & ! 92--102
         'HSC_g', 'HSC_r', 'HSC_i', 'HSC_z', 'HSC_Y' & ! 103--107
         /)


    DO i = 1, NWAVE_ALL
       j = i + NWAVE_ALL
       ssp%bandname(i) = trim(bandname_base(i))
       ssp%bandname(j) = trim(bandname_base(i))//c_rest
    ENDDO

    call InitializeGlobalFilterIntegers

    IF(param%nwave < TNWAVE_ALL) THEN
       ! --- check the number of given bandnames against 'param%nwave'
       charac = param%line_filters; i = 0
       DO WHILE(index(charac, ' ') /= 0 .and. len_trim(charac) > 0)
          ispace = index(charac, ' ')
          charac = charac(ispace+1:len_trim(charac))
          i = i + 1
       ENDDO
       IF(i /= param%nwave) THEN
          print '(A, 2I3)', &
               '# The Number of Characters in The Final Row of The Input '//&
               'File is Inappropriate!!!: ', i, param%nwave; stop
       ENDIF

       ! --- check the 1st given bandname against 'Br'
       charac = param%line_filters
       ispace = index(charac, ' ')
       bandname = charac(1:ispace-1)
!!$       print *, ispace, charac(1:ispace-1), bandname
       IF(len_trim(bandname) < 2 .or. bandname(1:1) /= 'B' &
            .or. bandname(2:2) /= 'r') THEN
          print "(A)", "# The 1st Charactor in The Final Row of The Input "//&
               "File Should Be 'Br'!!"; stop
       ENDIF
       param%iBband_r = 1           ! mag(1) = B-band magnitude in rest-frame
       param%iwave(1) = 1+NWAVE_ALL ! B-band in rest-frame

       ! --- check the 2nd given bandname against 'Vr'
       charac   = charac(ispace+1:len_trim(charac))
       ispace   = index(charac, ' ')
       bandname = charac(1:ispace-1)
!!$       print *, ispace, charac(1:ispace-1), bandname
       IF(len_trim(bandname) < 2 .or. bandname(1:1) /= 'V' &
            .or. bandname(2:2) /= 'r') THEN
          print "(A)", "# The 2nd Charactor in The Final Row of The Input "//&
               "File Should Be 'Vr'!!"; stop
       ENDIF
       param%iVband_r = 2           ! mag(2) = V-band magnitude in rest-frame
       param%iwave(2) = 3+NWAVE_ALL ! V-band in rest-frame

       ! --- read the other bandnames which are used in this run
       DO i = 3, param%nwave
          ! --- substitute the corresponding integer to 'param%iwave(i)'
          param%iwave(i) = 0; icheck = 0
          charac   = charac(ispace+1:len_trim(charac))
          ispace   = index(charac, ' ')
          bandname = charac(1:ispace-1)
          j = len_trim(bandname)
          IF(bandname(j:j) == 'r' .and. &
               (bandname /= 'Suprime_r' .and. bandname /= 'HSC_r')) THEN
             bandname = bandname(1:j-1)//c_rest
          ENDIF
          DO j = 1, TNWAVE_ALL
             IF(len_trim(bandname) == len_trim(ssp%bandname(j))) THEN
                IF(bandname == ssp%bandname(j)) THEN
                   param%iwave(i) = j; icheck = 1
!!$                   j = NWAVE_ALL
                   exit
                ENDIF
             ENDIF
          ENDDO
          IF(icheck == 0) THEN
             print "(A)", "# There is NO Prepared Bandname in sspfiles "//&
                  "corresponding to '"//trim(bandname)//"'!!"
             print '(A)', "# The bandname should be selected in the "//&
                  "following names"
             print '(A, $)', "# "
             DO icheck = 1, NWAVE_ALL
                IF(inode == 0) print '(I3, A)', icheck, ': '//&
                     trim(ssp%bandname(icheck))//' (obs-frame)'//&
                     trim(ssp%bandname(icheck))//'_r (rest-frame)'
             ENDDO
             print *; stop
          ENDIF
          IF(inode == 0) print *, i, param%iwave(i), bandname, &
               ssp%bandname(param%iwave(i))
          call SubstGlobalFilterIntegers(i)
       ENDDO
    ELSE ! the case of param%nwave = TNWAVE_ALL
       call SubstGlobalFilterIntegersForAll
    ENDIF
  END SUBROUTINE SubstFilterNumb
!!$============================================================================
  SUBROUTINE SetConstants
    DOUBLE PRECISION, PARAMETER :: MB_sol = +5.33d0 ! absolute magnitude of the
                                                    !  Sun in B-band [ABmag]
                                   ! http://www.ucolick.org/~cnaw/sun.html
    INTEGER :: ier, ilam
    CHARACTER(LEN=50) :: cerr_a = '# SetConstants: fail to allocate'
    CHARACTER(LEN=50) :: cerr_d = '# SetConstants: fail to deallocate'
    DOUBLE PRECISION :: f21, Lsol_B
    DOUBLE PRECISION, ALLOCATABLE :: const_delc(:)
    DOUBLE PRECISION :: Square, z2t ! functions
    ! --- the followings are added to calculate const%corr_wdim(:) by MARK (2018/Mar/02)
    DOUBLE PRECISION, PARAMETER :: hP     = 6.6260693d-27 ! Planck constant [ergs s]
    DOUBLE PRECISION, PARAMETER :: LightV = 2.99792458d+18 ! c [A/s]
    DOUBLE PRECISION, PARAMETER :: lam(NWAVE_wdim) = (/912.d0, 1216.d0, 1400.d0, &
         1500.d0, 1600.d0, 1700.d0, 2800.d0, 4861.d0, 6563.d0/)
    DOUBLE PRECISION :: tmp


    const%PI   = acos(-1.d0)
    const%PI2  = Square(const%PI)
    const%Zsun = 0.019d0 ! Solar metallicity
    const%corr = -5.d0 * log10(param%h)
    const%corr_wdim(4) = 81.90637d0
                         ! -2.5log[(lam^2/c)/4pi/(10pc)^2]-48.6 at lam = 1500 A
    DO ilam = 1, NWAVE_wdim
       IF(ilam == 1) THEN
          tmp = hP * LightV / lam(ilam) ! = hnu [erg] per photon
          tmp = tmp / lam(ilam) ! [erg/A]
          const%corr_wdim(ilam) = -2.5d0 * log10(tmp) + const%corr_wdim(4) &
                                  - 5.d0 * log10(lam(ilam) / lam(4))
       ELSEIF(ilam /= 4) THEN
          const%corr_wdim(ilam) = const%corr_wdim(4) - 5.d0 * log10(lam(ilam) / lam(4))
       ENDIF
    ENDDO

    const%nub = 6.82d+14   ! [Hz] for bolometric correction of Marconi+ 04.
    const%B2UV = 0.87      ! UV-mag (1450\AA) = B-mag + const%B2UV (for AGNs)
                           ! (Added by Shirakata 2018/Mar/05)

    ! for M_HI/L_B
    Lsol_B = 10.d0**(-0.4d0 * MB_sol) ! normalized B-band luminosity of the Sun
    !const%fracH = 0.75d0 ! Hydrogen mass fraction in cold gas mass
    const%fracH = 0.54d0 ! atomic Hydrogen mass fraction in cold gas mass
    const%ML = const%fracH * param%munit * Lsol_B
               ! conversion factor from Mcool/lum_B [10^14 Msun/3630 Jy]
               !  to [Msun/Lsun^B]

    ! for the function of "delc"
    allocate(const_delc(6), stat=ier); call CheckIerr(ier, &
         trim(cerr_a)//' const_delc')
    const_delc(1) = 1.00029d0;     const_delc(2) = 0.169218d0
    const_delc(3) = -0.0280635d0;  const_delc(4) = 0.00400883d0
    const_delc(5) = -0.00031472d0; const_delc(6) = 9.88907d-6
    f21 = 1.d0 / (const_delc(1) + (const_delc(2) + (const_delc(3) + &
         (const_delc(4) + (const_delc(5) + const_delc(6) * param%OMEGA_rat) &
         * param%OMEGA_rat) * param%OMEGA_rat) * param%OMEGA_rat) &
         * param%OMEGA_rat)
    const%delc0 = 3.d0 * (12.d0*const%PI)**(2.d0/3.d0) * f21 / 20.d0
    deallocate(const_delc, stat=ier); call CheckIerr(ier, &
         trim(cerr_d)//' const_delc')

    const%V1 = 1.d0 / 3.d0; const%V2 = 1.d0 / 6.d0 ! for the function of "CircVel"
    const%Rb = 2.175d+8 * param%h ! [kpc/h] for the function of "CalRbulge"

    const%tout    = z2t(param%zsp1) ! the age of the universe at output redshift
                                    !   in hubble time
    const%toutGyr = const%tout * param%th ! the age of the universe at output
                                          !   redshift in Gyr

    ! --- for calculation of mean SFR (added by MARK on 2014/Nov/04)
    const%t0 = const%t0yr / param%th_yr ! in hubble time

    IF(param%UVfb == 2) THEN ! Okamoto+08 UV feedback mode
       param%Vlow = 0.d0
       const%alpha_UV = 2.d0
       const%UV1 = 2.d0**(const%alpha_UV / 3.d0) - 1.d0
       const%UV2 = -3.d0 / const%alpha_UV
    ENDIF

    const%t_dyn0 = sqrt(3.d0 * const%PI / (32.d0 * 6.674d-11 * 1.88d-26 &
                                          * Square(param%h))) / const%yr2sec
    ! --- for calculation of the Bower AGN feedback mode
    IF(param%AGNFB_key == 3) THEN
       const%dMbh0  = 1.d-7 / (Square(const%c) * 1.d-1 * const%Msolar &
            * param%munit)
    ENDIF
  END SUBROUTINE SetConstants
!!$============================================================================
  SUBROUTINE WriteCaptions(ionum, nz,nz_end)
    INTEGER, INTENT(IN) :: ionum, nz, nz_end
    INTEGER :: i, i0, i1, i2, i3, i4, i5, i6, i7, i8, i9, i10, i11, i12, iout
    DOUBLE PRECISION :: redshift
    CHARACTER(LEN=20) :: type_mag(2) = (/CHARACTER(LEN=20) ::'[Vega mag]', '[AB mag]'/)


    redshift = param%zsp1 - 1.d0

    ! +---------------------+
    ! | Captions for *0.dat |
    ! +---------------------+
    IF(param%run_type /= 2) THEN
      write(ionum, '(A)') '# log file: the adopted parameter values in '//&
           'the calculation are written'
    ELSE
      IF(nz == 1) THEN
        write(ionum, '(A)') '# log file: the adopted parameter values in '//&
             'the calculation are written'
      ENDIF
    ENDIF



    ! +---------------------+
    ! | Captions for *1.dat |
    ! +---------------------+
    IF(param%run_type == 2) THEN ! all redshift
       iout = ionum + NFILE * nz_end + nnode * (nz-1) + 1
    ELSE ! other
       iout = ionum + NFILE + 1
    ENDIF

    write(iout, '(A, F8.5)') '# properties of each galaxy at z = ', redshift
    ! --- the following captions should be modified according to
    !       the expressions around LL166--L187 in "out_norm_l9a.f90"
    write(iout, '(A, F6.2, A, $)') '# (1)ID (2)ID of central galaxy '//&
         '(3,4)mor,mor_d[1:E,2:S0,3:S] (5)flag_c (6)flag_burst '//&
         '(7)flag_merger (8)flag_di (9)flag_ccut (10)iforest '//&
         '(11)ID(host) (12)mpi(gal) (13)mpi(host) (14)f_des (15)BT ratio '//&
         'in rest B-band (16,17,18,19)Mstarb,Mstard,Mcoolb,Mcoold[Msun] '//&
         '(20,21)Mhalo^{host,prog}[Msun] (22)redshift of first collapse '//&
         '(23,24)Vcirc^{host,prog}[km/s] (25,26)MZcb,MZcd[Msun] '//&
         '(27,28)z{g,t} (29)rb0[kpc/h] (30)Vbulge[km/s] (31)rd0[kpc/h] '//&
         '(32)Vdisk[km/s] (33)SB(Brest)[mag/arcsec^2] (34)m/L(Brest) '//&
         '(35)Mdyn[Msun] (36)Mbh[Msun] (37)MZc_rem[Msun] (38)SFR[Msun/yr] '//&
         '(39)mean SFR during the past ', const%t0yr * 1.d-6, ' Myr[Msun/yr] '//&
         '(40,41,42)mass-weighted mean stellar metallicity(bulge,disk,all: '//&
         '-99 for no stars) (43)bolometric lumq[erg/s] (44)taust[yr] (45)Lpeak[10^12 Lsun] '//&
         '(46)dMstar in burst[Msun] '
    i0 = 47
    IF(param%iNLyC /= 0 .or. param%iNLyC_r /= 0) THEN
       write(iout, '(A, I2, A, $)') '(', i0, ')L(Lya)(w/o dust)[erg/s] '
       i0 = i0 + 1
    ENDIF
    DO iband = 1, param%nwave
       i1 = i0 + 2 * (iband - 1)
       IF(i1 < 99) THEN
          write(iout, '(2(A, I2), A, $)') '(', i1, ',', i1+1, ')'
       ELSEIF(i1 == 99) THEN
          write(iout, '(A, I2, A, I3, A, $)') '(', i1, ',', i1+1, ')'
       ELSE
          write(iout, '(2(A, I3), A, $)') '(', i1, ',', i1+1, ')'
       ENDIF
       write(iout, '(A, $)') trim(ssp%bandname(param%iwave(iband)))
       IF(param%iwave(iband) == 50 .or. param%iwave(iband) == NWAVE_ALL+50) THEN
          ! NLyC
          write(iout, '(A, $)') '(w/o,w/ dust)[photons/s] '
       ELSEIF((param%iwave(iband) >= 51 .and. param%iwave(iband) <= 58) &
            .or. (param%iwave(iband) >= NWAVE_ALL+51 .and. &
                  param%iwave(iband) <= NWAVE_ALL+58)) THEN
          ! L1216, L1400, L1500, L1600, L1700, L2800, L4861, L6563
          write(iout, '(A, $)') '(w/o,w/ dust)[erg/s/A] '
       ELSE
          write(iout, '(A, $)') '-band mag(w/o,w/ dust)'//&
               trim(type_mag(param%type_mag))//' '
       ENDIF
    ENDDO
    write(iout, *)

    ! +---------------------+
    ! | Captions for *2.dat |
    ! +---------------------+
    iout = ionum + 2
    IF(param%run_type == 2) &
       iout = ionum + NFILE * (nz-1) + 2

    write(iout, '(A, F8.5)') '# LFs in each band at z = ', redshift
    i1 = param%nwave + 1
    i2 = i1 + param%nwave
    i3 = i2 + param%nwave
    i4 = i3 + param%nwave
    write(iout, '(7(A, I3), A)') '# (1)M-5log(h) '//&
         trim(type_mag(param%type_mag))//' (2-', i1, ', ', i1+1, '-', i2, &
         ')dn/dM[h^3/Mpc^3/mag] for burst+quiescent (w/o, w/ dust), [space] (', &
         i2+1, '-', i3, ',', i3+1, '-', i4, ')dn/dM[h^3/Mpc^3/mag] for '// &
         'burst (w/o, w/ dust)'
    write(iout, '(A, $)') '#'
    DO i = 1, param%nwave
       write(iout, '(4(A, I3), A, $)') ' (', i+1, ',', i+i1, ':', i+i2, ',', &
            i+i3, ')'//trim(ssp%bandname(param%iwave(i)))
    ENDDO
    write(iout, '(A)') ' (w/o, w/ dust: w/o, w/ dust for burst)'


    ! +---------------------+
    ! | Captions for *3.dat |
    ! +---------------------+
    iout = ionum + 3
    IF(param%run_type == 2) &
       iout = ionum + NFILE * (nz-1) + 3

    write(iout, '(A, F8.5)') &
         '# morphology-dependent LFs in each band at z = ', redshift
    i0 = 3 * param%nwave
    i1 = i0 + 1
    i2 = i1 + i0
    i3 = i2 + i0
    i4 = i3 + i0
    write(iout, '(7(A, I3), A)') '# (1)M-5log(h) '//&
         trim(type_mag(param%type_mag))//' (2-', i1, ', ', i1+1, '-', i2, &
         ')dn/dM[h^3/Mpc^3/mag] for burst+quiescent (w/o, w/ dust), [space] (', &
         i2+1, '-', i3, ',', i3+1, '-', i4, ')dn/dM[h^3/Mpc^3/mag] for '// &
         'burst (w/o, w/ dust)'
    i1  = param%nwave + 1
    i2  = i1  + param%nwave
    i3  = i2  + param%nwave
    i4  = i3  + param%nwave
    i5  = i4  + param%nwave
    i6  = i5  + param%nwave
    i7  = i6  + param%nwave
    i8  = i7  + param%nwave
    i9  = i8  + param%nwave
    i10 = i9  + param%nwave
    i11 = i10 + param%nwave
    i12 = i11 + param%nwave
    write(iout, '(23(A, I3), A)') &
         '# (2-',        i1, ',', i3+1, '-', i4, ':', i6+1, '-', i7, ',', &
          i9+1, '-', i10, ')E, '//&
         '(', i1+1, '-', i2, ',', i4+1, '-', i5, ':', i7+1, '-', i8, ',', &
         i10+1, '-', i11, ')S0, '//&
         '(', i2+1, '-', i3, ',', i5+1, '-', i6, ':', i8+1, '-', i9, ',', &
         i11+1, '-', i12, ')S '//&
         '(w/o, w/ dust for burst+quiescent: w/o, w/ dust for burst)'
    DO i = 1, param%nwave
       write(iout, '(12(A, I3), A)') '# (', &
            i+ 1, ',',  i+i1, ',',  i+i2, '[w/o dust for burst+quiescent]:', &
            i+i3, ',',  i+i4, ',',  i+i5, '[w/ dust for burst+quiescent]: ', &
            i+i6, ',',  i+i7, ',',  i+i8, '[w/o dust for burst]:', &
            i+i9, ',', i+i10, ',', i+i11, '[w/ dust for burst]) '//&
            trim(ssp%bandname(param%iwave(i)))
    ENDDO


    ! +---------------------+
    ! | Captions for *4.dat |
    ! +---------------------+
    iout = ionum + 4
    IF(param%run_type == 2) &
       iout = ionum + NFILE * (nz-1) + 4

    write(iout, '(A, F8.5)') '# Mass Functions at z = ', redshift
    i0 = 11
    i1 = i0 + 1
    i2 = i1 + i0
    write(iout, '(A)') '# index  0: all mass functions'
    DO i = 1, i0
       write(iout, '(A, I2, A)') '# index ', i, ': morphology-dependent '//&
            trim(c_mftype(i))//' function'
    ENDDO
    write(iout, '(A)') '# ------------------------------------------------------'
    write(iout, '(A)') '# index 0: all mass functions'
    write(iout, '(3(A, I2), A)') '# (1)log10[Mass/(Msun/h^2)] '//&
         '(2-', i1, ')dn/dlogM[h^3/Mpc^3/dex] for burst+quiescent, '//&
         '(', i1+1, '-', i2, ')dn/dlogM[h^3/Mpc^3/dex] for burst'
    write(iout, '(11(A,I2), A)') '# (2, ', i1+1, ')Mstar^bulge, (3, ', &
         i1+2, ')Mstar^disk, (4, ', i1+3, ')Mcold, (5, ', i1+4, &
         ')M_SMBH, (6, ', i1+5, ')Mhot, (7, ', i1+6, &
         ')Mstar^bulge+M_SMBH, (8, ', i1+7, ')Mstar^disk+Mcold, (9, ', i1+8, &
         ')Mgalaxy(=Mstar^bulge+M_SMBH+Mstar^disk+Mcold), (10, ', i1+9, &
         ')Mstar(=Mstar^bulge+Mstar^disk) (11, ', i1+10, &
         ')M_H (12, ', i1+11, ')Mhalo (burst+quiescent, burst)'


    ! +---------------------+
    ! | Captions for *5.dat |
    ! +---------------------+
    iout = ionum + 5
    IF(param%run_type == 2) &
       iout = ionum + NFILE * (nz-1) + 5

    write(iout, '(A, F8.5)') '# Hydrogen mass relative to B-band luminosity '//&
         'as a function of M_B-5logh at z = ', redshift
    write(iout, '(A)') '# (1)M_B-5logh'//trim(type_mag(param%type_mag))//&
         ' (2-25)M_H/L_B[Msun/Lsun^B]'
    write(iout, '(A)') '# (2-7)E (8-13)S0 (14-19)S (20-25)all'
    write(iout, '(A)') '# (2,8,14,20)mean (3,9,15,21)mean-sigma '// &
         '(4,10,16,22)mean+sigma'
    write(iout, '(A)') '# (5,11,17,23)median (6,12,18,24)25% percentile '//&
         '(7,13,19,25)75% percentile'


    ! +---------------------+
    ! | Captions for *6.dat |
    ! +---------------------+
    iout = ionum + 6
    IF(param%run_type == 2) &
       iout = ionum + NFILE * (nz-1) + 6

    write(iout, '(A, F8.5)') '# QSO LFs in each band at z = ', redshift
    write(iout, '(A)') '# (1)x (2--5)dn/dx[h^3/Mpc^3/dex]'
    write(iout, '(A)') '# (2)x=M_B-5log(h)[ABmag]'
    write(iout, '(A)') '# (3)x=M_UV-5log(h)[ABmag] (1450 \AA, intrinsic)'
    write(iout, '(A)') '# (4)x=M_UV-5log(h)[ABmag] (1450 \AA, observed)'
    write(iout, '(A)') '# (5)x=log10[L_X(2--10 keV)/(erg/s)]'
    write(iout, '(A)') '# (6)x=log10[L_X(0.5--2 keV)/(erg/s)]'
    write(iout, '(A)') '# (7)x=log10[L_bol/(erg/s)]'
    write(iout, '(A)') '# (8)x=log10[Mdot/Medd] (Lx > 1e+43 erg/s)'
    write(iout, '(A)') '# (9)x=log10[Lbol/Ledd] (Lx > 1e+43 erg/s)'
    write(iout, '(A)') '# (10)x=log10[M_bh/10^14 M_sun] (Lx > 1e+43 erg/s)'

    ! +---------------------+
    ! | Captions for *7.dat |
    ! +---------------------+
    IF(param%run_type == 2) THEN ! all redshift
       iout = ionum + NFILE * nz_end + nnode * (nz_end + nz - 1) + 1 
    ELSE ! other
       iout = ionum + NFILE + nnode + 1
    ENDIF

    write(iout, '(A, F8.5)') '# Quantities related to the AGNs activated '//&
         'at the output redshift zout = ', redshift
    write(iout, '(A)') '# (1)igal (2)flag_burst[1:burst,0:normal] '//&
         '(3)flag_di[1:DI,0:normal] (4)flag_merger[1:merger,0:normal] '//&
         '(5)morphology[1:S,2:S0,3:E] (6)flag_c '//&
         '(7,8)Mstard,Mstarb[Msun] (9,10)Mcoold,Mcoolb[Msun] '//&
         '(11)Mbh[Msun] (12,13)MZcd,MZcb[Msun] '//&
         '(14)Mhalo[Msun] (15,16)Vbulge,Vdisk[km/s] '//&
         '(17,18)rbulge,rdisk[kpc/h] (19,20)mag,mag_d '//&
         '(21,22) SFR, averaged SFR in past 10 Myr [Msun/yr] '//&
         '(23)Lbol[10^12 Lsun] at output time (24)Lbol[10^12 Lsun] of peak '//&
         '(25)Delta Macc [M_sun] (26)Mcold[Msun] at tstart '//&
         '(27)MZcold[Msun] at tstart (28)t_ad[Gyr] (29)t_gal[Gyr] '//&
         '(30)dMbhrem[10^14 Msun] (31,32)merging mass ratio (av/max) '//&
         '(33-36)Lx(2-10keV)[erg/s],Lbol[erg/s],mdot,lambda_edd'

    ! +---------------------+
    ! | Captions for *8.dat |
    ! +---------------------+
    IF(param%traceIDs == 1) THEN
       IF(param%run_type == 2) THEN ! all redshift
          iout = ionum + NFILE * (nz-1) + 8
       ELSE ! other
          iout = ionum + NFILE + nnode + nnode + 1
       ENDIF

       write(iout, '(A, F8.5)') '# Quantities related to the target galaxies '//&
            'at the output redshift zout = ', redshift
       write(iout, '(A)') '# (1)redshift (2)iforest (3)mpi '//&
            '(4)hstart (5)hfinal (6)flag_merger (7)flag_di (8)flag_ccut (9)BT '//&
            '(10-13)Mass[10^14 Msun]; Mcold, Mstard, Mstarb, Mbh (14)Lbolpeak[erg/s] '//&
            '(15)Lbol[erg/s] (16,17)Accretion rate normailzed by Eddington rate at tstart, tout '//&
            '(18,19)Bolometric AGN luminosity normalized by Eddington luminosity at tstart, tout '//&
            '(20)Mbh_init[10^14 Msun] (21)Macc_eff[10^14 Msun]'//&
            '(22,23)tstart, tacc[hubble time]'
      ENDIF
  END SUBROUTINE WriteCaptions
!!$============================================================================
  SUBROUTINE WriteInputParameters(number)
    INTEGER, INTENT(IN) :: number
    INTEGER :: i, ier
    CHARACTER(LEN=100) :: ref(0:NFILE) = &
         (/CHARACTER(LEN=100) :: &
         ' = logfile (THIS FILE)', &
         ' = physical properties of each galaxy', &
         ' = luminosity functions in each band', &
         ' = morphology-dependent luminosity functions in each band', &
         ' = mass functions & morphology-dependent mass functions', &
         ' = hydrogen mass relative to B-band luminosity as a function of M_B-5logh', &
         ' = QSO luminosity functions in each band', &
         ' = AGN condition', &
         ' = ID trace'/)
    CHARACTER(LEN=100) :: cerr = '# WriteInputParameters: fail to'
    CHARACTER(LEN=20)  :: c_runtype(3) = &
         (/CHARACTER(LEN=20)  ::' (a redshift)', ' (all redshifts)', ' (MCMC)'/)
    CHARACTER(LEN=100) :: c_UVfb(2)    = &
         (/CHARACTER(LEN=100) ::' = Vlow', ' = Okamoto+08 Mc(z)'/)
    CHARACTER(LEN=100) :: c_CoolFN(2)    = &
         (/CHARACTER(LEN=100) ::' = SD93', ' = Okamoto based cloudy'/)
    CHARACTER(LEN=100) :: c_dusttype(4) = &
         (/CHARACTER(LEN=100) :: &
           ' = MW-type extinction curve given by Pei (1992)', &
           ' = LMC-type extinction curve given by Pei (1992)', &
           ' = SMC-type extinction curve given by Pei (1992)', &
           ' = Calzetti-law (Calzetti+00)'/)
    CHARACTER(LEN=100) :: c_AGNFB(3) = &
         (/CHARACTER(LEN=100) ::' = Vcut', ' = Croton AGN-FB mode Croton+06)', &
           ' = Bower AGN-FB mode (Bower+06)'/)
    CHARACTER(LEN=100) :: c_SFmodel(2) = (/' = CSF', ' = DSF'/)


    IF(number == 1) THEN
       write(ionum, '(A)') '# +--------------------------------------------+'
       write(ionum, '(A)') '# |      nGC main body related parameters      |'
       write(ionum, '(A)') '# +--------------------------------------------+'
       write(ionum, '(A)') '# N-body run: '//trim(fbase_nbody)//trim(param%file_nbody)
       write(ionum, '(A, I2, A)') '# run type: ', param%run_type, &
            trim(c_runtype(param%run_type))
       write(ionum, '(A)') '# filenames used in the nGC main calculation:'
       write(ionum, '(A, I1, A)') '# --- ', 0, &
            ':'//trim(file_out(0))//trim(ref(0))
       IF(param%run_type /= 3) THEN
          DO i = 1, NFILE
             write(ionum, '(A, I1, A)') '# --- ', i, &
                  ':'//trim(file_out(i))//trim(ref(i))
          ENDDO
       ENDIF

       write(ionum, '(A, I10)')  '# Seed for random number: ', param%idum
       ! write(ionum, '(A, F8.4)') '# Input 1+z: ', param%zsp1
       write(ionum, '(A, F8.4)') '# Input 1+z: ', param%zsp1_input ! modified by MARK
                                                                   !  (2016/Jun/24)
       write(ionum, '(A, F8.4)') '# alpst: ', param%alpst
       write(ionum, '(A, F8.4)') '# Vst: ',   param%Vst
       write(ionum, '(2(A, F8.4), A)') '# Vhot: ', param%Vhot(1), &
            ' (starburst), ', param%Vhot(2), ' (quiescent)'
       write(ionum, '(2(A, F8.4), A)') '# alphot: ', param%alphot(1), &
            ' (starburst), ', param%alphot(2), ' (quiescent)'
       write(ionum, '(A, I2, A)') '# SFmodel: ', param%SFmodel, &
            trim(c_SFmodel(param%SFmodel))
       IF(param%SFmodel == 1) THEN ! CSF
          write(ionum, '(A, F8.4, A)') '# --- tau0st: ', param%tau0st, ' [Gyr]'
       ELSEIF(param%SFmodel == 2) THEN ! DSF
          write(ionum, '(A, F8.4, A)') '# --- eps_SF: ', param%eps_SF
       ELSEIF(param%SFmodel == 3) THEN ! Makiya14
          write(ionum, '(A, F8.4, A)') '# --- emax: ', param%emax, ' [1/Gyr]'
          write(ionum, '(A, F8.4, A)') '# --- emin: ', param%emin, ' [1/Gyr]'
          write(ionum, '(A, F8.4, A)') '# --- Zch: ',  param%Zch,  ' [Z_sun]'
          write(ionum, '(A, F8.4)') '# --- taud_th: ', param%taud_th
       ENDIF


       write(ionum, '(A, F8.4)') '# fmerge: ',  param%fmerge
       write(ionum, '(A, F8.4)') '# fmajor: ',  param%fmajor
       write(ionum, '(A, F8.4)') '# Krem: ',  param%Krem
       write(ionum, '(A, F8.4)') '# alp_ret: ', param%alp_ret

       ! --- logical parameters
       write(ionum, '(A, L)') '# off_collision: ',     param%off_collision
       write(ionum, '(A, L)') '# equal_mass_merger: ', param%equal_mass_merger
       write(ionum, '(A, L)') '# dynamical_response_bulge: ', &
            param%dyn_resp_bulge
       write(ionum, '(A, L)') '# dynamical_response_disk: ', &
            param%dyn_resp_disk
       write(ionum, '(A, L)') &
            '# dynamical_response_halo: ', param%dyn_resp_halo

       write(ionum, '(A, F8.4)') '# freheat: ', param%freheat
       write(ionum, '(A, F8.4)') '# fb: ', param%fb
       write(ionum, '(A, F8.4, A)')  '# Vlow: ', param%Vlow, ' [km/s]'
       write(ionum, '(A, G10.3)') '# tauV0: ',   param%tauV0
                                  !dust, original: 9.d+5
       write(ionum, '(A, I2, A)') '# exttype: ', param%exttype, &
            trim(c_dusttype(param%exttype))
       write(ionum, '(A, F8.4)') '# fdm: ', param%fdm
                                  !for Vcent, DM contribution
       write(ionum, '(A, F8.4)') '# fdiss: ', param%fdiss
                                  !dissipation fraction
       write(ionum, '(A, F8.4)') '# em: ',  param%em
       write(ionum, '(A, F8.4)') '# fdi: ', param%fdi
       write(ionum, '(A, F8.4)') '# alpha_tau: ',   param%alpha_tau !optical depth
       write(ionum, '(A, F8.4)') '# alpha_rad: ',   param%alpha_rad !disk size
       write(ionum, '(A, F8.4)') '# alpha_burst: ', param%alpha_burst
       write(ionum, '(A, F8.4, A)')  '# Mh0: ', param%Mh0, ' [10^14 Msun]'
       write(ionum, '(A, F8.4)') '# alpha_halo: ', param%alpha_halo
       write(ionum, '(A, F8.4)') '# f_ngal: ', param%f_ngal ! added by MARK (2016/Jun/29)
       write(ionum, '(A, I2, A)') '# UVfb: ',  param%UVfb, trim(c_UVfb(param%UVfb))
       IF(param%UVfb == 2) &
            write(ionum, '(A, F7.3)') '# --- zreion: ', param%zreion

       write(ionum, '(A, I2, A)') '# CoolFN: ', param%CoolFN, trim(c_CoolFN(param%CoolFN))
       ! --- SMBH related parameters
       write(ionum, '(A, F8.4)')  '# fbh(1): ',     param%fbh(1)
       write(ionum, '(A, F8.4)')  '# fbh(2): ',     param%fbh(2)
       write(ionum, '(A, G10.3)') '# Mbhseed: ', param%Mbhseed
       write(ionum, '(A, F8.4)')  '# eps_agn: ', param%eps_agn
       write(ionum, '(A, F8.4)') '# tacc_0: ', param%tacc_0
       write(ionum, '(A, F8.4)') '# tacc_1: ', param%tacc_1
       write(ionum, '(A, F8.4)')  '# acc: ',     param%acc
       write(ionum, '(A, F8.4)')  '# alpha_ad: ',     param%alpha_ad
       write(ionum, '(A, F8.4)')  '# gamma_ad: ',     param%gamma_ad
       write(ionum, '(A, F8.4, A)')  '# tad_0: ',     param%tad_0, '[Gyr]'
       write(ionum, '(A, I3)')    '# nwaveq: ',  param%nwaveq

       ! --- AGN feedback related parameters
       write(ionum, '(A, I2, A)') '# AGNFB_key: ', param%AGNFB_key, &
            trim(c_AGNFB(param%AGNFB_key))
       IF(param%AGNFB_key == 1) THEN ! Vcut
          write(ionum, '(A, G10.3, A)') '# --- Vcut: ', param%Vcut, ' [km/s]'
       ELSEIF(param%AGNFB_key == 2) THEN ! Croton
          write(ionum, '(A, G10.3)') '# --- kappa_croton: ', param%kappa_croton
          write(ionum, '(A, G10.3)') '# --- eta_croton: ', param%eta_croton
       ELSEIF(param%AGNFB_key == 3) THEN ! Bower
          write(ionum, '(A, F9.5)') '# --- alp_bower: ', param%alp_bower
          write(ionum, '(A, F9.5)') '# --- eps_bower: ', param%eps_bower
       ENDIF
       ! --- parameters (added by MARK)
       write(ionum, '(A, L)') '# ML: ',  param%ML
       write(ionum, '(A, L)') '# Mbh: ', param%Mbh
       write(ionum, '(A, L)') '# LAE: ', param%LAE
       write(ionum, '(A, L)') '# SFH: ', param%SFH
       write(ionum, '(A, L)') '# SN: ',  param%SN

       ! --- surface brightness limit related parameters (added by MARK, 2015/Apr/11)
       write(ionum, '(A, L)')    '# SB: ', param%SB
       IF(param%SB) THEN
          write(ionum, '(A, F9.5)') '# --- SBlimit: ', param%SBlimit
          write(ionum, '(A)') '# --- used band: '//&
               trim(ssp%bandname(param%iwave(SB%iband)))
       ENDIF

       write(ionum, '(A, L)')  '# DI: ',      param%DI
       write(ionum, '(A, I2)') '# ML_AGN: ',  param%ML_AGN
       write(ionum, '(A, F8.4)') '# Edd_max: ',   param%Edd_max
       write(ionum, '(A, I2)') '# BolCor: ',  param%BolCor
       write(ionum, '(A, L)')  '# delay: ',   param%delay
       write(ionum, '(A, L)')  '# run_all_ssp: ', param%run_all_ssp
       write(ionum, '(A, I3)') '# nwave: ',       param%nwave
       DO i = 1, param%nwave
          IF(param%nwave < 10) THEN
             write(ionum, '(A, I1, A, $)') '# --- ', i, ' = '
          ELSEIF(param%nwave < 100) THEN
             write(ionum, '(A, I2, A, $)') '# --- ', i, ' = '
          ELSE
             write(ionum, '(A, I3, A, $)') '# --- ', i, ' = '
          ENDIF
          write(ionum, '(A)') trim(ssp%bandname(param%iwave(i)))
       ENDDO
    ELSEIF(number == 2) THEN
       write(ionum, '(A, I3.3)') &
            '# index for N-body slice: ', param%izout
       write(ionum, '(A, F8.4)') '# used 1+z: ', param%zsp1
       write(ionum, '(A)') '# sspfile for starburst : '//trim(ssp%filename(1))
       write(ionum, '(A)') '# sspfile for quiescent : '//trim(ssp%filename(2))
       write(ionum, '(2(A, F8.4), A)') '# --- alpha for burst and quies: ', &
            ssp%alp(1), ' (starburst), ', ssp%alp(2), ' (quiescent)'
       write(ionum, '(2(A, F8.4), A)') '# --- R     for burst and quies: ', &
            ssp%R(1), ' (starburst), ', ssp%R(2), ' (quiescent)'
       write(ionum, '(2(A, F8.4), A)') '# --- p     for burst and quies: ', &
            ssp%p(1), ' (starburst), ', ssp%p(2), ' (quiescent)'
       write(ionum, '(A, I3)') '# loopz: ', param%loopz
   ELSEIF(number == 3) THEN
      write(ionum, '(A)') '# ============ MCMC status ============ #'  
      write(ionum, '(A, I8)') '# iMCMC_max: ', mcmc%iMCMC_max
      write(ionum, '(A, L, $)') '# alpst: ', mcmc%alpst
      IF(mcmc%alpst) &
         write(ionum, '(2(A, F8.4),$)') '  min: ', mcmc%alpst_min, &
                                    '    max: ', mcmc%alpst_max
      write(ionum, *)
      write(ionum, '(A, L, $)') '# Vst: ', mcmc%Vst
      IF(mcmc%Vst) &
         write(ionum, '(2(A, F8.4),$)') ' min: ', mcmc%Vst_min, &
                                    '   max: ', mcmc%Vst_max
      write(ionum, *)
      write(ionum, '(A, L, $)') '# tau0st: ', mcmc%tau0st
      IF(mcmc%tau0st) &
         write(ionum, '(2(A, F8.4),$)') ' min: ', mcmc%tau0st_min, &
                                    '   max: ', mcmc%tau0st_max
      write(ionum, *)
      write(ionum, '(A, L, $)') '# eps_SF: ', mcmc%eps_SF
      IF(mcmc%eps_SF) &
         write(ionum, '(2(A, F8.4),$)') ' min: ', mcmc%eps_SF_min, &
                                    '   max: ', mcmc%eps_SF_max
      write(ionum, *)
      write(ionum, '(A, L, $)') '# emax: ', mcmc%emax
      IF(mcmc%emax) &
         write(ionum, '(2(A, F8.4),$)') ' min: ', mcmc%emax_min, &
                                    '   max: ', mcmc%emax_max
      write(ionum, *)
      write(ionum, '(A, L, $)') '# emin: ', mcmc%emin
      IF(mcmc%emin) &
         write(ionum, '(2(A, F8.4),$)') ' min: ', mcmc%emin_min, &
                                    '   max: ', mcmc%emin_max
      write(ionum, *)
      write(ionum, '(A, L, $)') '# Zch: ', mcmc%Zch
      IF(mcmc%Zch) &
         write(ionum, '(2(A, F8.4),$)') ' min: ', mcmc%Zch_min, &
                                    '   max: ', mcmc%Zch_max
      write(ionum, *)
      write(ionum, '(A, L, $)') '# taud_th: ', mcmc%taud_th
      IF(mcmc%taud_th) &
         write(ionum, '(2(A, F8.4),$)') ' min: ', mcmc%taud_th_min, &
                                    '   max: ', mcmc%taud_th_max
      write(ionum, *)
      write(ionum, '(A, L, $)') '# Vhot(1): ', mcmc%Vhot(1)
      IF(mcmc%Vhot(1)) &
         write(ionum, '(2(A, F8.4),$)') ' min: ', mcmc%Vhot_min(1), &
                                    '   max: ', mcmc%Vhot_max(1)
      write(ionum, *)
      write(ionum, '(A, L, $)') '# Vhot(2): ', mcmc%Vhot(2)
      IF(mcmc%Vhot(2)) &
         write(ionum, '(2(A, F8.4),$)') ' min: ', mcmc%Vhot_min(2), &
                                    '   max: ', mcmc%Vhot_max(2)
      write(ionum, *)
      write(ionum, '(A, L, $)') '# alphot(1): ', mcmc%alphot(1)
      IF(mcmc%alphot(1)) &
         write(ionum, '(2(A, F8.4),$)') ' min: ', mcmc%alphot_min(1), &
                                    '   max: ', mcmc%alphot_max(1)
      write(ionum, *)
      write(ionum, '(A, L, $)') '# alphot(2): ', mcmc%alphot(2)
      IF(mcmc%alphot(2)) &
         write(ionum, '(2(A, F8.4),$)') ' min: ', mcmc%alphot_min(2), &
                                    '   max: ', mcmc%alphot_max(2)
      write(ionum, *)
      write(ionum, '(A, L, $)') '# alp_ret: ', mcmc%alp_ret
      IF(mcmc%alp_ret) &
         write(ionum, '(2(A, F8.4),$)') ' min: ', mcmc%alp_ret_min, &
                                    '   max: ', mcmc%alp_ret_max
      write(ionum, *)
      write(ionum, '(A, L, $)') '# fmerge: ', mcmc%fmerge
      IF(mcmc%fmerge) &
         write(ionum, '(2(A, F8.4),$)') ' min: ', mcmc%fmerge_min, &
                                    '   max: ', mcmc%fmerge_max
      write(ionum, *)
      write(ionum, '(A, L, $)') '# fmajor: ', mcmc%fmajor
      IF(mcmc%fmajor) &
         write(ionum, '(2(A, F8.4),$)') ' min: ', mcmc%fmajor_min, &
                                    '   max: ', mcmc%fmajor_max
      write(ionum, *)
      write(ionum, '(A, L, $)') '# Krem: ', mcmc%Krem
      IF(mcmc%Krem) &
         write(ionum, '(2(A, F8.4),$)') ' min: ', mcmc%Krem_min, &
                                    '   max: ', mcmc%Krem_max
      write(ionum, *)
      write(ionum, '(A, L, $)') '# fdiss: ', mcmc%fdiss
      IF(mcmc%fdiss) &
         write(ionum, '(2(A, F8.4),$)') ' min: ', mcmc%fdiss_min, &
                                    '   max: ', mcmc%fdiss_max
      write(ionum, *)
      write(ionum, '(A, L, $)') '# Mh0: ', mcmc%Mh0
      IF(mcmc%Mh0) &
         write(ionum, '(2(A, F8.4),$)') ' min: ', mcmc%Mh0_min, &
                                    '   max: ', mcmc%Mh0_max
      write(ionum, *)
      write(ionum, '(A, L, $)') '# alpha_halo: ', mcmc%alpha_halo
      IF(mcmc%alpha_halo) &
         write(ionum, '(2(A, F8.4),$)') ' min: ', mcmc%alpha_halo_min, &
                                    '   max: ', mcmc%alpha_halo_max
      write(ionum, *)
      write(ionum, '(A, L, $)') '# Vcut: ', mcmc%Vcut
      IF(mcmc%Vcut) &
         write(ionum, '(2(A, F8.4),$)') ' min: ', mcmc%Vcut_min, &
                                    '   max: ', mcmc%Vcut_max
      write(ionum, *)
      write(ionum, '(A, L, $)') '# kappa_croton: ', mcmc%kappa_croton
      IF(mcmc%kappa_croton) &
         write(ionum, '(2(A, F8.4),$)') ' min: ', mcmc%kappa_croton_min, &
                                    '   max: ', mcmc%kappa_croton_max
      write(ionum, *)
      write(ionum, '(A, L, $)') '# eta_croton: ', mcmc%eta_croton
      IF(mcmc%eta_croton) &
         write(ionum, '(2(A, F8.4),$)') ' min: ', mcmc%eta_croton_min, &
                                    '   max: ', mcmc%eta_croton_max
      write(ionum, *)
      write(ionum, '(A, L, $)') '# alp_bower: ', mcmc%alp_bower
      IF(mcmc%alp_bower) &
         write(ionum, '(2(A, F8.4),$)') ' min: ', mcmc%alp_bower_min, &
                                    '   max: ', mcmc%alp_bower_max
      write(ionum, *)
      write(ionum, '(A, L, $)') '# eps_bower: ', mcmc%eps_bower
      IF(mcmc%eps_bower) &
         write(ionum, '(2(A, F8.4),$)') ' min: ', mcmc%eps_bower_min, &
                                    '   max: ', mcmc%eps_bower_max
      write(ionum, *)
      write(ionum, '(A, L, $)') '# em: ', mcmc%em
      IF(mcmc%em) &
         write(ionum, '(2(A, F8.4),$)') ' min: ', mcmc%em_min, &
                                    '   max: ', mcmc%em_max
      write(ionum, *)
      write(ionum, '(A, L, $)') '# fdi: ', mcmc%fdi
      IF(mcmc%fdi) &
         write(ionum, '(2(A, F8.4),$)') ' min: ', mcmc%fdi_min, &
                                    '   max: ', mcmc%fdi_max
      write(ionum, *)
      write(ionum, '(A, L, $)') '# fbh(1): ', mcmc%fbh(1)
      IF(mcmc%fbh(1)) &
         write(ionum, '(2(A, F8.4),$)') ' min: ', mcmc%fbh_min(1), &
                                    '   max: ', mcmc%fbh_max(1)
      write(ionum, *)
      write(ionum, '(A, L, $)') '# fbh(2): ', mcmc%fbh(2)
      IF(mcmc%fbh(2)) &
         write(ionum, '(2(A, F8.4),$)') ' min: ', mcmc%fbh_min(2), &
                                    '   max: ', mcmc%fbh_max(2)
      write(ionum, *)
      write(ionum, '(A, L, $)') '# Mbhseed: ', mcmc%Mbhseed
      IF(mcmc%Mbhseed) &
         write(ionum, '(2(A, F8.4),$)') ' min: ', mcmc%Mbhseed_min, &
                                    '   max: ', mcmc%Mbhseed_max
      write(ionum, *)
      write(ionum, '(A, L, $)') '# alpha_ad: ', mcmc%alpha_ad
      IF(mcmc%alpha_ad) &
         write(ionum, '(2(A, F8.4),$)') ' min: ', mcmc%alpha_ad_min, &
                                    '   max: ', mcmc%alpha_ad_max
      write(ionum, *)
      write(ionum, '(A, L, $)') '# gamma_ad: ', mcmc%gamma_ad
      IF(mcmc%gamma_ad) &
         write(ionum, '(2(A, F8.4),$)') ' min: ', mcmc%gamma_ad_min, &
                                    '   max: ', mcmc%gamma_ad_max
      write(ionum, *)
      write(ionum, '(A, L, $)') '# Edd_max: ', mcmc%Edd_max
      IF(mcmc%Edd_max) &
         write(ionum, '(2(A, F8.4),$)') ' min: ', mcmc%Edd_max_min, &
                                    '   max: ', mcmc%Edd_max_max
      write(ionum, *)
      write(ionum, '(A)') '# ===================================== #'  
   ENDIF
  END SUBROUTINE WriteInputParameters
!!$============================================================================
  SUBROUTINE CreateOutputFileNames(run_all_ssp, run_redshift)
    LOGICAL, INTENT(IN) :: run_all_ssp
    INTEGER, INTENT(IN) :: run_redshift ! 1:run at a redshift
                                        ! 2:run at all redshift
                                        ! 3:MCMC
    CHARACTER(LEN=2) :: ci

    IF(run_all_ssp) THEN

    ELSE
       IF(inode == 0) print '(A)', '# Output File Names:'
       IF(run_redshift == 1) THEN ! run at a redshift
          DO i = 0, NFILE
             write(ci, '(I1.1)') i
             file_out(i) = trim(fbase)//trim(file_o)//'_'//ci
             file_out(i) = trim(trim(file_out(i))//'.dat')
             IF(inode == 0) print '(A)', '  --- '//trim(file_out(i))
          ENDDO

          ! for parallel run, added by Makiya (2015/09/27)
          write(ci, '(I2.2)') inode+1
          ! galaxy catalog
          file_catalog = trim(trim(fbase)//trim(file_o)//'_catalog_'//ci//'.dat')
          ! quasar catalog
          file_catalog_q = trim(trim(fbase)//trim(file_o)//'_catalog_q_'//ci//'.dat')
          ! target catalog
          IF(param%traceIDs == 1) &
             file_catalog_t = trim(trim(fbase)//trim(file_o)//'_catalog_t_'//ci//'.dat')
       ELSEIF(run_redshift == 2) THEN ! 2:run at all redshift
         ! i = 0 (input parameters)
           write(ci, '(I1.1)') 0
           file_out(i) = trim(fbase)//trim(file_o)//'_'//ci
           file_out(i) = trim(trim(file_out(i))//'.dat')
           print '(A)', ' --- '//trim(file_out(i))

          DO i = 1, NFILE
             write(ci, '(I1.1)') i
             file_out(i) = trim(fbase)//trim(file_o)//'_'//ci
             file_out(i) = trim(trim(file_out(i))//'_***.dat')
             print '(A)', ' --- '//trim(file_out(i))
          ENDDO

          ! for parallel run, added by Makiya (2015/09/27)
          write(ci, '(I2.2)') inode+1
          ! galaxy catalog
          file_catalog = trim(fbase)//trim(file_o)//'_catalog_'//ci
          file_catalog = trim(trim(file_catalog)//'_***.dat')
          ! quasar catalog
          file_catalog_q = trim(fbase)//trim(file_o)//'_catalog_q_'//ci
          file_catalog_q = trim(trim(file_catalog_q)//'_***.dat')
       ELSEIF(run_redshift == 3 .and. inode == 0) THEN ! 3:MCMC
          write(ci, '(I1.1)') 0
          file_out(0) = trim(fbase)//trim(file_o)//'_'//ci
          file_out(0) = trim(trim(file_out(0))//'.dat')
          print '(A)', '  --- '//trim(file_out(0))
       ENDIF
    ENDIF
  END SUBROUTINE CreateOutputFileNames
!!$============================================================================
  SUBROUTINE AllocateDistributionFunctions(nz_end)
    INTEGER, INTENT(IN) :: nz_end
    INTEGER :: iz, j, N1, N2, Nbin, ier
    CHARACTER(LEN=50) :: cerr &
         = '# AllocateDistributionFunctions: fail to allocate '

    ! DistributionFunction1
    allocate(lf(nz_end),   stat=ier); call CheckIerr(ier, trim(cerr)//' lf')
    allocate(lf_d(nz_end), stat=ier); call CheckIerr(ier, trim(cerr)//' lf_d')
    allocate(lf_q(nz_end), stat=ier); call CheckIerr(ier, trim(cerr)//' lf_q')
    allocate(mf(nz_end),   stat=ier); call CheckIerr(ier, trim(cerr)//' mf')
    ! DistributionFunction2
    allocate(ml(nz_end), stat=ier); call CheckIerr(ier, trim(cerr)//' ml')
    allocate(MbhMbulge, stat=ier); call CheckIerr(ier, trim(cerr)//' MbhMbulge')
    allocate(DiskScale, stat=ier); call CheckIerr(ier, trim(cerr)//' DiskScale')
    allocate(SphScale, stat=ier); call CheckIerr(ier, trim(cerr)//' SphScale')

    DO iz = 1, nz_end
       ! --- luminosity functions
       N1 = param%nwave; N2 = 4; Nbin = NbinLF
       lf(iz)%Nbin   = Nbin; lf(iz)%N1   = N1; lf(iz)%N2   = N2
       lf_d(iz)%Nbin = Nbin; lf_d(iz)%N1 = N1; lf_d(iz)%N2 = N2
       allocate(lf(iz)%bin(Nbin),   stat=ier); call CheckIerr(ier, &
            trim(cerr)//' lf(iz)%bin')
       allocate(lf_d(iz)%bin(Nbin), stat=ier); call CheckIerr(ier, &
            trim(cerr)//' lf_d(iz)%bin')
       allocate(lf(iz)%n(N1,N2,Nbin),   stat=ier); call CheckIerr(ier, &
            trim(cerr)//' lf(iz)%n')
       allocate(lf_d(iz)%n(N1,N2,Nbin), stat=ier); call CheckIerr(ier, &
            trim(cerr)//' lf_d(iz)%n')
       allocate(lf(iz)%n_brst(N1,N2,Nbin),   stat=ier); call CheckIerr(ier, &
            trim(cerr)//' lf(iz)%n_brst')
       allocate(lf_d(iz)%n_brst(N1,N2,Nbin), stat=ier); call CheckIerr(ier, &
            trim(cerr)//' lf_d(iz)%n_brst')
       ! lf(:)%n(i,j,k) and lf_d(:)%n(i,j,k)
       ! --- i:band, j:mor(1:E,2:S0,3:S,4:all), k:mag
       lf(iz)%step = stepLF; lf_d(iz)%step = stepLF
       lf(iz)%base = baseLF; lf_d(iz)%base = baseLF
       lf(iz)%invstep = 1.d0 / lf(iz)%step; lf_d(iz)%invstep = 1.d0 / lf_d(iz)%step
       DO j = 1, Nbin
          lf(iz)%bin(j)   = lf(iz)%step   * dble(j) + lf(iz)%base
          lf_d(iz)%bin(j) = lf_d(iz)%step * dble(j) + lf_d(iz)%base
       ENDDO
       lf(iz)%n(:,:,:)   = 0.d0; lf(iz)%n_brst(:,:,:)   = 0.d0
       lf_d(iz)%n(:,:,:) = 0.d0; lf_d(iz)%n_brst(:,:,:) = 0.d0

       ! --- QSO luminosity functions
       N1 = param%nwaveq; N2 = 1; Nbin = NbinLF
       lf_q(iz)%Nbin = Nbin; lf_q(iz)%N1 = N1; lf_q(iz)%N2 = N2
       allocate(lf_q(iz)%bin(Nbin),     stat=ier); call CheckIerr(ier, &
            trim(cerr)//' lf_q(iz)%bin')
       allocate(lf_q(iz)%n(N1,N2,Nbin), stat=ier); call CheckIerr(ier, &
            trim(cerr)//' lf_q(iz)%n')
       ! lf_q(:)%n(i,j,k)
       ! --- i:band, j:mor(1:E,2:S0,3:S,4:all), k:mag
       lf_q(iz)%step = stepLF; lf_q(iz)%base = baseLF
       lf_q(iz)%invstep = 1.d0 / lf_q(iz)%step
       DO j = 1, Nbin
          lf_q(iz)%bin(j) = lf_q(iz)%step * dble(j) + lf_q(iz)%base
       ENDDO
       lf_q(iz)%n(:,:,:) = 0.d0

       ! --- mass functions
       N1 = NtypeMF; N2 = 4; Nbin = NbinMF
       mf(iz)%Nbin = Nbin; mf(iz)%N1 = N1; mf(iz)%N2 = N2
       allocate(mf(iz)%bin(Nbin),       stat=ier); call CheckIerr(ier, &
            trim(cerr)//' mf(iz)%bin')
       allocate(mf(iz)%n(N1, N2, Nbin), stat=ier); call CheckIerr(ier, &
            trim(cerr)//' mf(iz)%n')
       allocate(mf(iz)%n_brst(N1, N2, Nbin), stat=ier); call CheckIerr(ier, &
            trim(cerr)//' mf(iz)%n_brst')
       ! mf(:)%n(i,j,k)
       ! --- i:MFtype(1:Mstar^bulge, 2:Mstar^disk, 3:Mcold, 4:MBH,
       !              5:Mhot, 6:Mstar^bulge+MBH, 7:Mstar^disk+Mcold,
       !              8:Mstar^bulge+Mstar^disk+MBH+Mcold,
       !              9:Mstar^bulge+Mstar^disk, 10:hydrogen gas, 11:Mhalo),
       !     j:mor(1:E,2:S0,3:S,4:all), k:mass
       mf(iz)%step = stepMF; mf(iz)%base = param%log10munit + 2.d0 * log10(param%h)
       mf(iz)%invstep = 1.d0 / mf(iz)%step
       DO j = 1, Nbin
          mf(iz)%bin(j) = mf(iz)%step * dble(j)
       ENDDO
       mf(iz)%n(:,:,:) = 0.d0; mf(iz)%n_brst(:,:,:) = 0.d0

       ! calculation of ML is executed only for the case of param%ML = .true.
       IF(param%ML) THEN
          ! --- M_HI / L_B as a function of L_B
          N1 = 2; N2 = 4; Nbin = NbinML
          ml(iz)%Nbin = Nbin; ml(iz)%N1 = N1; ml(iz)%N2 = N2
          allocate(ml(iz)%bin(Nbin),       stat=ier); call CheckIerr(ier, &
               trim(cerr)//' ml(iz)%bin')
          allocate(ml(iz)%n(N1,N2,Nbin),   stat=ier); call CheckIerr(ier, &
               trim(cerr)//' ml(iz)%n')
          allocate(ml(iz)%xn(N1,N2,Nbin),  stat=ier); call CheckIerr(ier, &
               trim(cerr)//' ml(iz)%xn')
          allocate(ml(iz)%xxn(N1,N2,Nbin), stat=ier); call CheckIerr(ier, &
               trim(cerr)//' ml(iz)%xxn')
          IF(param%n_forest > 1000) THEN
             allocate(ml(iz)%x(NgalMax2,N2,Nbin), stat=ier); call CheckIerr(ier, &
                  trim(cerr)//' ml(iz)%x')
          ELSE
             allocate(ml(iz)%x(NgalMax,N2,Nbin), stat=ier); call CheckIerr(ier, &
                  trim(cerr)//' ml(iz)%x')
          ENDIF

          ! ml(:)%n(i,j,k), ml(:)%xn(i,j,k), ml(:)%xxn(i,j,k)
          ! --- i:type(1:mean+sigma, 2:median+quartiles),
          !     j:mor(1:E,2:S0,3:S,4:all), k:M_B(bulge+disk)
          ml(iz)%step = stepML; ml(iz)%base = baseML
          ml(iz)%invstep = 1.d0 / ml(iz)%step
          DO j = 1, Nbin
             ml(iz)%bin(j) = ml(iz)%step * dble(j) + ml(iz)%base
          ENDDO
          ml(iz)%n(:,:,:) = 0.d0; ml(iz)%xn(:,:,:) = 0.d0
          ml(iz)%xxn(:,:,:) = 0.d0; ml(iz)%x(:,:,:) = 0.d0
       ENDIF
    ENDDO

    ! calculation of Mbh as a function of Mbulge is executed only for
    !   the case of param%Mbh = .true.
    IF(param%Mbh .or. param%run_type == 3) THEN
       N1 = 1; N2 = 4; Nbin = NbinML
       MbhMbulge%Nbin = Nbin; MbhMbulge%N1 = N1; MbhMbulge%N2 = N2
       allocate(MbhMbulge%bin(Nbin),          stat=ier); call &
            CheckIerr(ier, trim(cerr)//'MbhMbulge%bin')
       allocate(MbhMbulge%n(N1,N2,Nbin),      stat=ier); call &
            CheckIerr(ier, trim(cerr)//'MbhMbulge%n')
       allocate(MbhMbulge%xn(N1,N2,Nbin),     stat=ier); call &
            CheckIerr(ier, trim(cerr)//'MbhMbulge%xn')
       allocate(MbhMbulge%xxn(N1,N2,Nbin),    stat=ier); call &
            CheckIerr(ier, trim(cerr)//'MbhMbulge%xxn')
       ! MbhMbulge(:)%n(i,j,k), MbhMbulge(:)%xn(i,j,k), MbhMbulge(:)%xxn(i,j,k)
       ! --- i:type(1:mean+sigma),
       !     j:mor(1:E,2:S0,3:S,4:all), k:M_B(bulge+disk)
       MbhMbulge%step = stepMBhMbulge; MbhMbulge%base = baseMbhMbulge
       MbhMbulge%invstep = 1.d0 / MbhMbulge%step
       DO j = 1, Nbin
          MbhMbulge%bin(j) = MbhMbulge%step * dble(j) + MbhMbulge%base
       ENDDO
       MbhMbulge%n(:,:,:) = 0.d0; MbhMbulge%xn(:,:,:) = 0.d0
       MbhMbulge%xxn(:,:,:) = 0.d0
    ENDIF

    ! For MCMC. Scaling relations of disk & bulge
    IF(param%run_type == 3) THEN
       N1 = 1; N2 = 4; Nbin = NbinML
       DiskScale%Nbin = Nbin; DiskScale%N1 = N1; DiskScale%N2 = N2
       allocate(DiskScale%bin(Nbin),          stat=ier); call &
            CheckIerr(ier, trim(cerr)//'DiskScale%bin')
       allocate(DiskScale%n(N1,N2,Nbin),      stat=ier); call &
            CheckIerr(ier, trim(cerr)//'DiskScale%n')
       allocate(DiskScale%xn(N1,N2,Nbin),     stat=ier); call &
            CheckIerr(ier, trim(cerr)//'DiskScale%xn')
       allocate(DiskScale%xxn(N1,N2,Nbin),    stat=ier); call &
            CheckIerr(ier, trim(cerr)//'DiskScale%xxn')
       ! DiskScale(:)%n(i,j,k), DiskScale(:)%xn(i,j,k), DiskScale(:)%xxn(i,j,k)
       ! --- i:type(1:mean+sigma),
       !     j:mor(1:E,2:S0,3:S,4:all), k:M_B(bulge+disk)
       DiskScale%step = stepDiskScale; DiskScale%base = baseDiskScale
       DiskScale%invstep = 1.d0 / DiskScale%step
       DO j = 1, Nbin
          DiskScale%bin(j) = DiskScale%step * dble(j) + DiskScale%base
       ENDDO
       DiskScale%n(:,:,:) = 0.d0; DiskScale%xn(:,:,:) = 0.d0
       DiskScale%xxn(:,:,:) = 0.d0

       N1 = 1; N2 = 4; Nbin = NbinML
       SphScale%Nbin = Nbin; SphScale%N1 = N1; SphScale%N2 = N2
       allocate(SphScale%bin(Nbin),          stat=ier); call &
            CheckIerr(ier, trim(cerr)//'SphScale%bin')
       allocate(SphScale%n(N1,N2,Nbin),      stat=ier); call &
            CheckIerr(ier, trim(cerr)//'SphScale%n')
       allocate(SphScale%xn(N1,N2,Nbin),     stat=ier); call &
            CheckIerr(ier, trim(cerr)//'SphScale%xn')
       allocate(SphScale%xxn(N1,N2,Nbin),    stat=ier); call &
            CheckIerr(ier, trim(cerr)//'SphScale%xxn')
       ! SphScale(:)%n(i,j,k), SphScale(:)%xn(i,j,k), SphScale(:)%xxn(i,j,k)
       ! --- i:type(1:mean+sigma),
       !     j:mor(1:E,2:S0,3:E+S0,4:all), k:M_B(bulge+disk)
       SphScale%step = stepSphScale; SphScale%base = baseSphScale
       SphScale%invstep = 1.d0 / SphScale%step
       DO j = 1, Nbin
          SphScale%bin(j) = SphScale%step * dble(j) + SphScale%base
       ENDDO
       SphScale%n(:,:,:) = 0.d0; SphScale%xn(:,:,:) = 0.d0
       SphScale%xxn(:,:,:) = 0.d0
    ENDIF
  END SUBROUTINE AllocateDistributionFunctions
!!$============================================================================
  SUBROUTINE DeallocateDistributionFunctions(nz)
    INTEGER, INTENT(IN) :: nz
    ! for luminosity functions
    INTEGER :: iz, ier
    CHARACTER(LEN=60) :: cerr &
         = '# DeallocateDistributionFunctions: fail to deallocate '


    DO iz = 1, nz
       ! --- luminosity functions
       deallocate(lf(iz)%bin,   stat=ier); call CheckIerr(ier, &
            trim(cerr)//' lf(iz)%bin')
       deallocate(lf_d(iz)%bin, stat=ier); call CheckIerr(ier, &
            trim(cerr)//' lf_d(iz)%bin')
       deallocate(lf(iz)%n,     stat=ier); call CheckIerr(ier, &
            trim(cerr)//' lf(iz)%n')
       deallocate(lf_d(iz)%n,   stat=ier); call CheckIerr(ier, &
            trim(cerr)//' lf_d(iz)%n')
       deallocate(lf(iz)%n_brst,   stat=ier); call CheckIerr(ier, &
            trim(cerr)//' lf(iz)%n_brst')
       deallocate(lf_d(iz)%n_brst, stat=ier); call CheckIerr(ier, &
            trim(cerr)//' lf_d(iz)%n_brst')
       deallocate(lf_q(iz)%bin, stat=ier); call CheckIerr(ier, &
            trim(cerr)//' lf_q(iz)%bin') ! add by Makiya
       deallocate(lf_q(iz)%n,   stat=ier); call CheckIerr(ier, &
            trim(cerr)//' lf_d(iz)%n') ! add by Makiya

       ! --- mass functions
       deallocate(mf(iz)%bin,     stat=ier); call CheckIerr(ier, &
            trim(cerr)//' mf(iz)%bin')
       deallocate(mf(iz)%n,       stat=ier); call CheckIerr(ier, &
            trim(cerr)//' mf(iz)%n')
       deallocate(mf(iz)%n_brst, stat=ier); call CheckIerr(ier, &
            trim(cerr)//' mf(iz)%n_brst')

       ! --- M_HI / L_B as a function of L_B
       IF(param%ML) THEN
          deallocate(ml(iz)%bin, stat=ier); call CheckIerr(ier, &
               trim(cerr)//' ml(iz)%bin')
          deallocate(ml(iz)%n,   stat=ier); call CheckIerr(ier, &
               trim(cerr)//' ml(iz)%n')
          deallocate(ml(iz)%xn,  stat=ier); call CheckIerr(ier, &
               trim(cerr)//' ml(iz)%xn')
          deallocate(ml(iz)%xxn, stat=ier); call CheckIerr(ier, &
               trim(cerr)//' ml(iz)%xxn')
          deallocate(ml(iz)%x,   stat=ier); call CheckIerr(ier, &
               trim(cerr)//' ml(iz)%x')
       ENDIF
    ENDDO

    ! --- Mbh as a function of Mbulge added by Makiya 20140319
    IF(param%Mbh .or. param%run_type == 3) THEN
       deallocate(MbhMbulge%bin, stat=ier); call CheckIerr(ier, &
            trim(cerr)//' MbhMbulge%bin')
       deallocate(MbhMbulge%n,   stat=ier); call CheckIerr(ier, &
            trim(cerr)//' MbhMbulge%n')
       deallocate(MbhMbulge%xn,  stat=ier); call CheckIerr(ier, &
            trim(cerr)//' MbhMbulge%xn')
       deallocate(MbhMbulge%xxn, stat=ier); call CheckIerr(ier, &
            trim(cerr)//' MbhMbulge%xxn')
    ENDIF

    IF(param%run_type == 3) THEN
       deallocate(DiskScale%bin, stat=ier); call CheckIerr(ier, &
            trim(cerr)//' DiskScale%bin')
       deallocate(DiskScale%n,   stat=ier); call CheckIerr(ier, &
            trim(cerr)//' DiskScale%n')
       deallocate(DiskScale%xn,  stat=ier); call CheckIerr(ier, &
            trim(cerr)//' DiskScale%xn')
       deallocate(DiskScale%xxn, stat=ier); call CheckIerr(ier, &
            trim(cerr)//' DiskScale%xxn')

       deallocate(SphScale%bin, stat=ier); call CheckIerr(ier, &
            trim(cerr)//' SphScale%bin')
       deallocate(SphScale%n,   stat=ier); call CheckIerr(ier, &
            trim(cerr)//' SphScale%n')
       deallocate(SphScale%xn,  stat=ier); call CheckIerr(ier, &
            trim(cerr)//' SphScale%xn')
       deallocate(SphScale%xxn, stat=ier); call CheckIerr(ier, &
            trim(cerr)//' SphScale%xxn')
    ENDIF

    ! DistributionFunction1
    deallocate(lf,   stat=ier); call CheckIerr(ier, trim(cerr)// 'lf')
    deallocate(lf_d, stat=ier); call CheckIerr(ier, trim(cerr)//' lf_d')
    deallocate(lf_q, stat=ier); call CheckIerr(ier, trim(cerr)//' lf_q')
                                ! added by Makiya
    deallocate(mf,   stat=ier); call CheckIerr(ier, trim(cerr)//' mf')

    ! DistributionFunction2
    deallocate(ml, stat=ier); call CheckIerr(ier, trim(cerr)//' ml')
    deallocate(MbhMbulge, stat=ier); call CheckIerr(ier, trim(cerr)//' MbhMbulge')
    deallocate(DiskScale, stat=ier); call CheckIerr(ier, trim(cerr)//' DiskScale')
    deallocate(SphScale, stat=ier); call CheckIerr(ier, trim(cerr)//' SphScale')
  END SUBROUTINE DeallocateDistributionFunctions
!!$============================================================================
  SUBROUTINE ReadSSPFile(loopz, nz_end, nz)
    INTEGER, INTENT(IN) :: loopz, nz_end, nz
    INTEGER :: i, i_num_filter, j, k, l, ier, i_file
    CHARACTER(LEN=3)   :: ci
    CHARACTER(LEN=4)   :: c_num_filter
    CHARACTER(LEN=100) :: charact
    CHARACTER(LEN=500) :: filename
    CHARACTER(LEN=50)  :: cerr_o = '# ReadSSPFile: fail to open file '
    CHARACTER(LEN=50)  :: cerr_a = '# ReadSSPFile: fail to allocate '


    NSFR = NSFR_BASE(param%type_ssp); NCHEM = NCHEM_BASE(param%type_ssp)
    allocate(ssp%time(NSFR), stat=ier); call CheckIerr(ier, &
         trim(cerr_a)//'ssp%time')
    allocate(ssp%chem(NCHEM), stat=ier); call CheckIerr(ier, &
         trim(cerr_a)//'ssp%chem')
    allocate(ssp%lumi(2, TNWAVE_ALL, NSFR, NCHEM), stat=ier); call &
         CheckIerr(ier, trim(cerr_a)//'ssp%lumi')

    write(ci, '(I3.3)') mrgp%num_step - loopz
    DO k = 1, 2
       ! --- SSPs for observer-frame ---
       i = index(ssp%filename(k), '.dat') - 3
       ssp%filename(k)(i:i+2) = ci
       filename = trim(fbase_nbody)//trim(param%file_nbody)//'/sspfiles/'//&
            trim(ssp%filename(k))
       i_file = 1
       open(i_file, file=trim(filename), status='old', iostat=ier)
       call CheckIerr(ier, trim(cerr_o)//trim(filename))
       IF(inode == 0 .and. nz == 1 .and. param%run_type /= 3) THEN
          print '(A, I1, A)', &
                '# Reading a sspfile(', k, ') for observer-frame: '//&
                trim(ssp%filename(k))
       ENDIF
       read(i_file, *) ! Adopted population synthesis code
       read(i_file, *) ! How to determine the calculated redshifts
       read(i_file, *) ! redshift, mag. system, time and metallicity slices
       read(i_file, '(A)') header ! IMF

       IF(inode == 0 .and. param%run_type /= 3) THEN
          DO i = 1, NFILE
             ! --- write the adopted population synthesis model and IMF
             IF(k == 1) THEN
                IF(param%run_type /= 2) write(ionum+i, '(A)') &
                     '# '//base_ssp(param%type_ssp)(2:&
                     index(base_ssp(param%type_ssp), ' ')-1)//' SSP is used'
             ENDIF
             IF(param%run_type /= 2) write(ionum+i, '(A)') &
                  c_sftype(k)(1:index(c_sftype(k),'  '))//header(index(header,&
                  '# ')+1:index(header, '    ')-1)
          ENDDO

          IF(k == 1) THEN
             ! galaxy catalog
             IF(param%run_type == 1) THEN
                write(ionum+NFILE+1, '(A)') &
                  '# '//base_ssp(param%type_ssp)(2:&
                  index(base_ssp(param%type_ssp), ' ')-1)//' SSP is used'
             ELSE IF(param%run_type == 2) THEN
                write(ionum+NFILE*nz_end+nnode*(nz-1)+1, '(A)') &
                  '# '//base_ssp(param%type_ssp)(2:&
                  index(base_ssp(param%type_ssp), ' ')-1)//' SSP is used'
             ENDIF
             ! quasar catalog
             IF(param%run_type == 1) THEN
                write(ionum+NFILE+nnode+1, '(A)') &
                  '# '//base_ssp(param%type_ssp)(2:&
                  index(base_ssp(param%type_ssp), ' ')-1)//' SSP is used'
             ELSE IF(param%run_type == 2) THEN
                write(ionum+NFILE*nz_end+nnode*(nz_end+nz-1)+1, '(A)') &
                  '# '//base_ssp(param%type_ssp)(2:&
                  index(base_ssp(param%type_ssp), ' ')-1)//' SSP is used'
             ENDIF
          ENDIF
          ! galaxy catalog
          IF(param%run_type == 1) THEN
            write(ionum+NFILE+1, '(A)') &
               c_sftype(k)(1:index(c_sftype(k),'  '))//header(index(header,&
               '# ')+1:index(header, '    ')-1)
          ELSE IF(param%run_type == 2) THEN
            write(ionum+NFILE*nz_end+nnode*(nz-1)+1, '(A)') &
               c_sftype(k)(1:index(c_sftype(k),'  '))//header(index(header,&
               '# ')+1:index(header, '    ')-1)
          ENDIF
          ! quasar catalog
          IF(param%run_type == 1) THEN
             write(ionum+NFILE+nnode+1, '(A)') &
               c_sftype(k)(1:index(c_sftype(k),'  '))//header(index(header,&
               '# ')+1:index(header, '    ')-1)
          ELSE IF(param%run_type == 2) THEN
             write(ionum+NFILE*nz_end+nnode*(nz_end+nz-1)+1, '(A)') &
               c_sftype(k)(1:index(c_sftype(k),'  '))//header(index(header,&
               '# ')+1:index(header, '    ')-1)
          ENDIF
       ENDIF

       read(i_file, '(A)') charact ! the number of filters
       c_num_filter = charact(index(charact, ':')+1:len_trim(charact))
       read(c_num_filter, *) i_num_filter

       read(i_file, *) ! captions for each column
       read(i_file, *) ! captions for the unit of column (3)-(47)
       read(i_file, *) ! central wavelength is written in the following
       DO i = 1, i_num_filter
          ! read the central wavelength of each filter [um]
          read(i_file, '(A)') charact
          j = index(charact, ': ') + 1
          charact = charact(j:len_trim(charact))
          read(charact, *) ssp%lam_c(i)
       ENDDO

       DO j = 1, NCHEM
          read(i_file, *) ! '# index *: Z = *****'
          IF(j == 5) THEN ! solar metallicity
             read(i_file, *) ssp%alp(k), ssp%R(k), ssp%p(k)
             ssp%y(k) = ssp%p(k) / ssp%alp(k)
          ELSE
             read(i_file, *) x, x, x
          ENDIF

          DO i = 1, NSFR
!!$     t Z NWAVE_ALL*lumi
             read(i_file, *) ssp%time(i), & ! (1)time [Myr]
                  ssp%chem(j), &            ! (2)metallicity
                  (ssp%lumi(k,l,i,j), l=1,NWAVE_ALL)
!!$                    print *, (ssp%lumi(k,l,i,j), l=1,NWAVE_ALL)
          ENDDO
          read(i_file, *); read(i_file, *)
       ENDDO
       close(i_file)

       ! --- SSPs for rest-frame ---
       i = index(ssp%filename(k), '.dat') - 3
       ssp%filename(k)(i:i+2) = 'z00'
       filename = trim(fbase_nbody)//trim(param%file_nbody)//'/sspfiles/'//&
            trim(ssp%filename(k))
       i_file = 1
       open(i_file, file = trim(filename), status='old', iostat=ier)
       call CheckIerr(ier, trim(cerr_o)//trim(filename))
       IF(inode == 0 .and. nz == 1 .and. param%run_type /= 3) THEN
          print '(A, I1, A)', &
               '# Reading a sspfile(', k, ') for rest-frame: '//&
               trim(ssp%filename(k))
       ENDIF
       read(i_file, *) ! Adopted population synthesis code
       read(i_file, *) ! How to determine the calculated redshifts
       read(i_file, *) ! redshift, mag. system, time and metallicity slices
       read(i_file, '(A)') header ! IMF
       read(i_file, '(A)') charact ! the number of filters
       c_num_filter = charact(index(charact, ':')+1:len_trim(charact))
       read(c_num_filter, *) i_num_filter
       read(i_file, *) ! captions for each column
       read(i_file, *) ! captions for the unit of column (3)-(47)
       read(i_file, *) ! central wavelength is written in the following
       DO i = i_num_filter+1, 2*i_num_filter
          ! read the central wavelength of each filter [um]
          read(i_file, '(A)') charact
          j = index(charact, ': ') + 1
          charact = charact(j:len_trim(charact))
          read(charact, *) ssp%lam_c(i)
       ENDDO

       DO j = 1, NCHEM
          read(i_file, *) ! '# index *: Z = *****'
          IF(param%type_ssp == 2) THEN ! PEGASE
             read(i_file, *)
          ENDIF

          DO i = 1, NSFR
!!$     t Z NWAVE_ALL*lumi
             read(i_file, *) x, &! (1)time [Myr]
                  x, &      ! (2)metallicity
                  (ssp%lumi(k,l,i,j), l=NWAVE_ALL+1, TNWAVE_ALL)
!!$                    print *, (ssp%lumi(k,l,i,j), l=1,NWAVE_ALL)
          ENDDO
          read(i_file, *); read(i_file, *)
       ENDDO
       close(i_file)

       IF(loopz == 1 .and. param%run_type /= 3) THEN
          DO i = 1, param%nwave
             IF(inode == 0) print '(A, I3, A, F6.3, A)', '  --- ', i, ': '//&
                  trim(ssp%bandname(param%iwave(i)))//'-band (', &
                  ssp%lam_c(param%iwave(i)), '[um])'
          ENDDO
       ENDIF
    ENDDO
    DO k = 1, 2
       i = index(ssp%filename(k), '.dat') - 3
       ssp%filename(k)(i:i+2) = ci
    ENDDO
    ssp%time(:) = dble(ssp%time(:)) * 1.d-3 / param%th
                  ! [Myr] --> [hubble time]
  END SUBROUTINE ReadSSPFile
!!$============================================================================
  SUBROUTINE ReadMRGB(iforest,nz)
    INTEGER, INTENT(IN) :: iforest, nz
    INTEGER :: i, j, ier, i_file
    CHARACTER(LEN=500) :: file_mrgb
    CHARACTER(LEN=5)   :: ci
    CHARACTER(LEN=50)  :: cerr = '# ReadMRGB: fail to allocate '


    IF(param%n_forest >= 10000) THEN
       file_mrgb = trim(fbase_nbody)//trim(param%file_nbody)//'/'&
            //trim(param%file_nbody)//'.00000.mrgb'
       write(ci, '(I5.5)') iforest ! This should be I5.5 (.5 is important)

       i = index(file_mrgb, '.mrgb') - 5
       file_mrgb(i:i+4) = ci
       IF(nz == 1 .and. param%run_type /= 3) &
          print '(A, I3, 2(A, I5), A)', '# inode = ', inode, ' : (', iforest+1, &
               '/', param%n_forest, ') Reading a N-body Data File: '//trim(file_mrgb)
    ELSEIF(param%n_forest >= 1000) THEN
       file_mrgb = trim(fbase_nbody)//trim(param%file_nbody)//'/'&
            //trim(param%file_nbody)//'.0000.mrgb'
       write(ci, '(I4.4)') iforest ! This should be I4.4 (.4 is important)

       i = index(file_mrgb, '.mrgb') - 4
       file_mrgb(i:i+3) = ci
       IF(nz == 1 .and. param%run_type /= 3) &
          print '(A, I3, 2(A, I4), A)', '# inode = ', inode, ' : (', iforest+1, &
               '/', param%n_forest, ') Reading a N-body Data File: '//trim(file_mrgb)
    ELSE
       file_mrgb = trim(fbase_nbody)//trim(param%file_nbody)//'/'&
            //trim(param%file_nbody)//'.000.mrgb'
       write(ci, '(I3.3)') iforest ! This should be I3.3 (.3 is important)

       i = index(file_mrgb, '.mrgb') - 3
       file_mrgb(i:i+2) = ci
       IF(nz == 1 .and. param%run_type /= 3) &
          print '(A, I3, 2(A, I3), A)', '# inode = ', inode, ' : (', iforest+1, &
               '/', param%n_forest, ') Reading a N-body Data File: '//trim(file_mrgb)
    ENDIF

    i_file = 1
    open(i_file, file=trim(file_mrgb), status='old', form='unformatted', iostat=ier)
    call CheckIerr(ier, 'ReadMRGB: fail to open file: '//trim(file_mrgb))
    read(i_file) mrgp%num_step, mrgt%num_tot

    ! --- allocate and read the quantities in the structure named "mrgp"
    j = mrgp%num_step
    allocate(mrgp%time(j), stat=ier); call CheckIerr(ier, &
         trim(cerr)//' mrtp%num_time')
    read(i_file) (mrgp%time(i),    i = 1, j)
    mrgp%time(:) = mrgp%time(:) * 1.d-3 / param%th ! [Myr] --> [hubble time]
    allocate(mrgp%zp1ar(j), stat=ier); call CheckIerr(ier, &
         trim(cerr)//' mrgp%zp1ar')
    read(i_file) (mrgp%zp1ar(i),   i = 1, j) ! mrgp%zp1ar(i) decreases
                                             !  with increasing i
    allocate(mrgp%num_now(j), stat=ier); call CheckIerr(ier, &
         trim(cerr)//' mrgp%num_now')
    read(i_file) (mrgp%num_now(i), i = 1, j)
    allocate(mrgp%st_halo(j), stat=ier); call CheckIerr(ier, &
         trim(cerr)//' mrgp%st_halo')
    read(i_file) (mrgp%st_halo(i), i = 1, j)
    allocate(mrgp%num_tot_gal(j), stat=ier); call CheckIerr(ier, &
         trim(cerr)//' mrgp%num_tot_gal')
    mrgp%num_tot_gal(:) = 0


    ! --- allocate and read the quantities in the structure named "mrgt"
    j = mrgt%num_tot
    allocate(mrgt%f_des(0:j), stat=ier); call CheckIerr(ier, &
         trim(cerr)//' mrgt%f_des')
    read(i_file) (mrgt%f_des(i), i = 0, j)
    allocate(mrgt%n_des(0:j), stat=ier); call CheckIerr(ier, &
         trim(cerr)//' mrgt%n_des')
    read(i_file) (mrgt%n_des(i), i = 0, j)
    allocate(mrgt%f_prg(0:j), stat=ier); call CheckIerr(ier, &
         trim(cerr)//' mrgt%f_prg')
    read(i_file) (mrgt%f_prg(i), i = 0, j)
    allocate(mrgt%numb(0:j), stat=ier); call CheckIerr(ier, &
         trim(cerr)//' mrgt%numb')
    read(i_file) (mrgt%numb(i),  i = 0, j)
    allocate(mrgt%hori(0:j), stat=ier); call CheckIerr(ier, &
         trim(cerr)//' mrgt%hori')
    read(i_file) (mrgt%hori(i),  i = 0, j)
    allocate(mrgt%mpi(0:j), stat=ier); call CheckIerr(ier, &
         trim(cerr)//' mrgt%mpi')
    read(i_file) (mrgt%mpi(i),   i = 0, j)
    allocate(mrgt%mhalo(0:j), stat=ier); call CheckIerr(ier, &
         trim(cerr)//' mrgt%mhalo')
    read(i_file) (mrgt%mhalo(i), i = 0, j)
    mrgt%mhalo(:) = mrgt%mhalo(:) / param%munit ! [Msun] --> [10^14 Msun]
    allocate(mrgt%hstart(0:j),  stat=ier); call CheckIerr(ier, &
         trim(cerr)//' mrgt%hstart')
    mrgt%hstart(:) = 0
    allocate(mrgt%c_halo(0:j),  stat=ier); call CheckIerr(ier, &
         trim(cerr)//' mrgt%c_halo')
    mrgt%c_halo(:) = 0
    allocate(mrgt%c2_halo(0:j), stat=ier); call CheckIerr(ier, &
         trim(cerr)//' mrgt%c2_halo')
    mrgt%c2_halo(:) = 0
    allocate(mrgt%num_g(0:j),   stat=ier); call CheckIerr(ier, &
         trim(cerr)//' mrgt%num_g')
    mrgt%num_g(:) = 0
    close(i_file)

    IF(inode == 0 .and. nz == 1 .and. param%run_type /= 3) &
       print '(A, I7, A, I3)', '            --- Finish: num_tot = ', &
         mrgt%num_tot, ', num_step = ', mrgp%num_step

    mrgt%num_galarray = int(mrgt%num_tot * param%f_ngal) ! added by MARK (2016/Jun/29)
!!$    mrgt%num_galarray = mrgt%num_tot
!!$    mrgt%num_galarray = mrgt%num_tot / 10
!!$    mrgt%num_galarray = mrgt%num_tot / 5 ! for 140 Mpc/h w/ 4096^3
!!$    mrgt%num_galarray = mrgt%num_tot / 20
  END SUBROUTINE ReadMRGB
!!$============================================================================
  SUBROUTINE CheckNumGtot(endhalo, end_step)
    INTEGER, INTENT(IN) :: endhalo, end_step
    INTEGER :: num_g_tot

    num_g_tot = 0
    DO i = 1, mrgp%num_now(end_step)
       me = mrgp%st_halo(end_step) + i - 1
       num_g_tot = num_g_tot + mrgt%num_g(me)
    ENDDO
    print '(4(A, I8))', '# endhalo = ', endhalo, ', num_g_tot = ', num_g_tot, &
         ': num_tot_gal = ', mrgp%num_tot_gal(end_step), &
         ', end_step = ', end_step
  END SUBROUTINE CheckNumGtot
!!$============================================================================
  SUBROUTINE DeallocMRGB
!!$    INTEGER :: i
    INTEGER :: ier
    CHARACTER(LEN=50) :: cerr = '# DeallocMRGB: fail to deallocate '


    ! --- mrgp
    deallocate(mrgp%time,        stat=ier); call CheckIerr(ier, &
         trim(cerr)//' mrgp%time')
    deallocate(mrgp%zp1ar,       stat=ier); call CheckIerr(ier, &
         trim(cerr)//' mrgp%zp1ar')
    deallocate(mrgp%num_now,     stat=ier); call CheckIerr(ier, &
         trim(cerr)//' mrgp%num_now')
    deallocate(mrgp%st_halo,     stat=ier); call CheckIerr(ier, &
         trim(cerr)//' mrgp%st_halo')
    deallocate(mrgp%num_tot_gal, stat=ier); call CheckIerr(ier, &
         trim(cerr)//' mrgp%num_tot_gal')

    ! --- mrgt
    deallocate(mrgt%f_des,   stat=ier); call CheckIerr(ier, &
         trim(cerr)//' mrgt%f_des')
    deallocate(mrgt%n_des,   stat=ier); call CheckIerr(ier, &
         trim(cerr)//' mrgt%n_des')
    deallocate(mrgt%f_prg,   stat=ier); call CheckIerr(ier, &
         trim(cerr)//' mrgt%f_prg')
    deallocate(mrgt%numb,    stat=ier); call CheckIerr(ier, &
         trim(cerr)//' mrgt%numb')
    deallocate(mrgt%hori,    stat=ier); call CheckIerr(ier, &
         trim(cerr)//' mrgt%hori')
    deallocate(mrgt%mpi,     stat=ier); call CheckIerr(ier, &
         trim(cerr)//' mrgt%mpi')
    deallocate(mrgt%mhalo,   stat=ier); call CheckIerr(ier, &
         trim(cerr)//' mrgt%mhalo')
    deallocate(mrgt%hstart,  stat=ier); call CheckIerr(ier, &
         trim(cerr)//' mrgt%hstart')
    deallocate(mrgt%c_halo,  stat=ier); call CheckIerr(ier, &
         trim(cerr)//' mrgt%c_halo')
    deallocate(mrgt%c2_halo, stat=ier); call CheckIerr(ier, &
         trim(cerr)//' mrgt%c2_halo')
    deallocate(mrgt%num_g,   stat=ier); call CheckIerr(ier, &
         trim(cerr)//' mrgt%num_g')
  END SUBROUTINE DeallocMRGB
!!$============================================================================
  SUBROUTINE AllocateGalaxy
    use global_var
    implicit none
    INTEGER :: i, ngal
    CHARACTER*50 :: cerr = '# AllocateGalaxy: fail to allocate '


    ngal = mrgt%num_galarray
    allocate(gal(ngal),      stat=ier); call CheckIerr(ier, &
         trim(cerr)//' gal')
    allocate(gal_prev(ngal), stat=ier); call CheckIerr(ier, &
         trim(cerr)//' gal_prev')
    allocate(gal_next(ngal), stat=ier); call CheckIerr(ier, &
         trim(cerr)//' gal_next')
    DO i = 1, ngal
       ! --- allocate lumg in the structures of gal, gal_prev, and gal_next
       allocate(gal(i)%lumg(param%tnw),      stat=ier); call CheckIerr(ier, &
            trim(cerr)//' gal(i)%lumg')
       allocate(gal_prev(i)%lumg(param%tnw), stat=ier); call CheckIerr(ier, &
            trim(cerr)//' gal_prev(i)%lumg')
       allocate(gal_next(i)%lumg(param%tnw), stat=ier); call CheckIerr(ier, &
            trim(cerr)//' gal_next(i)%lumg')

       ! --- allocate lumq in the structures of gal, gal_prev, and gal_next
       allocate(gal(i)%lumq(param%nwaveq),      stat=ier); call CheckIerr(ier, &
            trim(cerr)//' gal(i)%lumq')
       allocate(gal_prev(i)%lumq(param%nwaveq), stat=ier); call CheckIerr(ier, &
            trim(cerr)//' gal_prev(i)%lumq')
       allocate(gal_next(i)%lumq(param%nwaveq), stat=ier); call CheckIerr(ier, &
            trim(cerr)//' gal_next(i)%lumq')
              ! --- allocate agn in the structures of gal, gal_prev, and gal_next
       allocate(gal(i)%agn(10),      stat=ier); call CheckIerr(ier, &
            trim(cerr)//' gal(i)%agn')
       allocate(gal_prev(i)%agn(10), stat=ier); call CheckIerr(ier, &
            trim(cerr)//' gal_prev(i)%agn')
       allocate(gal_next(i)%agn(10), stat=ier); call CheckIerr(ier, &
            trim(cerr)//' gal_next(i)%agn')
       IF(param%traceIDs == 1) THEN
          allocate(gal(i)%mem(10), stat=ier); call CheckIerr(ier, &
             trim(cerr)//' gal(i)%mem')
       ENDIF
    ENDDO
    call InitializeTypeGalaxy(gal,      mrgt%num_galarray)
    call InitializeTypeGalaxy(gal_prev, mrgt%num_galarray)
    call InitializeTypeGalaxy(gal_next, mrgt%num_galarray)
  END SUBROUTINE AllocateGalaxy
!!$============================================================================
  SUBROUTINE InitializeTypeGalaxy(gal, numg)
!!$    use global_var
    implicit none
    INTEGER, INTENT(IN) :: numg
    TYPE(galaxy), INTENT(INOUT) :: gal(numg)
    INTEGER :: igal


    DO igal = 1, numg
       gal(igal)%flag_burst = 0
       gal(igal)%flag_di    = 0
       gal(igal)%flag_merger= 0
       gal(igal)%flag_seed  = 0
       gal(igal)%IDhost     = 0
       gal(igal)%flag_c     = 0
       gal(igal)%hori       = 0
       gal(igal)%mpi        = 0
       gal(igal)%hstart     = 0
       gal(igal)%hfinal     = 0
       gal(igal)%id_cgal    = 0

       gal(igal)%Mhalo = 0.d0; gal(igal)%Mhot = 0.d0; gal(igal)%Mhotout = 0.d0
       gal(igal)%MZh = 0.d0
       gal(igal)%Vc = 0.d0
       gal(igal)%Mreheat = 0.d0
       gal(igal)%Telapse = 0.d0
       gal(igal)%Mratio  = 0.d0
       gal(igal)%Morg    = 0.d0; gal(igal)%Mhotorg = 0.d0; gal(igal)%MZhorg = 0.d0
       gal(igal)%Mstarb  = 0.d0; gal(igal)%MZb = 0.d0; gal(igal)%Mtb = 0.d0
       gal(igal)%Mstard  = 0.d0; gal(igal)%MZd = 0.d0; gal(igal)%Mtd = 0.d0
       gal(igal)%Mcoold   = 0.d0; gal(igal)%MZcd = 0.d0
       gal(igal)%Mcoolb   = 0.d0; gal(igal)%MZcb = 0.d0
       gal(igal)%Vbulge  = 0.d0; gal(igal)%rbulge = 0.d0
       gal(igal)%Vdisk   = 0.d0; gal(igal)%rdisk  = 0.d0
       gal(igal)%diskmas = 0.d0
       gal(igal)%density = 0.d0
       gal(igal)%clps    = 0.d0
       gal(igal)%Vcent   = 0.d0
       gal(igal)%MZc_rem = 0.d0
       gal(igal)%SFR     = 0.d0
       gal(igal)%mSFR    = 0.d0 ! added by MARK (2014/Nov/05)

       gal(igal)%MZg = 0.d0; gal(igal)%Mtg = 0.d0

       gal(igal)%lumg(:) = 0.d0

       gal(igal)%Zmassb = 0.d0; gal(igal)%Zmassd = 0.d0 ! added by MARK (2014/Nov/05)

       ! --- added by MARK (2017/Mar/16)
       gal(igal)%taust = 0.d0; gal(igal)%beta = 0.d0
       gal(igal)%dMstar_burst = 0.d0

       ! --- AGN related ---
       gal(igal)%Mbh = 0.d0; gal(igal)%lumq(:) = 0.d0
       gal(igal)%agn(:) = 0.d0

       ! --- begining of star formation redshift (by MN)
!!$         gal(igal)%z_sf = 0.d0
       gal(igal)%z_sf = 0

       ! --- add by MARK ---
       gal(igal)%Tmass = 0.d0; gal(igal)%Tlum_b = 0.d0; gal(igal)%Tlum_d = 0.d0
   ENDDO
  END SUBROUTINE InitializeTypeGalaxy
!!$============================================================================
  SUBROUTINE DeallocateGalaxy
    use global_var
    implicit none
    INTEGER :: i, ier
    CHARACTER*50 :: cerr = '# DeallocateGalaxy: fail to deallocate '

    DO i = 1, mrgt%num_galarray
       ! --- deallocate lumg in the structures of gal, gal_prev, and gal_next
       deallocate(gal(i)%lumg, stat=ier); call CheckIerr(ier, &
          trim(cerr)//' gal(i)%lumg')
       deallocate(gal_prev(i)%lumg, stat=ier); call CheckIerr(ier, &
          trim(cerr)//' gal_prev(i)%lumg')
       deallocate(gal_next(i)%lumg, stat=ier); call CheckIerr(ier, &
          trim(cerr)//' gal_next(i)%lumg')

       ! --- deallocate lumq in the structures of gal, gal_prev, and gal_next
       deallocate(gal(i)%lumq, stat=ier); call CheckIerr(ier, &
          trim(cerr)//' gal(i)%lumq')
       deallocate(gal_prev(i)%lumq, stat=ier); call CheckIerr(ier, &
          trim(cerr)//' gal_prev(i)%lumq')
       deallocate(gal_next(i)%lumq, stat=ier); call CheckIerr(ier, &
          trim(cerr)//' gal_next(i)%lumq')
       ! --- deallocate agn in the structures of gal, gal_prev, and gal_next
       deallocate(gal(i)%agn, stat=ier); call CheckIerr(ier, &
          trim(cerr)//' gal(i)%agn')
       deallocate(gal_prev(i)%agn, stat=ier); call CheckIerr(ier, &
          trim(cerr)//' gal_prev(i)%agn')
       deallocate(gal_next(i)%agn, stat=ier); call CheckIerr(ier, &
          trim(cerr)//' gal_next(i)%agn')
       IF(param%traceIDs == 1) THEN
          deallocate(gal(i)%mem, stat=ier); call CheckIerr(ier, &
             trim(cerr)//' gal(i)%mem')
       ENDIF
    ENDDO

    deallocate(gal,      stat=ier); call CheckIerr(ier, trim(cerr)//' gal')
    deallocate(gal_prev, stat=ier); call CheckIerr(ier, trim(cerr)//' gal_prev')
    deallocate(gal_next, stat=ier); call CheckIerr(ier, trim(cerr)//' gal_next')
  END SUBROUTINE DeallocateGalaxy
!!$============================================================================
  SUBROUTINE DealloctraceIDs
     INTEGER :: ier
     CHARACTER(LEN=50) :: cerr = '# DealloctraceIDs: fail to deallocate '

     deallocate(targ, stat=ier);  call CheckIerr(ier, trim(cerr)//'targ')
  END SUBROUTINE DealloctraceIDs
!!$============================================================================
  SUBROUTINE DeallocSSP
    INTEGER :: ier
    CHARACTER(LEN=50) :: cerr = '# DeallocSSP: fail to deallocate '

    deallocate(ssp%time, stat=ier); call CheckIerr(ier, trim(cerr)//'ssp%time')
    deallocate(ssp%chem, stat=ier); call CheckIerr(ier, trim(cerr)//'ssp%chem')
    deallocate(ssp%lumi, stat=ier); call CheckIerr(ier, trim(cerr)//'ssp%lumi')
  END SUBROUTINE DeallocSSP
!!$============================================================================
  SUBROUTINE InitializeGlobalFilterIntegers
    use global_var

    param%iBband  = 0; param%iUband  = 0; param%iVband  = 0
    param%iRcband = 0; param%iRJband = 0; param%iIcband = 0
    param%iIJband = 0; param%izband  = 0; param%iJband  = 0
    param%iHband  = 0; param%iKband  = 0; param%iKpband = 0
    param%iLband  = 0
    param%iSuprime_Bband = 0; param%iSuprime_gband  = 0; param%iSuprime_Vband  = 0
    param%iSuprime_rband = 0; param%iSuprime_RLband = 0; param%iSuprime_ipband = 0
    param%iSuprime_Iband = 0; param%iSuprime_zpband = 0; param%iSuprime_zp_redband = 0
    param%iobs5um = 0
    param%i2MASS_Jband = 0; param%i2MASS_Hband = 0; param%i2MASS_Ksband = 0
    param%iACS_F775Wband = 0; param%iACS_F850LPband = 0
    param%iCISCO_zband = 0; param%iCISCO_Jband  = 0; param%iCISCO_Hband = 0
    param%iCISCO_Kband = 0; param%iCISCO_Kpband = 0
    param%iGALEX_FUVband = 0; param%iGALEX_NUVband = 0
    param%iGOODS_RBband = 0; param%iGOODS_RGband = 0; param%iGOODS_RSband = 0
    param%iHST_F300wband = 0; param%iHST_F450wband = 0; param%iHST_F555wband = 0
    param%iHST_F606wband = 0; param%iHST_F702wband = 0; param%iHST_F814wband = 0
    param%iSDSS_upband = 0; param%iSDSS_gpband = 0; param%iSDSS_rpband = 0
    param%iSDSS_ipband = 0; param%iSDSS_zpband = 0
    param%iNLyC  = 0; param%iL1216 = 0; param%iL1400 = 0; param%iL1500 = 0
    param%iL1600 = 0; param%iL1700 = 0; param%iL2800 = 0; param%iL4861 = 0
    param%iL6563 = 0
    param%iIRACch1band = 0; param%iIRACch2band = 0; param%iIRACch3band = 0
    param%iIRACch4band = 0
    param%iWFCAM_zband = 0; param%iWFCAM_Yband = 0; param%iWFCAM_Jband = 0
    param%iWFCAM_Hband = 0; param%iWFCAM_Kband = 0
    param%iVIRCAM_Yband = 0; param%iVIRCAM_Jband = 0
    param%iVIRCAM_Hband = 0; param%iVIRCAM_Kband = 0
    param%iCIBER_Iband = 0; param%iCIBER_Hband = 0
    param%iNEWFIRM_J1band = 0; param%iNEWFIRM_J2band = 0; param%iNEWFIRM_J3band = 0
    param%iNEWFIRM_H1band = 0; param%iNEWFIRM_H2band = 0; param%iNEWFIRM_Ksband = 0
    param%iAKARI_N2band = 0; param%iAKARI_N3band  = 0; param%iAKARI_N4band  = 0
    param%iAKARI_S7band = 0; param%iAKARI_S9Wband = 0; param%iAKARI_S11band = 0
    param%iWISH_0band = 0; param%iWISH_1band = 0; param%iWISH_2band = 0
    param%iWISH_3band = 0; param%iWISH_4band = 0; param%iWISH_5band = 0
    param%iACS_F435Wband   = 0; param%iACS_F475Wband    = 0; param%iACS_F606Wband   = 0
    param%iACS_F814Wband   = 0; param%iNICMOS_F160Wband = 0; param%iWFPC2_F300Wband = 0
    param%iWFPC2_F450Wband = 0; param%iWFPC2_F555Wband  = 0; param%iWFC3_F125Wband  = 0
    param%iWFC3_F140Wband  = 0; param%iWFC3_F160Wband   = 0
    param%iHSC_gband = 0; param%iHSC_rband = 0; param%iHSC_iband = 0
    param%iHSC_zband = 0; param%iHSC_Yband = 0

    param%iBband_r  = 0; param%iUband_r  = 0; param%iVband_r  = 0
    param%iRcband_r = 0; param%iRJband_r = 0; param%iIcband_r = 0
    param%iIJband_r = 0; param%izband_r  = 0; param%iJband_r  = 0
    param%iHband_r  = 0; param%iKband_r  = 0; param%iKpband_r = 0
    param%iLband_r  = 0
    param%iSuprime_Bband_r  = 0; param%iSuprime_gband_r  = 0
    param%iSuprime_Vband_r  = 0; param%iSuprime_rband_r  = 0
    param%iSuprime_RLband_r = 0; param%iSuprime_ipband_r = 0
    param%iSuprime_Iband_r  = 0; param%iSuprime_zpband_r = 0
    param%iSuprime_zp_redband_r = 0
    param%iobs5um = 0
    param%i2MASS_Jband_r = 0; param%i2MASS_Hband_r = 0
    param%i2MASS_Ksband_r = 0
    param%iACS_F775Wband_r = 0; param%iACS_F850LPband_r = 0
    param%iCISCO_zband_r = 0; param%iCISCO_Jband_r = 0
    param%iCISCO_Hband_r = 0; param%iCISCO_Kband_r = 0
    param%iCISCO_Kpband_r = 0
    param%iGALEX_FUVband_r = 0; param%iGALEX_NUVband_r = 0
    param%iGOODS_RBband_r = 0; param%iGOODS_RGband_r = 0
    param%iGOODS_RSband_r = 0
    param%iHST_F300wband_r = 0; param%iHST_F450wband_r = 0
    param%iHST_F555wband_r = 0; param%iHST_F606wband_r = 0
    param%iHST_F702wband_r = 0; param%iHST_F814wband_r = 0
    param%iSDSS_upband_r = 0; param%iSDSS_gpband_r = 0
    param%iSDSS_rpband_r = 0; param%iSDSS_ipband_r = 0
    param%iSDSS_zpband_r = 0
    param%iNLyC  = 0; param%iL1216 = 0; param%iL1400 = 0
    param%iL1500 = 0; param%iL1600 = 0; param%iL1700 = 0
    param%iL2800 = 0; param%iL4861 = 0; param%iL6563 = 0
    param%iIRACch1band_r = 0; param%iIRACch2band_r = 0
    param%iIRACch3band_r = 0; param%iIRACch4band_r = 0
    param%iWFCAM_zband_r = 0; param%iWFCAM_Yband_r = 0
    param%iWFCAM_Jband_r = 0; param%iWFCAM_Hband_r = 0
    param%iWFCAM_Kband_r = 0
    param%iVIRCAM_Yband_r = 0; param%iVIRCAM_Jband_r = 0
    param%iVIRCAM_Hband_r = 0; param%iVIRCAM_Kband_r = 0
    param%iCIBER_Iband_r = 0; param%iCIBER_Hband_r = 0
    param%iNEWFIRM_J1band_r = 0; param%iNEWFIRM_J2band_r = 0
    param%iNEWFIRM_J3band_r = 0; param%iNEWFIRM_H1band_r = 0
    param%iNEWFIRM_H2band_r = 0; param%iNEWFIRM_Ksband_r = 0
    param%iAKARI_N2band_r = 0; param%iAKARI_N3band_r = 0
    param%iAKARI_N4band_r = 0; param%iAKARI_S7band_r = 0
    param%iAKARI_S9Wband_r = 0; param%iAKARI_S11band_r = 0
    param%iWISH_0band_r = 0; param%iWISH_1band_r = 0
    param%iWISH_2band_r = 0; param%iWISH_3band_r = 0
    param%iWISH_4band_r = 0; param%iWISH_5band_r = 0
    param%iACS_F435Wband_r    = 0; param%iACS_F475Wband_r   = 0
    param%iACS_F606Wband_r    = 0; param%iACS_F814Wband_r   = 0
    param%iNICMOS_F160Wband_r = 0; param%iWFPC2_F300Wband_r = 0
    param%iWFPC2_F450Wband_r  = 0; param%iWFPC2_F555Wband_r = 0
    param%iWFC3_F125Wband_r   = 0; param%iWFC3_F140Wband_r = 0
    param%iWFC3_F160Wband_r   = 0
    param%iHSC_gband_r = 0; param%iHSC_rband_r = 0; param%iHSC_iband_r = 0
    param%iHSC_zband_r = 0; param%iHSC_Yband_r = 0

    param%wdim_band = .false. ! added by MARK (2018/Mar/02)
  END SUBROUTINE InitializeGlobalFilterIntegers
!!$============================================================================
  SUBROUTINE SubstGlobalFilterIntegers(i_filter)
    ! --- substitute the integer of 'i_filter' to the corresponding integer
    !      in global_var of 'param%i**band' or 'param%i**band_r'
    use global_var

    INTEGER, INTENT(IN) :: i_filter
    INTEGER :: i_cons 

    i_cons = param%iwave(i_filter)
             ! ---  the consecutive number related to the bandname of
             !       gal(*)%lumg(i_filter)

    IF(i_cons <= NWAVE_ALL) THEN
       IF(i_cons <= 10) THEN
          IF(i_cons == 1) THEN
             param%iBband = i_filter
          ELSEIF(i_cons == 2) THEN
             param%iUband = i_filter
          ELSEIF(i_cons == 3) THEN
             param%iVband = i_filter
          ELSEIF(i_cons == 4) THEN
             param%iRcband = i_filter
          ELSEIF(i_cons == 5) THEN
             param%iRJband = i_filter
          ELSEIF(i_cons == 6) THEN
             param%iIcband = i_filter
          ELSEIF(i_cons == 7) THEN
             param%iIJband = i_filter
          ELSEIF(i_cons == 8) THEN
             param%izband = i_filter
          ELSEIF(i_cons == 9) THEN
             param%iJband = i_filter
          ELSEIF(i_cons == 10) THEN
             param%iHband = i_filter
          ENDIF
       ELSEIF(i_cons <= 20) THEN
          IF(i_cons == 11) THEN
             param%iKband = i_filter
          ELSEIF(i_cons == 12) THEN
             param%iKpband = i_filter
          ELSEIF(i_cons == 13) THEN
             param%iLband = i_filter
          ELSEIF(i_cons == 14) THEN
             param%iSuprime_Bband = i_filter
          ELSEIF(i_cons == 15) THEN
             param%iSuprime_gband = i_filter
          ELSEIF(i_cons == 16) THEN
             param%iSuprime_Vband = i_filter
          ELSEIF(i_cons == 17) THEN
             param%iSuprime_rband = i_filter
          ELSEIF(i_cons == 18) THEN
             param%iSuprime_RLband = i_filter
          ELSEIF(i_cons == 19) THEN
             param%iSuprime_ipband = i_filter
          ELSEIF(i_cons == 20) THEN
             param%iSuprime_Iband = i_filter
          ENDIF
       ELSEIF(i_cons <= 30) THEN
          IF(i_cons == 21) THEN
             param%iSuprime_zpband = i_filter
          ELSEIF(i_cons == 22) THEN
             param%iSuprime_zp_redband = i_filter
          ELSEIF(i_cons == 23) THEN
             param%iobs5um = i_filter
          ELSEIF(i_cons == 24) THEN
             param%i2MASS_Jband = i_filter
          ELSEIF(i_cons == 25) THEN
             param%i2MASS_Hband = i_filter
          ELSEIF(i_cons == 26) THEN
             param%i2MASS_Ksband = i_filter
          ELSEIF(i_cons == 27) THEN
             param%iACS_F775Wband = i_filter
          ELSEIF(i_cons == 28) THEN
             param%iACS_F850LPband = i_filter
          ELSEIF(i_cons == 29) THEN
             param%iCISCO_zband = i_filter
          ELSEIF(i_cons == 30) THEN
             param%iCISCO_Jband = i_filter
          ENDIF
       ELSEIF(i_cons <= 40) THEN
          IF(i_cons == 31) THEN
             param%iCISCO_Hband = i_filter
          ELSEIF(i_cons == 32) THEN
             param%iCISCO_Kband = i_filter
          ELSEIF(i_cons == 33) THEN
             param%iCISCO_Kpband = i_filter
          ELSEIF(i_cons == 34) THEN
             param%iGALEX_FUVband = i_filter
          ELSEIF(i_cons == 35) THEN
             param%iGALEX_NUVband = i_filter
          ELSEIF(i_cons == 36) THEN
             param%iGOODS_RBband = i_filter
          ELSEIF(i_cons == 37) THEN
             param%iGOODS_RGband = i_filter
          ELSEIF(i_cons == 38) THEN
             param%iGOODS_RSband = i_filter
          ELSEIF(i_cons == 39) THEN
             param%iHST_F300wband = i_filter
          ELSEIF(i_cons == 40) THEN
             param%iHST_F450wband = i_filter
          ENDIF
       ELSEIF(i_cons <= 50) THEN
          IF(i_cons == 41) THEN
             param%iHST_F555wband = i_filter
          ELSEIF(i_cons == 42) THEN
             param%iHST_F606wband = i_filter
          ELSEIF(i_cons == 43) THEN
             param%iHST_F702wband = i_filter
          ELSEIF(i_cons == 44) THEN
             param%iHST_F814wband = i_filter
          ELSEIF(i_cons == 45) THEN
             param%iSDSS_upband = i_filter
          ELSEIF(i_cons == 46) THEN
             param%iSDSS_gpband = i_filter
          ELSEIF(i_cons == 47) THEN
             param%iSDSS_rpband = i_filter
          ELSEIF(i_cons == 48) THEN
             param%iSDSS_ipband = i_filter
          ELSEIF(i_cons == 49) THEN
             param%iSDSS_zpband = i_filter
          ELSEIF(i_cons == 50) THEN
             param%iNLyC = i_filter
          ENDIF
       ELSEIF(i_cons <= 60) THEN
          IF(i_cons == 51) THEN
             param%iL1216 = i_filter
          ELSEIF(i_cons == 52) THEN
             param%iL1400 = i_filter
          ELSEIF(i_cons == 53) THEN
             param%iL1500 = i_filter
          ELSEIF(i_cons == 54) THEN
             param%iL1600 = i_filter
          ELSEIF(i_cons == 55) THEN
             param%iL1700 = i_filter
          ELSEIF(i_cons == 56) THEN
             param%iL2800 = i_filter
          ELSEIF(i_cons == 57) THEN
             param%iL4861 = i_filter
          ELSEIF(i_cons == 58) THEN
             param%iL6563 = i_filter
          ELSEIF(i_cons == 59) THEN
             param%iIRACch1band = i_filter
          ELSEIF(i_cons == 60) THEN
             param%iIRACch2band = i_filter
          ENDIF
       ELSEIF(i_cons <= 70) THEN
          IF(i_cons == 61) THEN
             param%iIRACch3band = i_filter
          ELSEIF(i_cons == 62) THEN
             param%iIRACch4band = i_filter
          ELSEIF(i_cons == 63) THEN
             param%iWFCAM_zband = i_filter
          ELSEIF(i_cons == 64) THEN
             param%iWFCAM_Yband = i_filter
          ELSEIF(i_cons == 65) THEN
             param%iWFCAM_Jband = i_filter
          ELSEIF(i_cons == 66) THEN
             param%iWFCAM_Hband = i_filter
          ELSEIF(i_cons == 67) THEN
             param%iWFCAM_Kband = i_filter
          ELSEIF(i_cons == 68) THEN
             param%iVIRCAM_Yband = i_filter
          ELSEIF(i_cons == 69) THEN
             param%iVIRCAM_Jband = i_filter
          ELSEIF(i_cons == 70) THEN
             param%iVIRCAM_Hband = i_filter
          ENDIF
       ELSEIF(i_cons <= 80) THEN
          IF(i_cons == 71) THEN
             param%iVIRCAM_Kband = i_filter
          ELSEIF(i_cons == 72) THEN
             param%iCIBER_Iband = i_filter
          ELSEIF(i_cons == 73) THEN
             param%iCIBER_Hband = i_filter
          ELSEIF(i_cons == 74) THEN
             param%iNEWFIRM_J1band = i_filter
          ELSEIF(i_cons == 75) THEN
             param%iNEWFIRM_J2band = i_filter
          ELSEIF(i_cons == 76) THEN
             param%iNEWFIRM_J3band = i_filter
          ELSEIF(i_cons == 77) THEN
             param%iNEWFIRM_H1band = i_filter
          ELSEIF(i_cons == 78) THEN
             param%iNEWFIRM_H2band = i_filter
          ELSEIF(i_cons == 79) THEN
             param%iNEWFIRM_Ksband = i_filter
          ELSEIF(i_cons == 80) THEN
             param%iAKARI_N2band = i_filter
          ENDIF
       ELSEIF(i_cons <= 90) THEN
          IF(i_cons == 81) THEN
             param%iAKARI_N3band = i_filter
          ELSEIF(i_cons == 82) THEN
             param%iAKARI_N4band = i_filter
          ELSEIF(i_cons == 83) THEN
             param%iAKARI_S7band = i_filter
          ELSEIF(i_cons == 84) THEN
             param%iAKARI_S9Wband = i_filter
          ELSEIF(i_cons == 85) THEN
             param%iAKARI_S11band = i_filter
          ELSEIF(i_cons == 86) THEN
             param%iWISH_0band = i_filter
          ELSEIF(i_cons == 87) THEN
             param%iWISH_1band = i_filter
          ELSEIF(i_cons == 88) THEN
             param%iWISH_2band = i_filter
          ELSEIF(i_cons == 89) THEN
             param%iWISH_3band = i_filter
          ELSEIF(i_cons == 90) THEN
             param%iWISH_4band = i_filter
          ENDIF
       ELSEIF(i_cons <= 100) THEN
          IF(i_cons == 91) THEN
             param%iWISH_5band = i_filter
          ELSEIF(i_cons == 92) THEN
             param%iACS_F435Wband = i_filter
          ELSEIF(i_cons == 93) THEN
             param%iACS_F475Wband = i_filter
          ELSEIF(i_cons == 94) THEN
             param%iACS_F606Wband = i_filter
          ELSEIF(i_cons == 95) THEN
             param%iACS_F814Wband = i_filter
          ELSEIF(i_cons == 96) THEN
             param%iNICMOS_F160Wband = i_filter
          ELSEIF(i_cons == 97) THEN
             param%iWFPC2_F300Wband = i_filter
          ELSEIF(i_cons == 98) THEN
             param%iWFPC2_F450Wband = i_filter
          ELSEIF(i_cons == 99) THEN
             param%iWFPC2_F555Wband = i_filter
          ELSEIF(i_cons == 100) THEN
             param%iWFC3_F125Wband = i_filter
          ENDIF
       ELSEIF(i_cons <= 110) THEN
          IF(i_cons == 101) THEN
             param%iWFC3_F140Wband = i_filter
          ELSEIF(i_cons == 102) THEN
             param%iWFC3_F160Wband = i_filter
          ELSEIF(i_cons == 103) THEN
             param%iHSC_gband = i_filter
          ELSEIF(i_cons == 104) THEN
             param%iHSC_rband = i_filter
          ELSEIF(i_cons == 105) THEN
             param%iHSC_iband = i_filter
          ELSEIF(i_cons == 106) THEN
             param%iHSC_zband = i_filter
          ELSEIF(i_cons == 107) THEN
             param%iHSC_Yband = i_filter
          ENDIF
       ENDIF
    ELSE ! i_cons > NWAVE_ALL
       IF(i_cons <= NWAVE_ALL+10) THEN
          IF(i_cons == NWAVE_ALL+1) THEN
             param%iBband_r = i_filter
          ELSEIF(i_cons == NWAVE_ALL+2) THEN
             param%iUband_r = i_filter
          ELSEIF(i_cons == NWAVE_ALL+3) THEN
             param%iVband_r = i_filter
          ELSEIF(i_cons == NWAVE_ALL+4) THEN
             param%iRcband_r = i_filter
          ELSEIF(i_cons == NWAVE_ALL+5) THEN
             param%iRJband_r = i_filter
          ELSEIF(i_cons == NWAVE_ALL+6) THEN
             param%iIcband_r = i_filter
          ELSEIF(i_cons == NWAVE_ALL+7) THEN
             param%iIJband_r = i_filter
          ELSEIF(i_cons == NWAVE_ALL+8) THEN
             param%izband_r = i_filter
          ELSEIF(i_cons == NWAVE_ALL+9) THEN
             param%iJband_r = i_filter
          ELSEIF(i_cons == NWAVE_ALL+10) THEN
             param%iHband_r = i_filter
          ENDIF
       ELSEIF(i_cons <= NWAVE_ALL+20) THEN
          IF(i_cons == NWAVE_ALL+11) THEN
             param%iKband_r = i_filter
          ELSEIF(i_cons == NWAVE_ALL+12) THEN
             param%iKpband_r = i_filter
          ELSEIF(i_cons == NWAVE_ALL+13) THEN
             param%iLband_r = i_filter
          ELSEIF(i_cons == NWAVE_ALL+14) THEN
             param%iSuprime_Bband_r = i_filter
          ELSEIF(i_cons == NWAVE_ALL+15) THEN
             param%iSuprime_gband_r = i_filter
          ELSEIF(i_cons == NWAVE_ALL+16) THEN
             param%iSuprime_Vband_r = i_filter
          ELSEIF(i_cons == NWAVE_ALL+17) THEN
             param%iSuprime_rband_r = i_filter
          ELSEIF(i_cons == NWAVE_ALL+18) THEN
             param%iSuprime_RLband_r = i_filter
          ELSEIF(i_cons == NWAVE_ALL+19) THEN
             param%iSuprime_ipband_r = i_filter
          ELSEIF(i_cons == NWAVE_ALL+20) THEN
             param%iSuprime_Iband_r = i_filter
          ENDIF
       ELSEIF(i_cons <= NWAVE_ALL+30) THEN
          IF(i_cons == NWAVE_ALL+21) THEN
             param%iSuprime_zpband_r = i_filter
          ELSEIF(i_cons == NWAVE_ALL+22) THEN
             param%iSuprime_zp_redband_r = i_filter
          ELSEIF(i_cons == NWAVE_ALL+23) THEN
             param%iobs5um_r = i_filter
          ELSEIF(i_cons == NWAVE_ALL+24) THEN
             param%i2MASS_Jband_r = i_filter
          ELSEIF(i_cons == NWAVE_ALL+25) THEN
             param%i2MASS_Hband_r = i_filter
          ELSEIF(i_cons == NWAVE_ALL+26) THEN
             param%i2MASS_Ksband_r = i_filter
          ELSEIF(i_cons == NWAVE_ALL+27) THEN
             param%iACS_F775Wband_r = i_filter
          ELSEIF(i_cons == NWAVE_ALL+28) THEN
             param%iACS_F850LPband_r = i_filter
          ELSEIF(i_cons == NWAVE_ALL+29) THEN
             param%iCISCO_zband_r = i_filter
          ELSEIF(i_cons == NWAVE_ALL+30) THEN
             param%iCISCO_Jband_r = i_filter
          ENDIF
       ELSEIF(i_cons <= NWAVE_ALL+40) THEN
          IF(i_cons == NWAVE_ALL+31) THEN
             param%iCISCO_Hband_r = i_filter
          ELSEIF(i_cons == NWAVE_ALL+32) THEN
             param%iCISCO_Kband_r = i_filter
          ELSEIF(i_cons == NWAVE_ALL+33) THEN
             param%iCISCO_Kpband_r = i_filter
          ELSEIF(i_cons == NWAVE_ALL+34) THEN
             param%iGALEX_FUVband_r = i_filter
          ELSEIF(i_cons == NWAVE_ALL+35) THEN
             param%iGALEX_NUVband_r = i_filter
          ELSEIF(i_cons == NWAVE_ALL+36) THEN
             param%iGOODS_RBband_r = i_filter
          ELSEIF(i_cons == NWAVE_ALL+37) THEN
             param%iGOODS_RGband_r = i_filter
          ELSEIF(i_cons == NWAVE_ALL+38) THEN
             param%iGOODS_RSband_r = i_filter
          ELSEIF(i_cons == NWAVE_ALL+39) THEN
             param%iHST_F300wband_r = i_filter
          ELSEIF(i_cons == NWAVE_ALL+40) THEN
             param%iHST_F450wband_r = i_filter
          ENDIF
       ELSEIF(i_cons <= NWAVE_ALL+50) THEN
          IF(i_cons == NWAVE_ALL+41) THEN
             param%iHST_F555wband_r = i_filter
          ELSEIF(i_cons == NWAVE_ALL+42) THEN
             param%iHST_F606wband_r = i_filter
          ELSEIF(i_cons == NWAVE_ALL+43) THEN
             param%iHST_F702wband_r = i_filter
          ELSEIF(i_cons == NWAVE_ALL+44) THEN
             param%iHST_F814wband_r = i_filter
          ELSEIF(i_cons == NWAVE_ALL+45) THEN
             param%iSDSS_upband_r = i_filter
          ELSEIF(i_cons == NWAVE_ALL+46) THEN
             param%iSDSS_gpband_r = i_filter
          ELSEIF(i_cons == NWAVE_ALL+47) THEN
             param%iSDSS_rpband_r = i_filter
          ELSEIF(i_cons == NWAVE_ALL+48) THEN
             param%iSDSS_ipband_r = i_filter
          ELSEIF(i_cons == NWAVE_ALL+49) THEN
             param%iSDSS_zpband_r = i_filter
          ELSEIF(i_cons == NWAVE_ALL+50) THEN
             param%iNLyC_r = i_filter
          ENDIF
       ELSEIF(i_cons <= NWAVE_ALL+60) THEN
          IF(i_cons == NWAVE_ALL+51) THEN
             param%iL1216_r = i_filter
          ELSEIF(i_cons == NWAVE_ALL+52) THEN
             param%iL1400_r = i_filter
          ELSEIF(i_cons == NWAVE_ALL+53) THEN
             param%iL1500_r = i_filter
          ELSEIF(i_cons == NWAVE_ALL+54) THEN
             param%iL1600_r = i_filter
          ELSEIF(i_cons == NWAVE_ALL+55) THEN
             param%iL1700_r = i_filter
          ELSEIF(i_cons == NWAVE_ALL+56) THEN
             param%iL2800_r = i_filter
          ELSEIF(i_cons == NWAVE_ALL+57) THEN
             param%iL4861_r = i_filter
          ELSEIF(i_cons == NWAVE_ALL+58) THEN
             param%iL6563_r = i_filter
          ELSEIF(i_cons == NWAVE_ALL+59) THEN
             param%iIRACch1band_r = i_filter
          ELSEIF(i_cons == NWAVE_ALL+60) THEN
             param%iIRACch2band_r = i_filter
          ENDIF
       ELSEIF(i_cons <= NWAVE_ALL+70) THEN
          IF(i_cons == NWAVE_ALL+61) THEN
             param%iIRACch3band_r = i_filter
          ELSEIF(i_cons == NWAVE_ALL+62) THEN
             param%iIRACch4band_r = i_filter
          ELSEIF(i_cons == NWAVE_ALL+63) THEN
             param%iWFCAM_zband_r = i_filter
          ELSEIF(i_cons == NWAVE_ALL+64) THEN
             param%iWFCAM_Yband_r = i_filter
          ELSEIF(i_cons == NWAVE_ALL+65) THEN
             param%iWFCAM_Jband_r = i_filter
          ELSEIF(i_cons == NWAVE_ALL+66) THEN
             param%iWFCAM_Hband_r = i_filter
          ELSEIF(i_cons == NWAVE_ALL+67) THEN
             param%iWFCAM_Kband_r = i_filter
          ELSEIF(i_cons == NWAVE_ALL+68) THEN
             param%iVIRCAM_Yband_r = i_filter
          ELSEIF(i_cons == NWAVE_ALL+69) THEN
             param%iVIRCAM_Jband_r = i_filter
          ELSEIF(i_cons == NWAVE_ALL+70) THEN
             param%iVIRCAM_Hband_r = i_filter
          ENDIF
       ELSEIF(i_cons <= NWAVE_ALL+80) THEN
          IF(i_cons == NWAVE_ALL+71) THEN
             param%iVIRCAM_Kband_r = i_filter
          ELSEIF(i_cons == NWAVE_ALL+72) THEN
             param%iCIBER_Iband_r = i_filter
          ELSEIF(i_cons == NWAVE_ALL+73) THEN
             param%iCIBER_Hband_r = i_filter
          ELSEIF(i_cons == NWAVE_ALL+74) THEN
             param%iNEWFIRM_J1band_r = i_filter
          ELSEIF(i_cons == NWAVE_ALL+75) THEN
             param%iNEWFIRM_J2band_r = i_filter
          ELSEIF(i_cons == NWAVE_ALL+76) THEN
             param%iNEWFIRM_J3band_r = i_filter
          ELSEIF(i_cons == NWAVE_ALL+77) THEN
             param%iNEWFIRM_H1band_r = i_filter
          ELSEIF(i_cons == NWAVE_ALL+78) THEN
             param%iNEWFIRM_H2band_r = i_filter
          ELSEIF(i_cons == NWAVE_ALL+79) THEN
             param%iNEWFIRM_Ksband_r = i_filter
          ELSEIF(i_cons == NWAVE_ALL+80) THEN
             param%iAKARI_N2band_r = i_filter
          ENDIF
       ELSEIF(i_cons <= NWAVE_ALL+90) THEN
          IF(i_cons == NWAVE_ALL+81) THEN
             param%iAKARI_N3band_r = i_filter
          ELSEIF(i_cons == NWAVE_ALL+82) THEN
             param%iAKARI_N4band_r = i_filter
          ELSEIF(i_cons == NWAVE_ALL+83) THEN
             param%iAKARI_S7band_r = i_filter
          ELSEIF(i_cons == NWAVE_ALL+84) THEN
             param%iAKARI_S9Wband_r = i_filter
          ELSEIF(i_cons == NWAVE_ALL+85) THEN
             param%iAKARI_S11band_r = i_filter
          ELSEIF(i_cons == NWAVE_ALL+86) THEN
             param%iWISH_0band_r = i_filter
          ELSEIF(i_cons == NWAVE_ALL+87) THEN
             param%iWISH_1band_r = i_filter
          ELSEIF(i_cons == NWAVE_ALL+88) THEN
             param%iWISH_2band_r = i_filter
          ELSEIF(i_cons == NWAVE_ALL+89) THEN
             param%iWISH_3band_r = i_filter
          ELSEIF(i_cons == NWAVE_ALL+90) THEN
             param%iWISH_4band_r = i_filter
          ENDIF
       ELSEIF(i_cons <= NWAVE_ALL+100) THEN
          IF(i_cons == NWAVE_ALL+91) THEN
             param%iWISH_5band_r = i_filter
          ELSEIF(i_cons == NWAVE_ALL+92) THEN
             param%iACS_F435Wband_r = i_filter
          ELSEIF(i_cons == NWAVE_ALL+93) THEN
             param%iACS_F475Wband_r = i_filter
          ELSEIF(i_cons == NWAVE_ALL+94) THEN
             param%iACS_F606Wband_r = i_filter
          ELSEIF(i_cons == NWAVE_ALL+95) THEN
             param%iACS_F814Wband_r = i_filter
          ELSEIF(i_cons == NWAVE_ALL+96) THEN
             param%iNICMOS_F160Wband_r = i_filter
          ELSEIF(i_cons == NWAVE_ALL+97) THEN
             param%iWFPC2_F300Wband_r = i_filter
          ELSEIF(i_cons == NWAVE_ALL+98) THEN
             param%iWFPC2_F450Wband_r = i_filter
          ELSEIF(i_cons == NWAVE_ALL+99) THEN
             param%iWFPC2_F555Wband_r = i_filter
          ELSEIF(i_cons == NWAVE_ALL+100) THEN
             param%iWFC3_F125Wband_r = i_filter
          ENDIF
       ELSEIF(i_cons <= NWAVE_ALL+110) THEN
          IF(i_cons == NWAVE_ALL+101) THEN
             param%iWFC3_F140Wband_r = i_filter
          ELSEIF(i_cons == NWAVE_ALL+102) THEN
             param%iWFC3_F160Wband_r = i_filter
          ELSEIF(i_cons == NWAVE_ALL+103) THEN
             param%iHSC_gband_r = i_filter
          ELSEIF(i_cons == NWAVE_ALL+104) THEN
             param%iHSC_rband_r = i_filter
          ELSEIF(i_cons == NWAVE_ALL+105) THEN
             param%iHSC_iband_r = i_filter
          ELSEIF(i_cons == NWAVE_ALL+106) THEN
             param%iHSC_zband_r = i_filter
          ELSEIF(i_cons == NWAVE_ALL+107) THEN
             param%iHSC_Yband_r = i_filter
          ENDIF
       ENDIF
    ENDIF
    ! --- added by MARK (2018/Mar/02)
    IF(param%wdim_band .eqv. .false. .and. &
         ((i_cons >= 50 .and. i_cons <= 58) .or. &
          (i_cons >= NWAVE_ALL+50 .and. i_cons <= NWAVE_ALL+58))) THEN
       ! NLyC, L{1216,1400,1500,1600,1700,2800,4861,6563}
       param%wdim_band = .true.
    ENDIF
  END SUBROUTINE SubstGlobalFilterIntegers
!!$============================================================================
  SUBROUTINE SubstGlobalFilterIntegersForAll
    use global_var

    INTEGER :: i

    DO i = 1, param%nwave
       param%iwave(i) = i
    ENDDO

    ! --- for observer-frame
    param%iBband  =  1; param%iUband  =  2; param%iVband  =  3
    param%iRcband =  4; param%iRJband =  5; param%iIcband =  6
    param%iIJband =  7; param%izband  =  8; param%iJband  =  9
    param%iHband  = 10; param%iKband  = 11; param%iKpband = 12
    param%iLband  = 13
    param%iSuprime_Bband  = 14; param%iSuprime_gband  = 15
    param%iSuprime_Vband  = 16; param%iSuprime_rband  = 17
    param%iSuprime_RLband = 18; param%iSuprime_ipband = 19
    param%iSuprime_Iband  = 20; param%iSuprime_zpband = 21
    param%iSuprime_zp_redband = 22
    param%iobs5um = 23
    param%i2MASS_Jband = 24; param%i2MASS_Hband = 25
    param%i2MASS_Ksband = 26
    param%iACS_F775Wband = 27; param%iACS_F850LPband = 28
    param%iCISCO_zband = 29; param%iCISCO_Jband = 30
    param%iCISCO_Hband = 31; param%iCISCO_Kband = 32
    param%iCISCO_Kpband = 33
    param%iGALEX_FUVband = 34; param%iGALEX_NUVband = 35
    param%iGOODS_RBband = 36; param%iGOODS_RGband = 37
    param%iGOODS_RSband = 38
    param%iHST_F300wband = 39; param%iHST_F450wband = 40
    param%iHST_F555wband = 41; param%iHST_F606wband = 42
    param%iHST_F702wband = 43; param%iHST_F814wband = 44
    param%iSDSS_upband = 45; param%iSDSS_gpband = 46
    param%iSDSS_rpband = 47; param%iSDSS_ipband = 48
    param%iSDSS_zpband = 49
    param%iNLyC  = 50; param%iL1216 = 51; param%iL1400 = 52
    param%iL1500 = 53; param%iL1600 = 54; param%iL1700 = 55
    param%iL2800 = 56; param%iL4861 = 57; param%iL6563 = 58
    param%iIRACch1band = 59; param%iIRACch2band = 60
    param%iIRACch3band = 61; param%iIRACch4band = 62
    param%iWFCAM_zband = 63; param%iWFCAM_Yband = 64
    param%iWFCAM_Jband = 65; param%iWFCAM_Hband = 66
    param%iWFCAM_Kband = 67
    param%iVIRCAM_Yband = 68; param%iVIRCAM_Jband = 69
    param%iVIRCAM_Hband = 70; param%iVIRCAM_Kband = 71
    param%iCIBER_Iband = 72; param%iCIBER_Hband = 73
    param%iNEWFIRM_J1band = 74; param%iNEWFIRM_J2band = 75
    param%iNEWFIRM_J3band = 76; param%iNEWFIRM_H1band = 77
    param%iNEWFIRM_H2band = 78; param%iNEWFIRM_Ksband = 79
    param%iAKARI_N2band = 80; param%iAKARI_N3band = 81
    param%iAKARI_N4band = 82; param%iAKARI_S7band = 83
    param%iAKARI_S9Wband = 84; param%iAKARI_S11band = 85
    param%iWISH_0band = 86; param%iWISH_1band = 87
    param%iWISH_2band = 88; param%iWISH_3band = 89
    param%iWISH_4band = 90; param%iWISH_5band = 91
    param%iACS_F435Wband    = 92; param%iACS_F475Wband   = 93
    param%iACS_F606Wband    = 94; param%iACS_F814Wband   = 95
    param%iNICMOS_F160Wband = 96; param%iWFPC2_F300Wband = 97
    param%iWFPC2_F450Wband  = 98; param%iWFPC2_F555Wband = 99
    param%iWFC3_F125Wband   = 100; param%iWFC3_F140Wband = 101
    param%iWFC3_F160Wband   = 102
    param%iHSC_gband = 103; param%iHSC_rband = 104; param%iHSC_iband = 105
    param%iHSC_zband = 106; param%iHSC_Yband = 107

    ! --- for rest-frame
    param%iBband_r  = NWAVE_ALL+1;  param%iUband_r  = NWAVE_ALL+2
    param%iVband_r  = NWAVE_ALL+3;  param%iRcband_r = NWAVE_ALL+4
    param%iRJband_r = NWAVE_ALL+5;  param%iIcband_r = NWAVE_ALL+6
    param%iIJband_r = NWAVE_ALL+7;  param%izband_r  = NWAVE_ALL+8
    param%iJband_r  = NWAVE_ALL+9;  param%iHband_r  = NWAVE_ALL+10
    param%iKband_r  = NWAVE_ALL+11; param%iKpband_r = NWAVE_ALL+12
    param%iLband_r  = NWAVE_ALL+13
    param%iSuprime_Bband_r  = NWAVE_ALL+14; param%iSuprime_gband_r  = NWAVE_ALL+15
    param%iSuprime_Vband_r  = NWAVE_ALL+16; param%iSuprime_rband_r  = NWAVE_ALL+17
    param%iSuprime_RLband_r = NWAVE_ALL+18; param%iSuprime_ipband_r = NWAVE_ALL+19
    param%iSuprime_Iband_r  = NWAVE_ALL+20; param%iSuprime_zpband_r = NWAVE_ALL+21
    param%iSuprime_zp_redband_r = NWAVE_ALL+22
    param%iobs5um = NWAVE_ALL+23
    param%i2MASS_Jband_r  = NWAVE_ALL+24; param%i2MASS_Hband_r = NWAVE_ALL+25
    param%i2MASS_Ksband_r = NWAVE_ALL+26
    param%iACS_F775Wband_r = NWAVE_ALL+27; param%iACS_F850LPband_r = NWAVE_ALL+28
    param%iCISCO_zband_r  = NWAVE_ALL+29; param%iCISCO_Jband_r = NWAVE_ALL+30
    param%iCISCO_Hband_r  = NWAVE_ALL+31; param%iCISCO_Kband_r = NWAVE_ALL+32
    param%iCISCO_Kpband_r = NWAVE_ALL+33
    param%iGALEX_FUVband_r = NWAVE_ALL+34; param%iGALEX_NUVband_r = NWAVE_ALL+35
    param%iGOODS_RBband_r = NWAVE_ALL+36; param%iGOODS_RGband_r = NWAVE_ALL+37
    param%iGOODS_RSband_r = NWAVE_ALL+38
    param%iHST_F300wband_r = NWAVE_ALL+39; param%iHST_F450wband_r = NWAVE_ALL+40
    param%iHST_F555wband_r = NWAVE_ALL+41; param%iHST_F606wband_r = NWAVE_ALL+42
    param%iHST_F702wband_r = NWAVE_ALL+43; param%iHST_F814wband_r = NWAVE_ALL+44
    param%iSDSS_upband_r = NWAVE_ALL+45; param%iSDSS_gpband_r = NWAVE_ALL+46
    param%iSDSS_rpband_r = NWAVE_ALL+47; param%iSDSS_ipband_r = NWAVE_ALL+48
    param%iSDSS_zpband_r = NWAVE_ALL+49
    param%iNLyC_r  = NWAVE_ALL+50; param%iL1216_r = NWAVE_ALL+51
    param%iL1400_r = NWAVE_ALL+52; param%iL1500_r = NWAVE_ALL+53
    param%iL1600_r = NWAVE_ALL+54; param%iL1700_r = NWAVE_ALL+55
    param%iL2800_r = NWAVE_ALL+56; param%iL4861_r = NWAVE_ALL+57
    param%iL6563_r = NWAVE_ALL+58
    param%iIRACch1band_r = NWAVE_ALL+59; param%iIRACch2band_r = NWAVE_ALL+60
    param%iIRACch3band_r = NWAVE_ALL+61; param%iIRACch4band_r = NWAVE_ALL+62
    param%iWFCAM_zband_r = NWAVE_ALL+63; param%iWFCAM_Yband_r = NWAVE_ALL+64
    param%iWFCAM_Jband_r = NWAVE_ALL+65; param%iWFCAM_Hband_r = NWAVE_ALL+66
    param%iWFCAM_Kband_r = NWAVE_ALL+67
    param%iVIRCAM_Yband_r = NWAVE_ALL+68; param%iVIRCAM_Jband_r = NWAVE_ALL+69
    param%iVIRCAM_Hband_r = NWAVE_ALL+70; param%iVIRCAM_Kband_r = NWAVE_ALL+71
    param%iCIBER_Iband_r = NWAVE_ALL+72; param%iCIBER_Hband_r = NWAVE_ALL+73
    param%iNEWFIRM_J1band_r = NWAVE_ALL+74; param%iNEWFIRM_J2band_r = NWAVE_ALL+75
    param%iNEWFIRM_J3band_r = NWAVE_ALL+76; param%iNEWFIRM_H1band_r = NWAVE_ALL+77
    param%iNEWFIRM_H2band_r = NWAVE_ALL+78; param%iNEWFIRM_Ksband_r = NWAVE_ALL+79
    param%iAKARI_N2band_r = NWAVE_ALL+80; param%iAKARI_N3band_r = NWAVE_ALL+81
    param%iAKARI_N4band_r = NWAVE_ALL+82; param%iAKARI_S7band_r = NWAVE_ALL+83
    param%iAKARI_S9Wband_r = NWAVE_ALL+84; param%iAKARI_S11band_r = NWAVE_ALL+85
    param%iWISH_0band_r = NWAVE_ALL+86; param%iWISH_1band_r = NWAVE_ALL+87
    param%iWISH_2band_r = NWAVE_ALL+88; param%iWISH_3band_r = NWAVE_ALL+89
    param%iWISH_4band_r = NWAVE_ALL+90; param%iWISH_5band_r = NWAVE_ALL+91
    param%iACS_F435Wband_r    = NWAVE_ALL+92; param%iACS_F475Wband_r   = NWAVE_ALL+93
    param%iACS_F606Wband_r    = NWAVE_ALL+94; param%iACS_F814Wband_r   = NWAVE_ALL+95
    param%iNICMOS_F160Wband_r = NWAVE_ALL+96; param%iWFPC2_F300Wband_r = NWAVE_ALL+97
    param%iWFPC2_F450Wband_r  = NWAVE_ALL+98; param%iWFPC2_F555Wband_r = NWAVE_ALL+99
    param%iWFC3_F125Wband_r   = NWAVE_ALL+100; param%iWFC3_F140Wband_r = NWAVE_ALL+101
    param%iWFC3_F160Wband_r   = NWAVE_ALL+102
    param%iHSC_gband_r = NWAVE_ALL+103; param%iHSC_rband_r = NWAVE_ALL+104
    param%iHSC_iband_r = NWAVE_ALL+105; param%iHSC_zband_r = NWAVE_ALL+106
    param%iHSC_Yband_r = NWAVE_ALL+107

    param%wdim_band = .true. ! added by MARK (2018/Mar/02)
  END SUBROUTINE SubstGlobalFilterIntegersForAll
!!$============================================================================
  SUBROUTINE sort(n, arr)
    ! Quicksort Algorithm from Numerical Recipes in Fortran 90
    INTEGER, INTENT(IN) :: n
    DOUBLE PRECISION, INTENT(INOUT) :: arr(n)
    INTEGER, PARAMETER :: M = 7, NSTACK = 50
    INTEGER :: i, ir, j, jstack, k, l, istack(NSTACK)
    DOUBLE PRECISION :: a, temp

    jstack = 0; l = 1; ir = n
1   IF(ir - l < M) THEN
       DO j = l+1, ir
          a = arr(j)
          DO i=j-1, l, -1
             IF(arr(i) <= a) goto 2
             arr(i+1) = arr(i)
          ENDDO
          i = l - 1
2         arr(i+1) = a
       ENDDO
       IF(jstack == 0) return
       ir = istack(jstack); l = istack(jstack-1)
       jstack = jstack - 2
    ELSE
       k = (l + ir) / 2
       temp = arr(k)
       arr(k) = arr(l+1)
       arr(l+1) = temp
       IF(arr(l) > arr(ir)) THEN
          temp = arr(l)
          arr(l) = arr(ir)
          arr(ir) = temp
       ENDIF
       IF(arr(l+1) > arr(ir)) THEN
          temp = arr(l+1)
          arr(l+1) = arr(ir)
          arr(ir) = temp
       ENDIF
       IF(arr(l) > arr(l+1)) THEN
          temp = arr(l)
          arr(l) = arr(l+1)
          arr(l+1) = temp
       ENDIF
       i = l + 1; j = ir; a = arr(l+1)
3      continue
       i = i + 1
       IF(arr(i) < a) goto 3
4      continue
       j = j - 1
       IF(arr(j) > a) goto 4
       IF(j < i) goto 5
       temp = arr(i)
       arr(i) = arr(j)
       arr(j) = temp
       goto 3
5      arr(l+1) = arr(j)
       arr(j) = a
       jstack = jstack + 2
       IF(jstack > NSTACK) THEN
          pause 'NSTACK too small in sort'
       ENDIF
       IF(ir - i + 1 >= j - l) THEN
          istack(jstack) = ir; istack(jstack-1) = i
          ir = j - 1
       ELSE
          istack(jstack) = j - 1; istack(jstack-1) = l
          l = i
       ENDIF
    ENDIF
    goto 1
  END SUBROUTINE sort
!!$============================================================================
  SUBROUTINE GatherResults(nz)
    use MCMCrelated
    implicit none

    INTEGER, INTENT(IN) :: nz
    INTEGER :: bin, i25, i50, i75
    DOUBLE PRECISION :: Square      
    DOUBLE PRECISION, ALLOCATABLE :: arr(:) ! for M_HI/L_B distribution
    INTEGER :: k
    CHARACTER(LEN=3)  :: ci
    INTEGER, ALLOCATABLE :: mlx_size(:), mlx_size_index(:)
    CHARACTER(LEN=50) :: cerr_a = '# GatherResults: fail to allocate'
    CHARACTER(LEN=50) :: cerr_d = '# GatherResults: fail to deallocate'

    !--- merge catalogs
    IF(param%run_type /= 3 .and. inode == 0) THEN
       IF(param%run_type == 1) THEN
          ! galaxy catalog
          catalog_all = trim(fbase)//trim(file_o)//'_1.dat'
          write(file_catalog(len_trim(file_catalog)-5:len_trim(file_catalog)-4), &
               '(A2)') '??'
          command = 'cat '//trim(file_catalog)//' > '//trim(catalog_all)
          call system(trim(command))

          ! quasar catalog
          catalog_all = trim(fbase)//trim(file_o)//'_7.dat'
          write(file_catalog_q(len_trim(file_catalog_q)-5:&
                                        len_trim(file_catalog_q)-4), '(A2)') '??'
          command = 'cat '//trim(file_catalog_q)//' > '//trim(catalog_all)
          call system(trim(command))

          ! target catalog 
          IF(param%traceIDs == 1) THEN
             catalog_all = trim(fbase)//trim(file_o)//'_8.dat'
             write(file_catalog_t(len_trim(file_catalog_t)-5:&
                                           len_trim(file_catalog_t)-4), '(A2)') '??'
             command = 'cat '//trim(file_catalog_t)//' > '//trim(catalog_all)
             call system(trim(command))
          ENDIF
       ELSE IF(param%run_type == 2) THEN
          ! galaxy catalog
          write(ci, '(I3.3)') mrgp%num_step - nz ! = param%izout - 1
          catalog_all = trim(fbase)//trim(file_o)//'_1_'//trim(ci)//'.dat'
          write(file_catalog(len_trim(file_catalog)-9:len_trim(file_catalog)-8), &
               '(A2)') '??'
          write(file_catalog(len_trim(file_catalog)-6:len_trim(file_catalog)-4), &
               '(A3)') trim(ci)
          command = 'cat '//trim(file_catalog)//' > '//trim(catalog_all)
          call system(trim(command))

          ! quasar catalog
          catalog_all = trim(fbase)//trim(file_o)//'_7_'//trim(ci)//'.dat'
          write(file_catalog_q(len_trim(file_catalog_q)-9:&
                                        len_trim(file_catalog_q)-8), '(A2)') '??'
          write(file_catalog_q(len_trim(file_catalog_q)-6:len_trim(file_catalog_q)-4), &
               '(A3)') trim(ci)
          command = 'cat '//trim(file_catalog_q)//' > '//trim(catalog_all)
          call system(trim(command))
       ENDIF
    ENDIF


    ! luminosity function
    lf_n_all(:,:,:) = 0.d0
    call MPI_allreduce(lf(nz)%n, lf_n_all, N1N2Nbin, MPI_DOUBLE, MPI_SUM, &
                       MPI_COMM_WORLD, istat)
    lf(nz)%n = lf_n_all

    lf_n_all(:,:,:) = 0.d0
    call MPI_allreduce(lf_d(nz)%n, lf_n_all, N1N2Nbin, MPI_DOUBLE, MPI_SUM, &
                       MPI_COMM_WORLD, istat)
    lf_d(nz)%n = lf_n_all

    lf_n_all(:,:,:) = 0.d0
    call MPI_allreduce(lf(nz)%n_brst, lf_n_all, N1N2Nbin, MPI_DOUBLE, MPI_SUM, &
                       MPI_COMM_WORLD, istat)
    lf(nz)%n_brst = lf_n_all

    lf_n_all(:,:,:) = 0.d0
    call MPI_allreduce(lf_d(nz)%n_brst, lf_n_all, N1N2Nbin, MPI_DOUBLE, MPI_SUM, &
                       MPI_COMM_WORLD, istat)
    lf_d(nz)%n_brst = lf_n_all

    lf_q_all(:,:,:) = 0.d0
    call MPI_allreduce(lf_q(nz)%n, lf_q_all, param%nwaveq*lf_q(nz)%Nbin, &
                       MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD, istat)
    lf_q(nz)%n = lf_q_all

    IF(nz == 1) THEN ! at z ~ 0
       ! M_BH vs M_bulge relation
       IF(param%Mbh .or. param%run_type == 3) THEN
          MbhMbulge_n_all(:,:,:) = 0.d0
          call MPI_allreduce(MbhMbulge%n, MbhMbulge_n_all, &
                             size(MbhMbulge%n), &
                             MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD, istat)
          MbhMbulge_xn_all(:,:,:) = 0.d0
          call MPI_allreduce(MbhMbulge%xn, MbhMbulge_xn_all, &
                             size(MbhMbulge%xn), &
                             MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD, istat)
          MbhMbulge_xxn_all(:,:,:) = 0.d0
          call MPI_allreduce(MbhMbulge%xxn, MbhMbulge_xxn_all, &
                             size(MbhMbulge%xxn), &
                             MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD, istat)
          MbhMbulge_bin = MbhMbulge%bin
          DO i = 1, MbhMbulge%Nbin
             DO j = 1, 4
                ! --- mean and sigma
                IF(MbhMbulge_n_all(1,j,i) > 0.d0) THEN
                   MbhMbulge_xn_all(1,j,i)  = MbhMbulge_xn_all(1,j,i) &
                                              / MbhMbulge_n_all(1,j,i)
                   MbhMbulge_xxn_all(1,j,i) = MbhMbulge_xxn_all(1,j,i) &
                                              / MbhMbulge_n_all(1,j,i)
                   IF(MbhMbulge_xxn_all(1,j,i) - Square(MbhMbulge_xn_all(1,j,i)) &
                        < 0.d0) THEN
                      MbhMbulge_xxn_all(1, j, i) = 0.d0
                   ELSE
                      MbhMbulge_xxn_all(1,j,i) &
                           = sqrt(MbhMbulge_xxn_all(1,j,i) &
                                  - Square(MbhMbulge_xn_all(1,j,i)))
                   ENDIF
                ELSE
                   MbhMbulge_xn_all(1,j,i) = 0.d0; MbhMbulge_xxn_all(1,j,i) = 0.d0
                ENDIF
             ENDDO
          ENDDO
       ENDIF

       IF(param%run_type == 3) THEN
          DiskScale_n_all(:,:,:) = 0.d0
          call MPI_allreduce(DiskScale%n, DiskScale_n_all, &
                             size(DiskScale%n), &
                             MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD, istat)
          DiskScale_xn_all(:,:,:) = 0.d0
          call MPI_allreduce(DiskScale%xn, DiskScale_xn_all, &
                             size(DiskScale%xn), &
                             MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD, istat)
          DiskScale_xxn_all(:,:,:) = 0.d0
          call MPI_allreduce(DiskScale%xxn, DiskScale_xxn_all, &
                             size(DiskScale%xxn), &
                             MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD, istat)
          DiskScale_bin = DiskScale%bin
          DO i = 1, DiskScale%Nbin
            DO j = 1, 4
                ! --- mean and sigma
                IF(DiskScale_n_all(1,j,i) > 0.d0) THEN
                   DiskScale_xn_all(1,j,i)  = DiskScale_xn_all(1,j,i) &
                                              / DiskScale_n_all(1,j,i)
                   DiskScale_xxn_all(1,j,i) = DiskScale_xxn_all(1,j,i) &
                                              / DiskScale_n_all(1,j,i)
                   IF(DiskScale_xxn_all(1,j,i) - Square(DiskScale_xn_all(1,j,i)) &
                        < 0.d0) THEN
                      DiskScale_xxn_all(1, j, i) = 0.d0
                   ELSE
                      DiskScale_xxn_all(1,j,i) &
                           = sqrt(DiskScale_xxn_all(1,j,i) &
                                  - Square(DiskScale_xn_all(1,j,i)))
                   ENDIF
                ELSE
                   DiskScale_xn_all(1,j,i) = 0.d0; DiskScale_xxn_all(1,j,i) = 0.d0
                ENDIF
             ENDDO
          ENDDO
    
          SphScale_n_all(:,:,:) = 0.d0
          call MPI_allreduce(SphScale%n, SphScale_n_all, &
                             size(SphScale%n), &
                             MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD, istat)
          SphScale_xn_all(:,:,:) = 0.d0
          call MPI_allreduce(SphScale%xn, SphScale_xn_all, &
                             size(SphScale%xn), &
                             MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD, istat)
          SphScale_xxn_all(:,:,:) = 0.d0
          call MPI_allreduce(SphScale%xxn, SphScale_xxn_all, &
                             size(SphScale%xxn), &
                             MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD, istat)
          SphScale_bin = SphScale%bin
          DO i = 1, SphScale%Nbin
             Do j = 1, 4
                ! --- mean and sigma
                IF(SphScale_n_all(1,j,i) > 0.d0) THEN
                   SphScale_xn_all(1,j,i)  = SphScale_xn_all(1,j,i) &
                                              / SphScale_n_all(1,j,i)
                   SphScale_xxn_all(1,j,i) = SphScale_xxn_all(1,j,i) &
                                              / SphScale_n_all(1,j,i)
                   IF(SphScale_xxn_all(1,j,i) - Square(SphScale_xn_all(1,j,i)) &
                        < 0.d0) THEN
                      SphScale_xxn_all(1, j, i) = 0.d0
                   ELSE
                      SphScale_xxn_all(1,j,i) &
                           = sqrt(SphScale_xxn_all(1,j,i) &
                                  - Square(SphScale_xn_all(1,j,i)))
                   ENDIF
                ELSE
                   SphScale_xn_all(1,j,i) = 0.d0; SphScale_xxn_all(1,j,i) = 0.d0
                ENDIF
             ENDDO
          ENDDO
       ENDIF
    ENDIF ! at z ~ 0

    ! M(H)/L(B)
    IF(param%ML) THEN
       ML_x_all(:,:,:) = 0.d0
       allocate(mlx_size(nnode), stat=ier); call CheckIerr(ier, &
            trim(cerr_a)//' mlx_size')
       allocate(mlx_size_index(nnode+1), stat=ier); call CheckIerr(ier, &
            trim(cerr_a)//' mlx_size_index')
       ! for median
       DO j = 1, N2
          DO i = 1, ml(nz)%Nbin
             mlx_size = 0; mlx_size_index = 0
             ! get the number of galaxies in each bin, for each node
             call MPI_allGather(int(ml(nz)%n(2,j,i)), 1, MPI_INTEGER, &
                                mlx_size,   1, MPI_INTEGER, &
                                MPI_COMM_WORLD, istat)
             mlx_size_index(1) = 0
             DO k = 1, nnode
                mlx_size_index(k+1) = mlx_size_index(k) + mlx_size(k)
             ENDDO
             ! gather the ml(nz)%x of each node into ML_x_all
             call MPI_GatherV(ml(nz)%x(:,j,i), mlx_size(inode+1), &
                              MPI_DOUBLE_PRECISION, &
                              ML_x_all(:,j,i), mlx_size, mlx_size_index, &
                              MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, istat)

          ENDDO
       ENDDO
       deallocate(mlx_size, stat=ier); call CheckIerr(ier, &
            trim(cerr_d)//' mlx_size')
       deallocate(mlx_size_index, stat=ier); call CheckIerr(ier, &
            trim(cerr_d)//' mlx_size_index')
       IF(inode == 0) THEN
         deallocate(ml(nz)%x, stat=ier); call CheckIerr(ier, &
              trim(cerr_d)//' ml(nz)%x')
         allocate(ml(nz)%x(NgalMax*nnode,N2,ml(nz)%Nbin), stat=ier); call &
              CheckIerr(ier, trim(cerr_a)//' ml(nz)%x')
         ml(nz)%x = ML_x_all
       ENDIF

       ML_n_all(:,:,:) = 0.d0
       call MPI_allreduce(ml(nz)%n, ML_n_all, size(ml(nz)%n), MPI_DOUBLE, &
                          MPI_SUM, MPI_COMM_WORLD, istat)
       ml(nz)%n = ML_n_all

       ML_n_all(:,:,:) = 0.d0
       call MPI_allreduce(ml(nz)%xn, ML_n_all, size(ml(nz)%xn), MPI_DOUBLE, &
                          MPI_SUM, MPI_COMM_WORLD, istat)
       ml(nz)%xn = ML_n_all

       ML_n_all(:,:,:) = 0.d0
       call MPI_allreduce(ml(nz)%xxn, ML_n_all, size(ml(nz)%xxn), MPI_DOUBLE, &
                          MPI_SUM, MPI_COMM_WORLD, istat)
       ml(nz)%xxn = ML_n_all

    ENDIF

    ! mass functions
    mf_n_all(:,:,:) = 0.d0
    call MPI_allreduce(mf(nz)%n, mf_n_all, size(mf(nz)%n), MPI_DOUBLE, &
                       MPI_SUM, MPI_COMM_WORLD, istat)
    mf(nz)%n = mf_n_all

    mf_n_all(:,:,:) = 0.d0
    call MPI_allreduce(mf(nz)%n_brst, mf_n_all, size(mf(nz)%n), MPI_DOUBLE, &
                       MPI_SUM, MPI_COMM_WORLD, istat)
    mf(nz)%n_brst = mf_n_all

    IF(nz == 1) THEN
       ! CSFR
       CSFR_all(2,:) = 0.d0
       call MPI_allreduce(CSFR(2,:), CSFR_all(2,:), size(CSFR(2,:)), MPI_DOUBLE, &
                          MPI_SUM, MPI_COMM_WORLD, istat)
       CSFR = CSFR_all

       SMF_z_all(:,:) = 0.d0
       call MPI_allreduce(SMF_z, SMF_z_all, size(SMF_z_all), MPI_DOUBLE, &
                          MPI_SUM, MPI_COMM_WORLD, istat)
       SMF_z = SMF_z_all
    ENDIF
  END SUBROUTINE GatherResults
!!$===========================================================================
!  SUBROUTINE WriteTmpResults(nz)
!    use MCMCrelated
!    INTEGER, INTENT(IN) :: nz
!    INTEGER :: ibin, iband, imor, ier
!    INTEGER :: iband0(3) = (/26, 47, 34/) ! K-, r-, FUV-bands
!    CHARACTER(LEN=50) :: cerr = '# WriteTmpResults: fail to open file '
!    CHARACTER(LEN=100) :: fname
!    DOUBLE PRECISION :: CSFR0
!
!    iout = 97
!    IF(inode == 0) THEN
!       fname = 'LF_'//trim(file_o)//'.dat'
!       open(iout, file = trim(fname), status = 'replace', iostat=ier); call &
!            CheckIerr(ier, trim(cerr)//' '//trim(fname))
!       DO ibin = 1, 100
!          write(iout, '(F11.4, $)') lf(nz)%bin(ibin)
!          DO iband = 1, 3
!             write(iout, '(3F11.4, $)') &
!                  (lf_mcmc(iband0(iband), imor, ibin) * lf(nz)%invstep, &
!                  imor=1,3)
!          ENDDO
!          write(iout, *)
!       ENDDO
!       close(iout)
!
!       fname = 'MbhMbulge_'//trim(file_o)//'.dat'
!       open(iout, file = trim(fname), status = 'replace', iostat=ier); call &
!            CheckIerr(ier, trim(cerr)//' '//trim(fname))
!       DO ibin = 1, 50
!          write(iout, *) MbhMbulge_bin(ibin), &
!               MbhMbulge_xn_all(1,4,ibin), &
!               MbhMbulge_xxn_all(1,4,ibin)
!       ENDDO
!       close(iout)
!
!       fname = 'DiskScale_'//trim(file_o)//'.dat'
!       open(iout, file = trim(fname), status = 'replace', iostat=ier); call &
!            CheckIerr(ier, trim(cerr)//' '//trim(fname))
!       DO ibin = 1, 50
!          write(iout, *) DiskScale_bin(ibin), &
!               DiskScale_xn_all(1,3,ibin), &
!               DiskScale_xxn_all(1,3,ibin)
!       ENDDO
!       close(iout)
!
!       fname = 'SphScale_'//trim(file_o)//'.dat'
!       open(iout, file = trim(fname), status = 'replace', iostat=ier); call &
!            CheckIerr(ier, trim(cerr)//' '//trim(fname))
!       DO ibin = 1, 50
!          write(iout, *) SphScale_bin(ibin), &
!               SphScale_xn_all(1,1,ibin), &
!               SphScale_xxn_all(1,1,ibin)
!       ENDDO
!       close(iout)
!
!       fname = 'MHI_'//trim(file_o)//'.dat'
!       open(iout, file = trim(fname), status = 'replace', iostat=ier); call &
!            CheckIerr(ier, trim(cerr)//' '//trim(fname))
!       DO ibin = 1, mf(nz)%Nbin
!          write(iout, *) mf(nz)%bin(ibin), mf(nz)%invstep &
!               * (mf_mcmc(10,1,ibin) + mf_mcmc(10,2,ibin) + mf_mcmc(10,3,ibin))
!       ENDDO
!       close(iout)
!
!       fname = 'Ms_'//trim(file_o)//'.dat'
!       open(iout, file = trim(fname), status = 'replace', iostat=ier); call &
!            CheckIerr(ier, trim(cerr)//' '//trim(fname))
!       DO ibin = 1, mf(nz)%Nbin
!          write(iout, *) mf(nz)%bin(ibin), &
!               (mf_mcmc(9,1,ibin) + mf_mcmc(9,2,ibin)) * mf(nz)%invstep, &
!               mf_mcmc(9,3,ibin) * mf(nz)%invstep
!       ENDDO
!       close(iout)
!
!       fname = 'CSFR_'//trim(file_o)//'.dat'
!       open(iout, file = trim(fname), status = 'replace', iostat=ier); call &
!            CheckIerr(ier, trim(cerr)//' '//trim(fname))
!       CSFR0 = 1.65d0 / 1.d+6 * param%munit / param%th_yr
!       DO ibin = 1, mrgp%num_step ! CSFR_stepnum
!          write(iout, '(2F11.4)') CSFR_all(1,ibin), CSFR0 * CSFR_all(2,ibin)
!       ENDDO
!       close(iout)
!    ENDIF
!
!    fname = 'SMBH_'//trim(file_o)//'.dat'
!    open(iout, file = trim(fname), status = 'replace', iostat=ier); call &
!         CheckIerr(ier, trim(cerr)//' '//trim(fname))
!    DO ibin = 1, mf(nz)%Nbin
!       write(iout, *) mf(nz)%bin(ibin), &
!            (mf_mcmc(4,1,ibin) + mf_mcmc(4,2,ibin)) * mf(nz)%invstep, &
!            mf_mcmc(4,3,ibin) * mf(nz)%invstep
!    ENDDO
!    close(iout)
!
!    fname = 'SMF_z2p0_'//trim(file_o)//'.dat'
!    open(iout, file = trim(fname), status = 'replace', iostat=ier); call &
!         CheckIerr(ier, trim(cerr)//' '//trim(fname))
!    DO ibin = 1, mf(nz)%Nbin
!       write(iout, *) mf(nz)%bin(ibin), &
!            SMF_z_all(24,ibin) * inv_V * mf(nz)%invstep
!    ENDDO
!    close(iout)
!
!  END SUBROUTINE WriteTmpResults
!!$============================================================================
!  SUBROUTINE WriteMeanResults(nz, iterate_num)
!     use MCMCrelated
!     INTEGER, INTENT(IN) :: nz
!     INTEGER, INTENT(IN) :: iterate_num
!     DOUBLE PRECISION :: mean(3), var(3)
!
!     IF(inode == 0) THEN
!        open(97, file='LF_'//trim(file_o)//'.dat', status = 'replace')
!        DO j = 1, 100
!           mean(1) = lf_n_mean(3,j)/iterate_num*lf(nz)%invstep
!           mean(2) = lf_n_mean(4,j)/iterate_num*lf(nz)%invstep
!           var(1)  = abs(lf_n_var(3,j)/iterate_num*lf(nz)%invstep**2-mean(1)**2)
!           var(2)  = abs(lf_n_var(4,j)/iterate_num*lf(nz)%invstep**2-mean(2)**2)
!           write(97, '(7E15.4e3)') lf(nz)%bin(j), &
!                mean(1), sqrt(var(1)), mean(2), sqrt(var(2)) ! K, r
!        ENDDO
!        close(97)
!
!        open(97, file = 'MHI_'//trim(file_o)//'.dat', status = 'replace')
!        DO j = 1, mf(nz)%Nbin
!           mean(3) = mf_n_mean(10,j)/iterate_num*mf(nz)%invstep
!           var(3) = abs(mf_n_var(10,j)/iterate_num*mf(nz)%invstep**2-mean(3)**2)
!           write(97, '(7E15.4e3)') mf(nz)%bin(j), mean(3), sqrt(var(3)) ! HI MF
!        ENDDO
!        close(97)
!     ENDIF
!  END SUBROUTINE
!!$============================================================================
  SUBROUTINE ReadPetSBcor
    INTEGER, PARAMETER :: N_HEADER = 2
    INTEGER :: ier, irow, irow_mh, iRratio, iBT, iread = 1, i
    CHARACTER(LEN=100) :: fname = 'inpdata/Pet2TotRatio.dat'
    CHARACTER(LEN=500) :: buf
    CHARACTER(LEN=50)  :: cerr_o = '# ReadPetSBcor: fail to open file'
    CHARACTER(LEN=50)  :: cerr_a = '# ReadPetSBcor: fail to allocate'
    DOUBLE PRECISION   :: tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8


    ! SB%mu0 = 37.815d0 + 2.5d0 * log10(2.d0 / const%PI) + 5.d0 * log10(param%h)
    SB%mu0 = (1.d+3 / param%h) / 10.d0 & ! converting [kpc/h] to [rad] at 10 pc
             * (180.d0 / const%PI) * 60.d0 * 60.d0 ! converting [rad] to [arcsec]
    SB%mu0 = 2.5d0 * log10(const%PI) + 5.d0 * log10(SB%mu0)
             ! = 2.5*log10[surface area of r = 1 kpc/h at 10 pc in arcsec^2]
    SB%mu0 = SB%mu0 + 5.d0 * log10(param%h) ! converting m-5logh to m
    SB%N_Rratiop1 = SB%N_Rratio + 1; SB%N_BTp1 = SB%N_BT + 1
    allocate(SB%Pet2TotRatio(3, SB%N_Rratiop1, SB%N_BTp1), stat=ier); call &
         CheckIerr(ier, trim(cerr_a)//' SB%Pet2TotRatio')
    SB%pet2TotRatio(:,:,:) = 0.d0
    open(iread, file = trim(fname), status = 'old', iostat = ier); call &
         CheckIerr(ier, trim(cerr_o)//' '//trim(fname))
    ier = 0; irow = 1
    DO WHILE(ier == 0)
       read(iread, '(A)', iostat=ier) buf
       IF(ier == 0) THEN
          IF(irow > N_HEADER) THEN
             irow_mh = irow - N_HEADER
             iRratio = irow_mh / SB%N_BTp1 + 1; iBT = mod(irow_mh, SB%N_BTp1)
             IF(iBT == 0) THEN
                iBT = SB%N_BTp1; iRratio = iRratio - 1
             ENDIF
             read(buf, *) tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, &
                  SB%Pet2TotRatio(1, iRratio, iBT), & ! -2.5log10(L_pet/L_tot)
                  tmp7, tmp8, &
                  SB%Pet2TotRatio(2, iRratio, iBT), & ! R50(L_pet)/Re(bulge)
                  SB%Pet2TotRatio(3, iRratio, iBT)    ! R50(L_pet)/Re(disk)
!!$               print '(3I10, 3G10.3)', irow_mh, iRratio, iBT, &
!!$                    (SB%Pet2TotRatio(i, iRratio, iBT), i=1,3)
          ENDIF
          irow = irow + 1
       ENDIF
    ENDDO
    close(iread)
  END SUBROUTINE ReadPetSBcor
!!$============================================================================
  SUBROUTINE ReadCST
    INTEGER, PARAMETER :: N_HEADER = 2
    INTEGER :: ier, irow, iz, iM, iread = 1
    CHARACTER(LEN=100) :: fname = 'inpdata/cst_Prada12.dat'
    CHARACTER(LEN=500) :: buf
    CHARACTER(LEN=50)  :: cerr_a = '# ReadCST: fail to allocate'
    CHARACTER(LEN=50)  :: cerr_o = '# ReadCST: fail to open file'
    DOUBLE PRECISION   :: tmp1, tmp2
    INTEGER :: Nz, NM


    Nz = 50
    NM = 100
    allocate(CP%CST(Nz, NM), stat=ier); call CheckIerr(ier, trim(cerr_a)//' CST')
    CP%CST(:,:) = 0.d0
    open(iread, file = trim(fname), status = 'old', iostat = ier); call &
         CheckIerr(ier, trim(cerr_o)//' '//trim(fname))
    ier = 0; irow = 1
    DO WHILE(ier == 0)
       read(iread, '(A)', iostat=ier) buf
       IF(ier == 0) THEN
          IF(irow > N_HEADER) THEN
             iz = (irow-N_HEADER)/NM+1
             iM = mod(irow-N_HEADER, NM)
             IF(iM == 0) THEN
                iM = NM; iz = iz -1
             ENDIF
             read(buf, *) tmp1, tmp2, CP%CST(iz, iM)
          ENDIF
          irow = irow + 1
       ENDIF
    ENDDO
    close(iread)
  END SUBROUTINE ReadCST
!!$============================================================================
END PROGRAM main_nugc
!!$============================================================================
DOUBLE PRECISION FUNCTION CalBulgeReff(rb) RESULT(rb0)
  use global_var
  implicit none
  DOUBLE PRECISION, INTENT(IN) :: rb

  rb0 = param%fb * 0.744d0 * rb !kpc/h
END FUNCTION CalBulgeReff
!!$============================================================================
DOUBLE PRECISION FUNCTION CalDiskReff(rd, rb0) RESULT(rd0)
  use global_var
  implicit none
  DOUBLE PRECISION, INTENT(IN) :: rd, rb0
  DOUBLE PRECISION :: lamH
  REAL :: gasdev ! function

  ! Mo et al. 1998
  !lamH = 0.05d0 * exp(dble(gasdev(param%idum)) * 0.5d0) ! spin parameter

  ! Bett et al. 2007
  lamH = 0.042d0 * exp(dble(gasdev(param%idum)) * 0.26d0) ! spin parameter

!!$       rd0 = fd * lamH * gal(i)%rdisk / sqrt(2.d0) / 1.68d0
  rd0 = lamH * rd * 1.68d0 / sqrt(2.d0) ! [kpc/h]
        ! Fall79; Fall & Efstathiou80; Fall83
  rd0 = rd0 * param%zsp1**param%alpha_rad
  IF(rd0 <= 0.d0) rd0 = rb0
END FUNCTION CalDiskReff
!!$============================================================================
INTEGER FUNCTION DetMorType(Lbulge, Ldisk) RESULT(mor)
  implicit none
  DOUBLE PRECISION, INTENT(IN) :: Lbulge, Ldisk
  DOUBLE PRECISION, PARAMETER :: EPS = 0.d0
  DOUBLE PRECISION :: bt

  IF(Lbulge + Ldisk <= EPS) THEN
     print '(A, G8.1, A)', '(DetMorType) Lbulge + Ldisk < EPS = ', EPS, '!!'
     print '(2(A, G10.3))',  '             --- Lbulge = ', Lbulge, ', Ldisk = ', Ldisk
     bt = 0; mor = 0 ! undefined
  ELSE
     bt = Lbulge / (Lbulge + Ldisk); mor = 3 ! S
  ENDIF
  IF(bt > 0.6d0) THEN
     mor = 1 ! E
  ELSEIF(bt > 0.4d0) THEN
     mor = 2 ! S0
  ENDIF
END FUNCTION DetMorType
!!$============================================================================
SUBROUTINE CheckIONum(ionum, charac)
  INTEGER, INTENT(IN) :: ionum
  CHARACTER(LEN=*), INTENT(IN) :: charac
  INTEGER, PARAMETER :: Nsave_max = 10000
  INTEGER :: i, icallp1
  INTEGER, SAVE :: icall = 0, ionum_used(Nsave_max)
  CHARACTER(LEN=500), SAVE :: c_ionum_used(Nsave_max)

  icallp1 = icall + 1
  IF(icall == 0) THEN
     ionum_used(icallp1) = ionum; c_ionum_used(icallp1) = trim(charac)
  ELSE
     DO i = 1, icall
        IF(ionum == ionum_used(i)) THEN
           print '(A, I5, A)', '# The requested io number of ', ionum, &
               ' in '//trim(charac)//'is already used in '// &
               trim(c_ionum_used(i))//' !!'; exit
        ENDIF
     ENDDO
     ionum_used(icallp1) = ionum; c_ionum_used(icallp1) = trim(charac)
  ENDIF
!!$  --- to check the used I/O number w/ the subroutine name in which
!!$        the I/O number is used
!!$  print '(A, X, 2I5, X, A)', '(CheckIONum)', icall, ionum_used(icallp1), &
!!$       trim(c_ionum_used(icallp1))
  icall = icallp1
END SUBROUTINE CheckIONum
!!$============================================================================
SUBROUTINE CheckIerr(ier, cerr)
  implicit none
  INTEGER, INTENT(IN) :: ier
  CHARACTER(LEN=*), INTENT(IN) :: cerr

  IF(ier /= 0) THEN
     print '(A, I4)', trim(cerr)//', stat = ', ier; stop
  ENDIF
END SUBROUTINE CheckIerr
!!$============================================================================
DOUBLE PRECISION FUNCTION AbsDiff(xbase, x) RESULT(diff)
  implicit none
  DOUBLE PRECISION, INTENT(IN) :: xbase, x

  diff = x
  IF(xbase /= 0.d0) diff = 1.d0 - x / xbase
  IF(diff < 0.d0) diff = -diff
END FUNCTION AbsDiff
!!$============================================================================
DOUBLE PRECISION FUNCTION ObsFracAGN_UV(Lbol, zp1)
   use global_var
   DOUBLE PRECISION, INTENT(IN) :: Lbol, zp1
   DOUBLE PRECISION :: coeff_hs(4) = (/0.17d0, -0.085d0, 0.073d0, -0.085d0/) 
   DOUBLE PRECISION :: coeff_ph(2) = (/1.243d0, 0.066d0/) 
   DOUBLE PRECISION :: a,b
   IF(param%ObsFrac == 1) THEN ! Shirakata et al. (2018)
      a = coeff_hs(1) * zp1 ** coeff_hs(2)
      b = coeff_hs(3) * zp1 ** coeff_hs(4)
   ELSE IF(param%ObsFrac == 2) THEN ! Hopkins et al. (2007)
      a = coeff_ph(1)
      b = coeff_ph(2)
   ENDIF
   ObsFracAGN_UV = a * (Lbol*1.d-46) ** b
END FUNCTION ObsFracAGN_UV
!!$ ========================================================================
