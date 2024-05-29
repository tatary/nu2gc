! ========================================================================
SUBROUTINE SetAndOpenOutputFilesForLAE(fbase, file_o, run_redshift, nnode, inode)
  use LAErelated
  implicit none

  CHARACTER(LEN=*), INTENT(IN) :: fbase, file_o
  INTEGER, INTENT(IN) :: run_redshift ! 1:run at a redshift, 2:run at all redshift
  INTEGER, INTENT(IN) :: nnode, inode ! for parallel
  INTEGER :: i, ii, j, ier
  CHARACTER(LEN=100) :: fbaseLAE
  CHARACTER(LEN=50) :: cerr = '# (LAE)SetAndOpenOutputFilesForLAE: fail to '
  CHARACTER(LEN=2)  :: ci


  allocate(fnameLAE(paramLAE%nfile), stat=ier)
  IF(ier /= 0) THEN
     print *, trim(cerr)//'allocate stat=', ier; stop
  ENDIF

  IF(inode == 0) print '(A)', '# (LAE)Output File Names:'
  fbaseLAE = trim(fbase)//trim(file_o)//'_'
  DO i = 1, paramLAE%nfile
     ii = i + paramLAE%base
     call CheckIONum(ii, 'SetAndOpenOutputFilesForLAE')
     IF(i == 1) THEN
        fnameLAE(i) = trim(fbaseLAE)//'intLyaLF.dat';  iLAE%LyaLF = ii
     ELSEIF(i == 2) THEN
        fnameLAE(i) = trim(fbaseLAE)//'contLFs.dat';   iLAE%contLF = ii
     ELSEIF(i == 3) THEN
        fnameLAE(i) = trim(fbaseLAE)//'NgasDist.dat';  iLAE%dist(1) = ii
     ELSEIF(i == 4) THEN
        fnameLAE(i) = trim(fbaseLAE)//'ZcDist.dat';    iLAE%dist(2) = ii
     ELSEIF(i == 5) THEN
        fnameLAE(i) = trim(fbaseLAE)//'NcZcDist.dat';  iLAE%dist(3) = ii
     ELSEIF(i == 6) THEN
        fnameLAE(i) = trim(fbaseLAE)//'MstarDist.dat'; iLAE%dist(4) = ii
     ELSEIF(i == 7) THEN
        fnameLAE(i) = trim(fbaseLAE)//'MhostDist.dat'; iLAE%dist(5) = ii
     ELSEIF(i == 8) THEN
        fnameLAE(i) = trim(fbaseLAE)//'ReffDist.dat';  iLAE%dist(6) = ii
     ELSEIF(i == 9) THEN
        fnameLAE(i) = trim(fbaseLAE)//'SFRDist.dat';   iLAE%dist(7) = ii
     ELSEIF(i == 10) THEN
        fnameLAE(i) = trim(fbaseLAE)//'TeffDist.dat';  iLAE%dist(8) = ii
     ELSEIF(i == 11) THEN
        fnameLAE(i) = trim(fbaseLAE)//'LyaNgasDist.dat'; iLAE%distLya(1) = ii
     ELSEIF(i == 12) THEN
        fnameLAE(i) = trim(fbaseLAE)//'LyaZcDist.dat'; iLAE%distLya(2) = ii
     ELSEIF(i == 13) THEN
        fnameLAE(i) = trim(fbaseLAE)//'LyaMstarDist.dat'; iLAE%distLya(3) = ii
     ELSEIF(i == 14) THEN
        fnameLAE(i) = trim(fbaseLAE)//'LyaNcZcDist.dat'; iLAE%distLya(4) = ii
     ELSEIF(i == 15) THEN
        fnameLAE(i) = trim(fbaseLAE)//'LyAbright.dat'; iLAE%eachLAE = ii
     ELSEIF(i == 16) THEN
        fnameLAE(i) = trim(fbaseLAE)//'UVbright.dat'; iLAE%each = ii
     ENDIF
!!$     j = index(fnameLAE(i), ' ') - 1
!!$     fnameLAE(i) = fnameLAE(i)(1:j)

     IF(run_redshift == 1) THEN ! run at a redshift
        open(ii, file=trim(fnameLAE(i)), status='unknown', iostat=ier)
        call CheckIerr(ier, trim(cerr)//'open file: '//trim(fnameLAE(i)))
     ELSEIF(run_redshift == 2) THEN ! run at all redshift
        j = index(fnameLAE(i), '.dat') - 1
        fnameLAE(i) = fnameLAE(i)(1:j)//'_***.dat'
     ENDIF

     IF(inode == 0) print '(A)', '  --- '//trim(fnameLAE(i))
  ENDDO

  ! for parallel run (added by MARK on 2017/Mar/16)
  IF(run_redshift == 1) THEN
     write(ci, '(I2.2)') inode+1
     fnameLAE_catalog_LyA = trim(fbaseLAE)//'LyAbright_catalog_'//ci//'.dat'
     fnameLAE_catalog_UV  = trim(fbaseLAE)//'UVbright_catalog_'//ci//'.dat'

     i = paramLAE%base + paramLAE%nfile + inode + 1
     call CheckIONum(i, 'SetAndOpenOutputFilesForLAE')
     open(i, file=trim(fnameLAE_catalog_LyA), status='unknown', iostat=ier)
     call CheckIerr(ier, trim(cerr)//'open file: '//trim(fnameLAE_catalog_LyA))

     i = paramLAE%base + paramLAE%nfile + nnode + inode + 1
     call CheckIONum(i, 'SetAndOpenOutputFilesForLAE')
     open(i, file=trim(fnameLAE_catalog_UV), status='unknown', iostat=ier)
     call CheckIerr(ier, trim(cerr)//'open file: '//trim(fnameLAE_catalog_UV))
  ENDIF

  ! --- for Makiya related (2011/01/31)
  ! iLAE%Makiya = paramLAE%base + paramLAE%nfile + 1
  ! call SetAndOpenOutputFileForMakiya(fbaseLAE, fend)
END SUBROUTINE SetAndOpenOutputFilesForLAE
! ========================================================================
SUBROUTINE AllocLAE(Ngal, nwave)
  use LAErelated
  implicit none
  INTEGER, INTENT(IN) :: Ngal, nwave
  INTEGER :: ier, tnwave
  CHARACTER(LEN=50) :: cerr = '# (LAE)AllocLAE: fail to allocate'

  allocate(allgal(Ngal),      stat=ier); call CheckIerr(ier, trim(cerr)//' allgal')
  allocate(allgal_prev(Ngal), stat=ier); call CheckIerr(ier, trim(cerr)//' allgal_prev')
  allocate(allgal_next(Ngal), stat=ier); call CheckIerr(ier, trim(cerr)//' allgal_next')

  tnwave = 2 * nwave
  call InitializeTypeAllGalaxiesForLAE(allgal,      Ngal, tnwave)
  call InitializeTypeAllGalaxiesForLAE(allgal_prev, Ngal, tnwave)
  call InitializeTypeAllGalaxiesForLAE(allgal_next, Ngal, tnwave)
CONTAINS
  ! ========================================================================
  SUBROUTINE InitializeTypeAllGalaxiesForLAE(allgal, numg, tnwave)
!    use LAErelated
    implicit none
    INTEGER, INTENT(IN) :: numg, tnwave
    TYPE(PhysicalQuantitiesForAllGalaxies), INTENT(INOUT) :: allgal(numg)
    INTEGER :: i
    CHARACTER(LEN=100) :: cerr = '# (LAE)InitializeTypeAllGalaxiesForLAE: fail to allocate'

    allgal(:)%num_b  = 0
    allgal(:)%Zc     = 0.d0; allgal(:)%Zc_pre = 0.d0; allgal(:)%Mc_pre = 0.d0
    allgal(:)%Ms_pre = 0.d0; allgal(:)%rb_pre = 0.d0; allgal(:)%Vb_pre = 0.d0
    DO i = 1, numg
       allgal(i)%SFR(:) = 0.d0; allgal(i)%beta(:) = 0.d0; allgal(i)%taueff(:,:) = 0.d0
       allocate(allgal(i)%lumg_pre(tnwave), stat=ier); call &
            CheckIerr(ier, trim(cerr)//' allgal(i)%lumg_pre')
       allgal(i)%lumg_pre(:) = 0.d0
    ENDDO
  END SUBROUTINE InitializeTypeAllGalaxiesForLAE
END SUBROUTINE AllocLAE
! ========================================================================
SUBROUTINE DeallocLAE
  use LAErelated
  implicit none
  INTEGER :: ier
  CHARACTER(LEN=50) :: cerr = '# (LAE)DeallocLAE: fail to deallocate'

  deallocate(allgal,      stat=ier); call CheckIerr(ier, trim(cerr)//' allgal')
  deallocate(allgal_prev, stat=ier); call CheckIerr(ier, trim(cerr)//' allgal_prev')
  deallocate(allgal_next, stat=ier); call CheckIerr(ier, trim(cerr)//' allgal_next')
END SUBROUTINE DeallocLAE
! ========================================================================
SUBROUTINE AllocateDistributionFunctions_LAE(nz_end, nwave)
  use LAErelated
  implicit none
  INTEGER, INTENT(IN) :: nz_end, nwave
  INTEGER :: i, ier, N1, N2, N3, Nbin
  CHARACTER(LEN=60) :: cerr = '# (LAE)AllocateDistributionFunctions_LAE: fail to allocate'

  ! --- Luminosity Functions
  allocate(LF_LAE(nz_end),  stat=ier); call CheckIerr(ier, trim(cerr)//' LF_LAE')
  allocate(LFd_LAE(nz_end), stat=ier); call CheckIerr(ier, trim(cerr)//' LFd_LAE')
  DO i = 1, nz_end
     N1 = 2; N2 = 4; N3 = nwave; Nbin = 200
     LF_LAE(i)%Nbin  = Nbin; LF_LAE(i)%N1  = N1; LF_LAE(i)%N2  = N2; LF_LAE(i)%N3  = N3
     LFd_LAE(i)%Nbin = Nbin; LFd_LAE(i)%N1 = N1; LFd_LAE(i)%N2 = N2; LFd_LAE(i)%N3 = N3
     allocate(LF_LAE(i)%bin(Nbin),  stat=ier); call CheckIerr(ier, &
          trim(cerr)//' LF_LAE(i)%bin')
     allocate(LFd_LAE(i)%bin(Nbin), stat=ier); call CheckIerr(ier, &
          trim(cerr)//' LFd_LAE(i)%bin')
     allocate(LF_LAE(i)%n(N1,N2,N3,Nbin),  stat=ier); call CheckIerr(ier, &
          trim(cerr)//' LF_LAE(i)%n')
     allocate(LFd_LAE(i)%n(N1,N2,N3,Nbin), stat=ier); call CheckIerr(ier, &
          trim(cerr)//' LFd_LAE(i)%n')
  ENDDO

  ! --- CSFRD
  allocate(CSFRD(nz_end, 2, 3), stat=ier); call CheckIerr(ier, trim(cerr)//' CSFRD')
  CSFRD(:,:,:) = 0.d0
END SUBROUTINE AllocateDistributionFunctions_LAE
! ========================================================================
SUBROUTINE DeallocateDistributionFunctions_LAE(nz_end)
  use LAErelated
  implicit none
  INTEGER, INTENT(IN) :: nz_end
  INTEGER :: iz, ier
  CHARACTER(LEN=50) :: cerr = '# (LAE)DeallocLAE: fail to deallocate'

  DO iz = 1, nz_end
     deallocate(LF_LAE(iz)%bin,  stat=ier); call CheckIerr(ier, &
          trim(cerr)//' LF_LAE(iz)%bin')
     deallocate(LFd_LAE(iz)%bin, stat=ier); call CheckIerr(ier, &
          trim(cerr)//' LFd_LAE(iz)%bin')
     deallocate(LF_LAE(iz)%n,  stat=ier); call CheckIerr(ier, trim(cerr)//' LF_LAE(iz)%n')
     deallocate(LFd_LAE(iz)%n, stat=ier); call CheckIerr(ier, trim(cerr)//' LFd_LAE(iz)%n')
  ENDDO
  deallocate(LF_LAE,  stat=ier); call CheckIerr(ier, trim(cerr)//' LF_LAE')
  deallocate(LFd_LAE, stat=ier); call CheckIerr(ier, trim(cerr)//' LFd_LAE')

  deallocate(CSFRD,   stat=ier); call CheckIerr(ier, trim(cerr)//' CSFRD')
END SUBROUTINE DeallocateDistributionFunctions_LAE
! ========================================================================
SUBROUTINE AllocSchaererArrays
  use LAErelated
  implicit none
  INTEGER :: i, ier, NT, NZ, NW
  CHARACTER(LEN=60) :: cerr = '# (LAE)AllocSchaererArrays: fail to allocate S03base(i)%'

  DO i = 1, 2
     NT = NT_S03; NZ = NZ_S03(i); NW = NWAVE_S03
     allocate(S03base(i)%time(NT),      stat=ier); call CheckIerr(ier, trim(cerr)//'time')
     allocate(S03base(i)%metal(NZ),     stat=ier); call CheckIerr(ier, trim(cerr)//'metal')
     allocate(S03base(i)%Lis(NW,NZ,NT), stat=ier); call CheckIerr(ier, trim(cerr)//'Lis')
     allocate(S03base(i)%Lcs(NW,NZ,NT), stat=ier); call CheckIerr(ier, trim(cerr)//'Lcs')

     ! --- initialize
     IF(i == 1) THEN
        S03base(i)%tstep = 1.d+5; S03base(i)%base = 9.d+4
     ELSEIF(i == 2) THEN
        S03base(i)%tstep = 1.d+6; S03base(i)%base = 9.9d+5
     ENDIF
     S03base(i)%invtstep  = 1.d0 / S03base(i)%tstep
     S03base(i)%base_iS03 = S03base(i)%base * S03base(i)%invtstep

     S03base(i)%time(:)    = 0.d0; S03base(i)%metal(:)   = 0.d0
     S03base(i)%Lis(:,:,:) = 0.d0; S03base(i)%Lcs(:,:,:) = 0.d0
  ENDDO
END SUBROUTINE AllocSchaererArrays
! ========================================================================
SUBROUTINE DeallocSchaererArrays
  use LAErelated
  implicit none
  INTEGER :: i, ier
  CHARACTER(LEN=70) :: cerr = '# (LAE)DeallocSchaererArrays: fail to deallocate S03base(i)%'

  DO i = 1, 2
     deallocate(S03base(i)%time,  stat=ier); call CheckIerr(ier, trim(cerr)//'time')
     deallocate(S03base(i)%metal, stat=ier); call CheckIerr(ier, trim(cerr)//'metal')
     deallocate(S03base(i)%Lis, stat=ier); call CheckIerr(ier, trim(cerr)//'Lis')
     deallocate(S03base(i)%Lcs, stat=ier); call CheckIerr(ier, trim(cerr)//'Lcs')
  ENDDO
END SUBROUTINE DeallocSchaererArrays
! ========================================================================
SUBROUTINE ReadSchaererSSP
  use LAErelated
  implicit none
  INTEGER :: i, j, jj, k, kend, loopend, istype
  CHARACTER(LEN=100) :: fbase, fname, fbase0 = 'inpdata/Schaerer03/'
  CHARACTER(LEN=50)  :: fmiddle(NZ_S03(1)+NZ_S03(2)) = &
       (/CHARACTER(LEN=50):: 'pop3_', 'e-70_', 'e-50_', 'e0004_', 'e001_', &
         'e004_', 'e008_', 'e020_', 'e040_', 'pop3_', 'e0004_', 'e004_', 'e020_'/)
  DOUBLE PRECISION, PARAMETER :: EPS = 1.d-8
  DOUBLE PRECISION :: x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12
  DOUBLE PRECISION :: x, NZ1

  DO i = 1, NT_S03
     S03base(1)%time(i) = S03base(1)%tstep * dble(i) - S03base(1)%base
     S03base(2)%time(i) = S03base(2)%tstep * dble(i) - S03base(2)%base
  ENDDO
  S03base(1)%metal(1) =  0.d0; S03base(2)%metal(1) = 0.d0
  S03base(1)%metal(2) = 1.d-7
  S03base(1)%metal(3) = 1.d-5
  S03base(1)%metal(4) = 4.d-4; S03base(2)%metal(2) = 4.d-4
  S03base(1)%metal(5) = 1.d-3
  S03base(1)%metal(6) = 4.d-3; S03base(2)%metal(3) = 4.d-3
  S03base(1)%metal(7) = 8.d-3
  S03base(1)%metal(8) = 2.d-2; S03base(2)%metal(4) = 2.d-2
  S03base(1)%metal(9) = 4.d-2

  NZ1 = NZ_S03(1); loopend = NZ_S03(1) + NZ_S03(2)
  DO i = 1, 3
     IF(i == 1) THEN
        fbase = trim(fbase0)//'QH_'
     ELSEIF(i == 2) THEN
        fbase = trim(fbase0)//'L_UV_'
     ELSEIF(i == 3) THEN
        fbase = trim(fbase0)//'Lcont_'
     ENDIF

     DO j = 1, loopend
        fname = trim(fbase)//trim(fmiddle(j))
        IF(j <= NZ1) THEN
           istype = 1; jj = j
           call ReadLineNumbers(istype, fname, kend)
           open(1, file = trim(fname)//'is1.dat', status = 'old')
           open(2, file = trim(fname)//'cs1.dat', status = 'old')
        ELSE
           istype = 2; jj = j - NZ1
           call ReadLineNumbers(istype, fname, kend)
           open(1, file = trim(fname)//'is5.dat', status = 'old')
           open(2, file = trim(fname)//'cs5.dat', status = 'old')
        ENDIF

        IF(i == 1) THEN
           ! read N_LyC data
           read(1, *); read(2, *) ! skip the header
           DO k = 1, kend
              x1 = 0.d0; x2 = 0.d0
              read(1, *) x, x, x1; read(2, *) x, x, x2
              ! (1)Ncolumn (2)log[age/yr] (3)QH[s^-1]
              ! Lis(1,:,:), Lcs(1,:,:) --- L(Lya) [erg/s]
              ! L(Lya)/(erg/s) = [h_P*nu_Lya/erg]*[N(Lya)/QH]*[QH/s^-1]
              !                = [fLya/erg] * [QH/s^-1]
              S03base(istype)%Lis(1,jj,k) = constLAE%fLya * 10.d0**x1
              S03base(istype)%Lcs(1,jj,k) = constLAE%fLya * 10.d0**x2
           ENDDO
           close(1); close(2)
        ELSEIF(i == 2) THEN
           ! read L1500 and L2800 data
           DO k = 1, 5 ! skip the header
              read(1, *); read(2, *)
           ENDDO
           DO k = 1, kend
              x1 = 0.d0; x2 = 0.d0; x3 = 0.d0; x4 = 0.d0; x5 = 0.d0; x6 = 0.d0
              read(1, *) x, x, x1, x2, x5; read(2, *) x, x, x3, x4, x6
              ! (1)Ncolumn (2)log[age/yr]
              ! (3,4)log[Llam(1500)/(erg/s/A)] ave. over +- 20A,
              !      150A band, respectively
              ! (5)log[Llam(2800)/(erg/s/A)] ave. over +- 280A band
              ! Lis(2,:,:), Lcs(2,:,:) --- Llam(1500) [erg/s/A]
              ! Lis(3,:,:), Lcs(3,:,:) --- Llam(2800) [erg/s/A]
              S03base(istype)%Lis(2,jj,k) = 10.d0**x2
              S03base(istype)%Lcs(2,jj,k) = 10.d0**x4
              S03base(istype)%Lis(3,jj,k) = 10.d0**x5
              S03base(istype)%Lcs(3,jj,k) = 10.d0**x6
           ENDDO
           close(1); close(2)
        ELSEIF(i == 3) THEN
           ! read L6563 and L4861 data
           read(1, *); read(2, *) ! skip the header
           DO k = 1, kend
              x1  = 0.d0;  x2 = 0.d0; x3 = 0.d0; x4 = 0.d0;  x5 = 0.d0
              x6  = 0.d0;  x7 = 0.d0; x8 = 0.d0; x9 = 0.d0; x10 = 0.d0
              x11 = 0.d0; x12 = 0.d0
              read(1, *) x, x, x1, x2, x3, x4, x9,  x10
              read(2, *) x, x, x5, x6, x7, x8, x11, x12
              ! (1)Ncolum    (2)log[age/year] (3)x{1,5}=L(Hb)[erg/s]
              ! (4)x{2,6}=EW(Ha)[A]  (5)x{3,7}=L(Ha)/L(Hb) (6)x{4,8}=EW(Hb)[A]
              ! (7)x{5,9}=EW(Lya)[A] (8)x{6,10}=L(Lya)/L(Hb)
              ! Lis(4,:,:), Lcs(4,:,:) --- Llam(6543) [erg/s/A]
              ! Lis(5,:,:), Lcs(5,:,:) --- Llam(4861) [erg/s/A]
              ! Lis(6,:,:), Lcs(6,:,:) --- Llam(1216) [erg/s/A]
              ! Llam(6563) = [L(Ha)/L(Hb)*L(Hb)]/EW(Ha): x3*x1/x2,
              !                                          x7*x5/x6
              ! Llam(4861) = L(Hb)/EW(Hb): x1/x4, x5/x8
              ! Llam(1216) = [L(Lya)/L(Hb)*L(Hb)]/EW(Lya): x10*x1/x9,
              !                                            x12*x5/x11
              IF(x2 < EPS) THEN ! x2 = EW(Ha)
                 S03base(istype)%Lis(4,jj,k) = 0.d0
              ELSE
                 S03base(istype)%Lis(4,jj,k) = x3 * x1 / x2
              ENDIF
              S03base(istype)%Lcs(4,jj,k) = x7 * x5 / x6
              IF(x4 < EPS) THEN ! x4 = EW(Hb)
                 S03base(istype)%Lis(5,jj,k) = 0.d0
              ELSE
                 S03base(istype)%Lis(5,jj,k) = x1 / x4
              ENDIF
              S03base(istype)%Lcs(5,jj,k) = x5 / x8
              IF(x9 < EPS) THEN ! x9 = EW(Lya)
                 S03base(istype)%Lis(6,jj,k) = 0.d0
              ELSE
                 S03base(istype)%Lis(6,jj,k) = x10 * x1 / x9
              ENDIF
              S03base(istype)%Lcs(6,jj,k) = x12 * x5 / x11
           ENDDO
           close(1); close(2)
        ENDIF
     ENDDO
  ENDDO
CONTAINS
! -----------------------------------------------------------------------
  SUBROUTINE ReadLineNumbers(istype, filename, linenum)
    CHARACTER(LEN=100), INTENT(IN) :: filename
    INTEGER, INTENT(IN) :: istype
    INTEGER, INTENT(INOUT) :: linenum
    INTEGER :: error
    CHARACTER(LEN=100) :: buf

    IF(istype == 1) THEN
       open(1, file = trim(filename)//'is1.dat', status = 'old')
    ELSEIF(istype == 2) THEN
       open(1, file = trim(filename)//'is5.dat', status = 'old')
    ENDIF
    linenum = 0; error = 0
    DO WHILE(error == 0)
       read(1, *, iostat=error) buf
       IF(error == 0 .and. buf(1:1) /= '#') linenum = linenum + 1
    ENDDO
    close(1)
  END SUBROUTINE ReadLineNumbers
! -----------------------------------------------------------------------
END SUBROUTINE ReadSchaererSSP
! ========================================================================
SUBROUTINE CorrectionSchaererSSP(inode)
  use LAErelated
  implicit none
  INTEGER, INTENT(IN) :: inode ! for parallel
  INTEGER, PARAMETER :: flag_correction = 0
                        ! 0: correction made for all (original)
                        ! 1:    only for metal-rich (Z>=Zcrit)
!!$        DOUBLE PRECISION, PARAMETER :: corr_QH = 0.35666d0
!!$        DOUBLE PRECISION, PARAMETER :: corr_UV = 0.38163d0
  DOUBLE PRECISION, PARAMETER :: corr = 0.392544d0 
                                 ! int_1^100[m*phi(m)]dm
                                 !  where phi(m)=phi0*m^-alpha
                                 !  for ml <= m <= mu (ml=0.1,mu=100) and
                                 !  alpha = 2.35
  DOUBLE PRECISION, PARAMETER :: Zcrit = 2.d-7 ! Zcrit=10^-5*Zsun=2e-7
  INTEGER :: i, NZ

  IF(inode == 0) print '(A, I1, $)', &
       '# (LAE)flag_correction = ', flag_correction
  IF(flag_correction == 0) THEN ! original
     IF(inode == 0) write(*, *)
     DO i = 1, 2
        S03base(i)%Lis(:,:,:) = corr * S03base(i)%Lis(:,:,:)
        S03base(i)%Lcs(:,:,:) = corr * S03base(i)%Lcs(:,:,:)
     ENDDO
  ELSEIF(flag_correction == 1) THEN ! only for metal-rich (Z>=10^-5)
     IF(inode == 0) print '(A, G9.2, A)', &
          ': correction made only for Z>=Zcrit (=', Zcrit, ')'
     NZ = NZ_S03(1)
     DO i = 1, NZ
        IF(S03base(1)%metal(i) >= Zcrit) THEN
           S03base(1)%Lis(:,i,:) = corr * S03base(1)%Lis(:,i,:)
           S03base(1)%Lcs(:,i,:) = corr * S03base(1)%Lcs(:,i,:)
        ENDIF
     ENDDO
     NZ = NZ_S03(2)
     DO i = 1, NZ
        IF(S03base(2)%metal(i) >= Zcrit) THEN
           S03base(2)%Lis(:,i,:) = corr * S03base(2)%Lis(:,i,:)
           S03base(2)%Lcs(:,i,:) = corr * S03base(2)%Lcs(:,i,:)
        ENDIF
     ENDDO
  ENDIF
END SUBROUTINE CorrectionSchaererSSP
! ========================================================================
SUBROUTINE InitializeArraysForLAE
  use global_var; use LAErelated
  implicit none
  INTEGER :: i, iband, ier
  CHARACTER(LEN=60) :: cerr = '# (LAE)InitializeArraysForLAE: fail to allocate LAE%'
  DOUBLE PRECISION :: tmp
  DOUBLE PRECISION :: Square ! function

  ! 1:Ngas[cm^-2], 2:Zc[Zsun], 3:NcZc[Zsun/cm^2], 4:Mstar[Msun],
  ! 5:Mhost[Msun], 6:reff[kpc/h], 7:SFR[Msun/yr], 8:taueff[yr]
  DOUBLE PRECISION, PARAMETER :: base_distLAE(n_distLAE) &
       = (/10.d0, -10.d0, 5.d0, 0.d0, 0.d0, -10.d0, -10.d0, 0.d0/)
  DOUBLE PRECISION, PARAMETER :: step_distLAE = 0.1d0 ! NLF1=200

  ! 1:Lya-Ngas, 2:Lya-Z, 3:Lya-Mstar, 4:Lya-NZ
  DOUBLE PRECISION, PARAMETER :: base_distLya(n_distLya) &
       = (/17.5d0, -4.5d0, 2.d0, 13.d0/)
  DOUBLE PRECISION, PARAMETER :: step_distLya(n_distLya) &
       = (/0.5d0, 0.5d0, 1.d0, 1.d0/) ! NLF2 = 10

  constLAE%Zsun  = const%Zsun

  constLAE%Ngas = 1.25d+14 * Square(param%h)
                  ! --- a conversion factor from [M_sun h^2 kpc^-2] to [cm^-2]
                  !      : 1.0 [M_sun h^2 kpc^-2] * 1.989e33 [g/M_sun] 
                  !          / 1.673e-24 [g] / (3.08e21/h [cm/ h^-1 kpc])^2
                  !      = 1.253e+14 [cm^-2]
  constLAE%Ngas = constLAE%Ngas * 0.5d0 / const%PI

  tmp = 1.d-14 * Square(param%h) / constLAE%Ngas
        ! --- a conversion factor from [cm^-2] to [10^14 Msun/kpc^2]
        !     : 1.0 [cm^-2] / const%Ngas [cm^-2/(Msun h^2/kpc^2)] * 10^-14
  constLAE%OptV = tmp * param%tauV0 / const%Zsun
  constLAE%OptB = tmp * param%tauV0 * xi(param%iBband_r) / const%Zsun

  IF(paramLAE%exttype == 1) THEN ! screen dust geometry
     IF(inode == 0) print '(A)', '# (LAE)screen dust geometry is assumed'
     constLAE%AV0 = 2.5d0 * log10(exp(1.d0))
     constLAE%AB0 = constLAE%AV0 * xi(param%iBband_r)
  ELSEIF(paramLAE%exttype == 2) THEN ! slab dust geometry
     IF(inode == 0) print '(A)', '# (LAE)slab dust geometry is assumed'
  ENDIF

  constLAE%Tdyn = 9.7775d+8 / param%h
                  ! conversion factor from [kpc/h/(km/s)] into [yr]
  constLAE%Tdyn = paramLAE%Rtburst * constLAE%Tdyn
                  ! convert the time duration of a burst star formation from the
                  !  dynamical timescale 'tdyn' into 'paramLAE%Rtburst * tdyn'

  constLAE%convH2 = Square(param%h) ! conversion factor from [erg/s] into [erg/s/h^2]

  iLAE%numLya = 0; iLAE%numUV = 0

  allocate(LAE%flag_mag(param%nwave), stat=ier); call CheckIerr(ier, &
       trim(cerr)//'flag_mag')
  LAE%flag_mag(:) = 0
  DO iband = 1, param%nwave
     IF((param%iwave(iband) >= 50 .and. param%iwave(iband) <= 58) .or. &
          (param%iwave(iband) >= 50+param%nwave .and. &
           param%iwave(iband) <= 58+param%nwave)) LAE%flag_mag(iband) = 1
  ENDDO
  allocate(LAE%lumg(param%tnw),    stat=ier); call CheckIerr(ier, trim(cerr)//'lumg')
  allocate(LAE%lumg_d(param%tnw),  stat=ier); call CheckIerr(ier, trim(cerr)//'lumg_d')
  allocate(LAE%mag(param%nwave),   stat=ier); call CheckIerr(ier, trim(cerr)//'mag')
  allocate(LAE%mag_d(param%nwave), stat=ier); call CheckIerr(ier, trim(cerr)//'mag_d')
  LAE%lumg(:) = 0.d0; LAE%lumg_d(:) = 0.d0; LAE%mag(:) = 0.d0; LAE%mag_d(:) = 0.d0
  LAE%LLya = 0.d0; LAE%LS03(:) = 0.d0; LAE%LS03d(:) = 0.d0; LAE%LS03pre(:) = 0.d0
  LAE%magS03(:) = 0.d0; LAE%magS03d(:) = 0.d0; xi_S03(:) = 0.d0
  ! --- modified by MARK (2018/Mar/02)
!!$  corrS03(2) = 81.90637d0 + const%corr
!!$               ! -2.5log[(lam^2/c)/4pi/(10pc)^2]-48.6 at lam = 1500 A
!!$  corrS03(1) = corrS03(2) - 5.d0 * log10(1.21567d0/1.5d0)
!!$               ! -2.5log[(lam^2/c)/4pi/(10pc)^2]-48.6 at lam = 1216 A
!!$  corrS03(3) = corrS03(2) - 5.d0 * log10(2.8d0/1.5d0)
!!$               ! -2.5log[(lam^2/c)/4pi/(10pc)^2]-48.6 at lam = 2800 A
!!$  corrS03(4) = corrS03(2) - 5.d0 * log10(6.563d0/1.5d0)
!!$               ! -2.5log[(lam^2/c)/4pi/(10pc)^2]-48.6 at lam = 6563 A
!!$  corrS03(5) = corrS03(2) - 5.d0 * log10(4.861d0/1.5d0)
!!$               ! -2.5log[(lam^2/c)/4pi/(10pc)^2]-48.6 at lam = 4861 A
!!$  corrS03(6) = corrS03(1) ! at lam = 1216 A
  corrS03(1) = const%corr_wdim(1) ! lam =  912A (unused)
  corrS03(2) = const%corr_wdim(4) ! lam = 1500A
  corrS03(3) = const%corr_wdim(7) ! lam = 2800A
  corrS03(4) = const%corr_wdim(9) ! lam = 6563A
  corrS03(5) = const%corr_wdim(8) ! lam = 4861A
  corrS03(6) = const%corr_wdim(2) ! lam = 1216A
  DO iband = 1, NWAVE_S03
     corrS03(iband) = corrS03(iband) + const%corr
  ENDDO
 
  LyaLF%count(:,:,:) = 0.d0; LyaLF%wcount(:,:,:) = 0.d0
  DO i = 1, n_distLAE
     distLAE(i)%count(:,:,:) = 0.d0; distLAE(i)%wcount(:,:,:) = 0.d0
     distLAE(i)%base = base_distLAE(i); distLAE(i)%step = step_distLAE
     distLAE(i)%invstep = 1.d0 / distLAE(i)%step
  ENDDO
  DO i = 1, n_distLya
     distLya(i)%count(:,:,:) = 0.d0
     distLya(i)%base = base_distLya(i); distLya(i)%step = step_distLya(i)
     distLya(i)%invstep = 1.d0 / distLya(i)%step
  ENDDO
  DO i = 1, NWAVE_S03
     LF_S03(i)%count(:,:,:) = 0.d0; LF_S03(i)%countd(:,:,:) = 0.d0
  ENDDO
END SUBROUTINE InitializeArraysForLAE
! ========================================================================
SUBROUTINE optdepthForLAE(ilog)
  use global_var; use LAErelated
  implicit none
  INTEGER, INTENT(IN) :: ilog
  INTEGER :: iband
  DOUBLE PRECISION :: x, xiV, lamV = 0.55121d0 ! (central wavelength of V-band [um])
  DOUBLE PRECISION :: xp(NWAVE_S03) = (/&
       0.912d0, 0.15d0, 0.28d0, 0.6563d0, 0.4861d0, 0.121567d0/)
  DOUBLE PRECISION :: Calzetti00, Pei92 ! functions

  ! --- determine xiV first
  IF(param%exttype == 4) THEN
     xiV = Calzetti00(lamV)
  ELSE
     xiV = Pei92(lamV)
  ENDIF

  DO iband = 1, NWAVE_S03
     x = xp(iband)
     IF(param%exttype == 4) THEN
        xi_S03(iband) = Calzetti00(x)
     ELSE
        xi_S03(iband) = Pei92(x)
     ENDIF
     IF(inode == 0) THEN
        print '(A, F8.5, A, F10.4)', '   --- xi_S03(lam_rest=', xp(iband), &
             '[um]) =', xi_S03(iband)
        write(ilog, '(A, F8.5, A, F10.4)') '# --- xi_S03(lam_rest=', xp(iband), &
             '[um]) =', xi_S03(iband)
     ENDIF

     xi_S03(iband) = xi_S03(iband) / xiV
  ENDDO
END SUBROUTINE optdepthForLAE
! ========================================================================
SUBROUTINE EraseGalaxyForLAE(id)
  use LAErelated
  INTEGER, INTENT(IN) :: id

  allgal(id)%num_b  = 0;    allgal(id)%flag_c  = 0
  allgal(id)%SFR(:) = 0.d0; allgal(id)%beta(:) = 0.d0
  allgal(id)%Zc     = 0.d0; allgal(id)%Zc_pre  = 0.d0
  allgal(id)%Mc_pre = 0.d0; allgal(id)%Ms_pre  = 0.d0
  allgal(id)%rb_pre = 0.d0; allgal(id)%Vb_pre  = 0.d0
  allgal(id)%taueff(:,:) = 0.d0
  allgal(id)%lumg_pre(:) = 0.d0
END SUBROUTINE EraseGalaxyForLAE
! ========================================================================
SUBROUTINE CalAllGalForLAE(b_or_q, id, Mcool0, zp1pc)
  use global_var; use LAErelated
  implicit none
  INTEGER, INTENT(IN) :: b_or_q, id ! b_or_q = 1:starburst, 2:quiescent
  DOUBLE PRECISION, INTENT(IN) :: Mcool0 ![Msun]
  DOUBLE PRECISION, INTENT(IN) :: zp1pc
  DOUBLE PRECISION, PARAMETER :: EPS = 1.d-10
  DOUBLE PRECISION :: taueff
  INTEGER :: nb

  IF(abs(zp1pc - param%zsp1) <= EPS .and. tlife > 0.d0) THEN
     allgal(id)%flag_c = gal(id)%flag_c
     IF(b_or_q == 2) THEN ! quiescent
        allgal(id)%SFR(b_or_q) = allgal(id)%SFR(b_or_q) &
             + (dMstar/ssp%alp(b_or_q) * param%munit) / tlife_yr
        allgal(id)%Zc = Zc; allgal(id)%beta(b_or_q) = beta
        allgal(id)%taueff(b_or_q, 1) = taust * param%th_yr / (ssp%alp(b_or_q) + beta)
     ELSEIF(b_or_q == 1) THEN ! starburst
        ! --- the time duration of the starburst is assumed
        !      to be the dynamical time-scale
        allgal(id)%num_b = allgal(id)%num_b + 1
        nb = allgal(id)%num_b; taueff = constLAE%Tdyn * gal(id)%rbulge / gal(id)%Vbulge
        IF(nb <= 2) allgal(id)%taueff(b_or_q, nb) = taueff
                    ! dynamical timescale [yr] estimated by bulge radius [kpc/h]
                    !  and velocity dispersion of bulge [km/s]
        allgal(id)%SFR(b_or_q) = allgal(id)%SFR(b_or_q) &
             + dMstar/ssp%alp(b_or_q) * param%munit / (taueff * (ssp%alp(b_or_q) + beta))
        IF(nb == 1) THEN
           allgal(id)%Zc      = Zc ! = 0.d0
           allgal(id)%Zc_pre  = Zc0
           allgal(id)%Mc_pre  = Mcool0 ! [M_sun]
           allgal(id)%Ms_pre  = Ms0
           allgal(id)%rb_pre  = gal(id)%rbulge
           allgal(id)%beta(b_or_q) = beta
           allgal(id)%Vb_pre       = gal(id)%Vbulge
           allgal(id)%lumg_pre(:)  = gal(id)%lumg(:)
        ENDIF
     ENDIF
  ENDIF
END SUBROUTINE CalAllGalForLAE
! ========================================================================
SUBROUTINE MainCalcForLAE(id, nz, iforest, Mhost, weight, rb0, rd0, mor, mord, Mstar)
  use global_var; use LAErelated
  implicit none
  INTEGER, INTENT(IN) :: id, nz, iforest
  INTEGER, INTENT(IN) :: mor, mord
  DOUBLE PRECISION, INTENT(IN) :: Mhost, Mstar ! [param%munit]
  DOUBLE PRECISION, INTENT(IN) :: weight ! [h^3/Mpc^3]
  DOUBLE PRECISION, INTENT(IN) :: rb0, rd0 ! [kpc/h]
  INTEGER :: iloop, bin, type, morLAE, mordLAE
  INTEGER :: iband
  LOGICAL :: logic_dust
  DOUBLE PRECISION :: zp1pre
  DOUBLE PRECISION :: lumtmp(param%tnw), lumtmp_d(param%tnw)
  DOUBLE PRECISION :: lumtmp2(param%nwave), mz, mt
  DOUBLE PRECISION :: rb0_pre, t_brst_ratio, lya_phase, fdecline
  DOUBLE PRECISION :: SFRini, taust_inv, tGWh, temp
  ! --- functions
  INTEGER :: CalBin
  DOUBLE PRECISION :: z2t, CalNgas, CalBulgeReff


!!$   call corrIGMopt(IGM_corr, mag, mag_d) ! include the IGM opacity effect
!!$  zp1pre     = (1.d0 + param%delz00) * param%zsp1
  zp1pre = mrgp%zp1ar(param%izout-1)
  tlife = z2t(param%zsp1) - z2t(zp1pre)
  tlife = tlife * param%th_yr ! [hubble time] ==> [yr]

  LAE%id = id; LAE%flag_c = gal(id)%flag_c; morLAE = mor; mordLAE = mord
  LAE%Vc = gal(id)%Vcent; LAE%Vcc = gal(gal(id)%id_cgal)%Vcent ! w/ dyn. resp.
!!$  LAE%Vc = gal(id)%Vc; LAE%Vcc = gal(gal(id)%id_cgal)%Vc ! w/o dyn. resp.
  LAE%Mhost = Mhost * param%munit; LAE%Mprog = gal(id)%Mhalo * param%munit
  rb0_pre = CalBulgeReff(allgal(id)%rb_pre)
!!$  IF(LAE%flag_c == 1) print '(I2, 2I10, 4G10.3)', LAE%flag_c, LAE%id, &
!!$       gal(LAE%id)%id_cgal, allgal(LAE%id)%SFR(1), allgal(LAE%id)%SFR(2), &
!!$       gal(LAE%id)%SFR, gal(LAE%id)%Mcool * param%munit
!-------------------------------------------------------------------------------
  IF(allgal(LAE%id)%SFR(1) > 0.d0 .and. allgal(LAE%id)%SFR(2) > 0.d0) THEN
     !---------------------------------------------------------
     !     for multiple-merging galaxies
     !---------------------------------------------------------
     LAE%b_or_q = 1; LAE%sftype = 3
     LAE%Ms_pre = allgal(LAE%id)%Ms_pre * param%munit ! [Msun]
     LAE%Vel    = allgal(LAE%id)%Vb_pre

     ! -------------------------------------------------------------
     !  add the contribution from quiescent star-formation activity
     !   into the disk luminosity of gal(i)%lumg(param%nwp1:param%tnw)
     !   if quiescent SF occurred after starburst
     !  add the contribution from quiescent star-formation activity
     !   into L(Lya) and L_lam(UV)
     !  the detail of this treatment does NOT significantly contribute
     !   to the total L(Lya) and L_lam(UV)
     ! -------------------------------------------------------------
     IF(gal(LAE%id)%Mcoold > 0.d0) THEN ! quiescent SF occurred after starburst
                                       !  (starburst --> quiescent)
        type  = 0
        LAE%Mgas_q = gal(LAE%id)%Mcoold * param%munit
                     ! gas mass remained after quiescent SF [Msun]
        LAE%McZc_q = LAE%Mgas_q * allgal(LAE%id)%Zc
                     ! cold metal mass remained after quiescent SF [Msun]
        LAE%Ms_q   = allgal(LAE%id)%SFR(2) * ssp%alp(2) * tlife
                     ! mass of the stars produced in quiescent SF [Msun]
        beta  = allgal(LAE%id)%beta(2)
        taust = allgal(LAE%id)%taueff(2, 1) * (ssp%alp(2) + beta) ![yr]
        Zc0   = allgal(LAE%id)%Zc - ssp%p(2) * tlife / taust
        taust = taust / param%th_yr ! [hubble time]
        call lum(param%zsp1, zp1pre, lumtmp2, mz, mt, 0.d0)
        temp = tlife / allgal(LAE%id)%taueff(2, 1)
        temp = gal(LAE%id)%Mcoold * exp(temp) * param%munit ! Mc0 [Msun]
        LAE%lumg(param%nwp1:param%tnw) = lumtmp2(1:param%nwave) * temp

        temp = Zc0 / const%Zsun
!!$          temp = allgal(LAE%id)%Zc / const%Zsun
     ELSE ! the quiescent SF occurred prior to starburst (quiescent --> starburst)
        type = 3
        ! metallicity is estimated as metal_eff prior to the starburst
        !  therefore, QH_is5 and LAE%LS03(1)_is is used to estimate
        !  L(Lya) and L(UV)
!!$        gal(LAE%id)%lumg(param%nwp1:param%tnw) = 0.d0
        LAE%lumg(param%nwp1:param%tnw) = 0.d0

        temp = allgal(LAE%id)%Zc_pre / const%Zsun
     ENDIF
     call CalPreLS03(tlife, temp, allgal(LAE%id)%taueff(2, 1), &
                     allgal(LAE%id)%SFR(2))

     call InitializeQuantitiesForBurst
     ! --- for Lya-bright galaxies
     ID_iloop1: DO iloop = 1, paramLAE%nloop
        call CalTelapsSFRandLS03(lya_phase)
        IF(type /= 3) THEN ! starburst --> quiescent
           LAE%SFR = LAE%SFR + allgal(LAE%id)%SFR(2)
        ENDIF

        bin = 0
        IF(LAE%LLya > 0.d0) bin = CalBin(log10(LAE%LLya), 10.d0, 30.d0)
        IF(bin <= NLF1 .and. bin >= 1) THEN
           call CalMultipleMerger(logic_dust)

           call SubroutineSetForLyaBrightGalaxies
        ENDIF
     ENDDO ID_iloop1

     ! --- for UV-bright galaxies
     LAE%weightAll = weight / dble(paramLAE%nloop)
     ID_iloop2: DO iloop = 1, paramLAE%nloop
        call CalTelapsSFRandLS03(tlife)
        IF(type /= 3) THEN ! starburst --> quiescent
           LAE%SFR = LAE%SFR + allgal(LAE%id)%SFR(2)
        ENDIF

        call CalMultipleMerger(logic_dust)
        call SubroutineSetForUVBrightGalaxies
     ENDDO ID_iloop2
  ELSEIF(allgal(LAE%id)%SFR(2) > 0.d0) THEN
     !------------------------------------------------------
     !     for quiescently star forming galaxies
     !------------------------------------------------------
     LAE%b_or_q   = 2; LAE%sftype = 1; beta = allgal(LAE%id)%beta(LAE%b_or_q)
     LAE%beta     = beta; LAE%ab = ssp%alp(LAE%b_or_q) + beta
     LAE%weight   = weight; LAE%weightAll = weight ! [h^3/Mpc^3]
     LAE%SFR      = allgal(LAE%id)%SFR(LAE%b_or_q) ! [Msun/yr]
     LAE%taust_yr = allgal(LAE%id)%taueff(LAE%b_or_q, 1) ! [yr]
     LAE%tel      = tlife; LAE%twind = 1.d+11
     LAE%Mstar    = Mstar * param%munit ! [Msun]
     LAE%Mgas     = gal(LAE%id)%Mcoold * param%munit ! [Msun]
     LAE%Mgas_pre = LAE%Mgas * exp(LAE%tel/LAE%taust_yr)
     LAE%Ms_pre   = LAE%Mstar - (ssp%alp(LAE%b_or_q) / LAE%ab) &
                                * (LAE%Mgas_pre - LAE%Mgas)
     IF(LAE%Ms_pre < 0.d0) LAE%Ms_pre = 0.d0
     LAE%Zc       = allgal(LAE%id)%Zc / const%Zsun ! [Zsun]
     LAE%mag(:)   = mag(:); LAE%mag_d(:) = mag_d(:)
     DO iband = 1, param%nwave
        IF(LAE%flag_mag(iband) == 1) THEN
           ! LAE%mag(:), LAE%mag_d(:)
           ! = log10[NLyC/(photons/s/h^2) or log10[Llam/(erg/s/A/h^2)]
           LAE%mag(iband)   = -0.4d0 * LAE%mag(iband)
           LAE%mag_d(iband) = -0.4d0 * LAE%mag_d(iband)
        ENDIF
     ENDDO
     IF(mordLAE == 1 .or. mordLAE == 2) THEN ! E or S0
        LAE%Vel = gal(LAE%id)%Vbulge; LAE%reff = rb0
     ELSE ! S
        LAE%Vel = gal(LAE%id)%Vdisk; LAE%reff = rd0
     ENDIF
     LAE%Ngas = CalNgas(LAE%Mgas, LAE%reff) ! [cm^-2]
     LAE%NZ   = LAE%Ngas * LAE%Zc * const%Zsun ! [cm^-2]
     LAE%Tmass = gal(LAE%id)%Tmass * param%th_yr / Mstar ! mass-weighted age [yr]

     taust  = LAE%taust_yr / param%th_yr ! tau_eff [hubble time]
     temp   = LAE%tel / LAE%taust_yr
     SFRini = LAE%SFR * temp / (1.d0 - exp(-temp))
              ! translate a mean SFR of 'LAE%SFR' in the time step
              !  into the initial SFR of 'SFRini'

     call CalLS03(LAE%sftype, LAE%tel, tlife, tlife, SFRini, LAE%Zc, &
                  LAE%LS03pre, LAE%LLya, LAE%NLyC)
     call extcont(LAE%NZ, LAE%LS03, LAE%LS03d)

     bin = 0
     IF(LAE%LLya > 0.d0) bin = CalBin(log10(LAE%LLya), 10.d0, 30.d0)
     ! --- for Lya-bright galaxies
     IF(bin <= NLF1 .and. bin >= 1) THEN
        call SubroutineSetForLyaBrightGalaxies
     ENDIF
     logic_dust = .true.
     ! --- for UV-bright galaxies
     call SubroutineSetForUVBrightGalaxies
  ELSEIF(allgal(LAE%id)%SFR(1) > 0.d0) THEN
     ! ------------------------------------------------------
     !     for bursting star forming galaxy
     ! ------------------------------------------------------
     LAE%b_or_q = 1; LAE%sftype = 2
     LAE%Ms_pre = allgal(LAE%id)%Ms_pre * param%munit ! [Msun]
     LAE%Vel    = allgal(LAE%id)%Vb_pre

     morLAE = 1; mordLAE = 1 ! morphology of starburst galaxy is always E

     call InitializeQuantitiesForBurst

     ! --- for Lya-bright galaxies
     ID_iloop3: DO iloop = 1, paramLAE%nloop
        call CalTelapsSFRandLS03(lya_phase)

        bin = 0
        IF(LAE%LLya > 0.d0) bin = CalBin(log10(LAE%LLya), 10.d0, 30.d0)
        IF(bin <= NLF1 .and. bin >= 1) THEN
           call CalStarburst(logic_dust)
           call SubroutineSetForLyaBrightGalaxies
        ENDIF
     ENDDO ID_iloop3

     ! --- for UV-bright galaxies
     LAE%weightAll = weight / dble(paramLAE%nloop)
     ID_iloop4: DO iloop = 1, paramLAE%nloop
        call CalTelapsSFRandLS03(tlife)

        call CalStarburst(logic_dust)
        call SubroutineSetForUVBrightGalaxies
     ENDDO ID_iloop4
  ENDIF
! -----------------------------------------------------------------------
CONTAINS
! -----------------------------------------------------------------------
  SUBROUTINE InitializeQuantitiesForBurst
    DOUBLE PRECISION :: CalTbrstRatio, CalLyaPhase, CalT_GW ! functions

    t_brst_ratio = CalTbrstRatio(allgal(LAE%id)%num_b,&
                                 allgal(LAE%id)%taueff(LAE%b_or_q,1),&
                                 allgal(LAE%id)%taueff(LAE%b_or_q,2),&
                                 tlife, LAE%id, allgal(LAE%id)%flag_c, &
                                 LAE%Mhost, weight)
    lya_phase = CalLyaPhase(t_brst_ratio, tlife)

    LAE%Mgas_pre = allgal(LAE%id)%Mc_pre ! [Msun]
    LAE%Zc_pre   = allgal(LAE%id)%Zc_pre / const%Zsun ! [Zsun]
    Zc0          = allgal(LAE%id)%Zc_pre ! metallicity (not in [Zsun])
    LAE%reff     = rb0_pre ! [kpc/h]
    beta         = allgal(LAE%id)%beta(LAE%b_or_q)
    LAE%beta     = beta; LAE%ab = ssp%alp(LAE%b_or_q) + beta
    LAE%Ms_b     = LAE%Mgas_pre * ssp%alp(LAE%b_or_q) / LAE%ab
    LAE%tweight  = weight * (lya_phase / tlife)    ! [h^3/Mpc^3]
                   ! the number density of galaxies is reduced by the factor
                   !   of 'lya_phase/tlife' according to the condition
                   !  that the galaxies are bright enough to be observed
    LAE%weight   = LAE%tweight / dble(paramLAE%nloop)
                   ! weight of a bursting galaxy [h^3/Mpc^3]
    IF(paramLAE%type_b == 1) THEN ! in the case of instantaneous burst
       taust     = 1.d-8 ! tau_star [hubble time]
       taust_inv = param%th_yr / taust ! [yr^-1]
       SFRini    = LAE%Mgas_pre * ssp%alp(LAE%b_or_q) / LAE%ab
    ELSEIF(paramLAE%type_b == 2) THEN ! in the case of exponential burst
       LAE%taust_yr = t_brst_ratio * tlife ! taustar[yr]
       taust        = LAE%taust_yr / param%th_yr ! taustar[hubble time]
       taust_inv    = 1.d0 / LAE%taust_yr ! [yr^-1]
       SFRini       = LAE%Mgas_pre * taust_inv
       LAE%twind    = CalT_GW(LAE%taust_yr, ssp%alp(LAE%b_or_q), LAE%beta)
       tGWh = LAE%twind / param%th_yr ! tGW [hubble time]
    ENDIF
    call CalBaseVarUsingMetal(LAE%Zc_pre)
  END SUBROUTINE InitializeQuantitiesForBurst
! -----------------------------------------------------------------------
  SUBROUTINE CalTelapsSFRandLS03(time)
    DOUBLE PRECISION, INTENT(IN) :: time
    DOUBLE PRECISION :: CalElapsedTime, CalFdecline, CalSFRburst

    LAE%tel  = CalElapsedTime(param%idum, time) ! [yr]
    fdecline = CalFdecline(paramLAE%type_b, LAE%tel, LAE%taust_yr)
    LAE%SFR  = CalSFRburst(LAE%tel, LAE%twind, SFRini, fdecline)
    call CalLS03(LAE%sftype, LAE%tel, LAE%twind, taust_inv, SFRini, &
                 LAE%Zc_pre, LAE%LS03pre, LAE%LLya, LAE%NLyC)
  END SUBROUTINE CalTelapsSFRandLS03
! -----------------------------------------------------------------------
  SUBROUTINE CalMultipleMerger(logic_dust)
    LOGICAL, INTENT(INOUT) :: logic_dust
    ! --- functions
    INTEGER :: DetMorType
    DOUBLE PRECISION :: CalNgas, Mstar_brst, CalTmass

    IF(paramLAE%type_b == 2 .and. LAE%tel < LAE%twind) THEN
       ! in the case of exponential burst and pre-outflow phase
       LAE%Mgas  = LAE%Mgas_pre * fdecline
       LAE%Mstar = Mstar_brst(LAE%Mgas_pre, LAE%Ms_pre, LAE%Mgas)
       IF(type /= 3) THEN ! starburst --> quiescent
          IF(gal(LAE%id)%Mcoold*param%munit > LAE%Mgas) THEN
             ! minor merger or cooling supplied sufficient cold gas
             type = 1
          ELSE
             type = 2
          ENDIF
          LAE%Mstar = LAE%Mstar + LAE%Ms_q
       ENDIF
    ELSE ! in the case of instantaneous burst OR in the case of exponential burst
         !                                        and outflow or post-outflow phases
       LAE%Mstar = LAE%Ms_pre + LAE%Ms_b; LAE%Mgas = 0.d0
       IF(type /= 3) THEN ! starburst --> quiescent
          type = 1; LAE%Mstar = LAE%Mstar + LAE%Ms_q
       ENDIF
    ENDIF
    LAE%Tmass = CalTmass(gal(LAE%id)%Tmass*param%th_yr*param%munit, LAE%tel, &
                         LAE%twind, LAE%taust_yr, LAE%Mgas_pre, &
                         ssp%alp(LAE%b_or_q)/LAE%beta) / LAE%Mstar
                ! mass-weighted age [yr]

    temp = (tlife - LAE%tel) / param%th_yr
           ! the onset time of star formation from 'tout' - 'tlife' [Gyr/h]
    IF(paramLAE%type_b == 1) THEN ! in the case of instantaneous burst
       call lum(param%zsp1, zp1pre, lumtmp2, mz, mt, temp)
    ELSEIF(paramLAE%type_b == 2) THEN ! in the case of exponential burst
       call lum_burst_LAE(param%zsp1, zp1pre, taust, lumtmp2, mz, mt, temp, tGWh, &
                          allgal(LAE%id)%Zc_pre)
    ENDIF
!!$    gal(LAE%id)%lumg(1:param%nwave) = allgal(LAE%id)%lumg_pre(1:param%nwave) &
!!$                                      + lumtmp2(1:param%nwave) * LAE%Mgas_pre
!!$    lumtmp(:) = gal(LAE%id)%lumg(:)

    LAE%lumg(1:param%nwave) = allgal(LAE%id)%lumg_pre(1:param%nwave) &
                              + lumtmp2(1:param%nwave) * LAE%Mgas_pre
    lumtmp(:) = LAE%lumg(:)


    IF(type == 3) THEN ! quiescent --> starburst
       morLAE = 1
    ELSE ! starburst --> quiescent
       morLAE = DetMorType(lumtmp(param%iBband_r), &
                           lumtmp(param%iBband_r+param%nwave))
    ENDIF

    IF(type /= 3) THEN ! starburst --> quiescent
       logic_dust = .true.
       IF(LAE%Mgas > 0.d0) THEN
          LAE%Ngas = CalNgas(LAE%Mgas, rb0_pre) ! [cm^-2]
          temp     = LAE%Zc_pre * const%Zsun ! metallicity
          LAE%NZ   = LAE%Ngas * temp ! [cm^-2]
          temp     = LAE%Mgas * temp + LAE%McZc_q
          LAE%Mgas = LAE%Mgas + LAE%Mgas_q
          LAE%Zc   = temp / (LAE%Mgas * const%Zsun) ! [Zsun]
       ELSE
          LAE%Ngas = 0.d0; LAE%NZ = 0.d0
          LAE%Zc   = allgal(LAE%id)%Zc / const%Zsun ! [Zsun]
          LAE%Mgas = gal(LAE%id)%Mcoold * param%munit
       ENDIF

       IF(morLAE == 1 .or. morLAE == 2) THEN
          LAE%reff = rb0
       ELSE
          LAE%reff = rd0
       ENDIF
       temp     = CalNgas(gal(LAE%id)%Mcoold*param%munit, LAE%reff)
       LAE%Ngas = LAE%Ngas + temp ! [cm^-2]
       LAE%NZ   = LAE%NZ   + temp * allgal(LAE%id)%Zc ! [cm^-2]
    ELSE ! quiescent --> starburst
       IF(LAE%Mgas > 0.d0) THEN
          logic_dust  = .true.; LAE%Ngas = CalNgas(LAE%Mgas, rb0_pre) ! [cm^-2]
          LAE%Zc = LAE%Zc_pre ! [Zsun]
          LAE%NZ = LAE%Ngas * LAE%Zc * const%Zsun ! [cm^-2]
       ELSE
          logic_dust  = .false.; LAE%Ngas = 0.d0; LAE%NZ = 0.d0; LAE%Zc = 0.d0
       ENDIF
    ENDIF

    IF(logic_dust) THEN
       call optNZ(param%nwave, lumtmp, lumtmp_d, LAE%NZ)
       call extcont(LAE%NZ, LAE%LS03, LAE%LS03d)
       mordLAE = DetMorType(lumtmp_d(param%iBband_r), &
                            lumtmp_d(param%iBband_r+param%nwave))
    ELSE
       lumtmp_d(:) = lumtmp(:); LAE%LS03d(:) = LAE%LS03(:); mordLAE = 1
    ENDIF
    call CalGalMag(logic_dust, lumtmp, lumtmp_d, LAE%mag, LAE%mag_d, param%nwave)
    DO iband = 1, param%nwave
       IF(LAE%flag_mag(iband) == 1) THEN
          ! LAE%mag(:), LAE%mag_d(:)
          ! = log10[NLyC/(photons/s/h^2) or log10[Llam/(erg/s/A/h^2)]
          LAE%mag(iband)   = -0.4d0 * LAE%mag(iband)
          LAE%mag_d(iband) = -0.4d0 * LAE%mag_d(iband)
        ENDIF
    ENDDO
  END SUBROUTINE CalMultipleMerger
! -----------------------------------------------------------------------
  SUBROUTINE CalStarburst(logic_dust)
    LOGICAL, INTENT(INOUT) :: logic_dust
    DOUBLE PRECISION :: CalNgas, Mstar_brst, CalTmass

    IF(paramLAE%type_b == 2 .and. LAE%tel < LAE%twind) THEN
       ! in the case of exponential burst and pre-outflow phase
       LAE%Mgas  = LAE%Mgas_pre * fdecline
       LAE%Mstar = Mstar_brst(LAE%Mgas_pre, LAE%Ms_pre, LAE%Mgas)
       LAE%Ngas  = CalNgas(LAE%Mgas, LAE%reff) ! [cm^-2]
    ELSE ! in the case of instantaneous burst
         ! OR in the case of exponential burst and outflow or post-outflow phases
       LAE%Mstar = LAE%Ms_pre + LAE%Ms_b
       LAE%Mgas  = 0.d0; LAE%Ngas = 0.d0
    ENDIF
    LAE%Tmass = CalTmass(gal(LAE%id)%Tmass*param%th_yr*param%munit, LAE%tel, &
                         LAE%twind, LAE%taust_yr, LAE%Mgas_pre, &
                         ssp%alp(LAE%b_or_q)/LAE%beta) / LAE%Mstar
                ! mass-weighted age [yr]

    LAE%Zc = LAE%Zc_pre ! [Zsun]
    LAE%NZ = LAE%Ngas * LAE%Zc_pre * const%Zsun ! [cm^-2]

    temp = (tlife - LAE%tel) / param%th_yr
           ! the onset time of star formation from 'tout-tlife' [Gyr/h]
    IF(paramLAE%type_b == 1) THEN ! in the case of instantaneous burst
       call lum(param%zsp1, zp1pre, lumtmp2, mz, mt, temp)
    ELSEIF(paramLAE%type_b == 2) THEN ! in the case of exponential burst
       call lum_burst_LAE(param%zsp1, zp1pre, taust, lumtmp2, mz, mt, temp, tGWh, &
                          allgal(LAE%id)%Zc_pre)
    ENDIF
!!$    gal(LAE%id)%lumg(1:param%nwave) = allgal(LAE%id)%lumg_pre(1:param%nwave) &
!!$                                      + lumtmp2(1:param%nwave) * LAE%Mgas_pre
!!$    lumtmp(:) = gal(LAE%id)%lumg(:)
    LAE%lumg(1:param%nwave) = allgal(LAE%id)%lumg_pre(1:param%nwave) &
                              + lumtmp2(1:param%nwave) * LAE%Mgas_pre
    lumtmp(:) = LAE%lumg(:)

    IF(LAE%Mgas > 0.d0) THEN
       logic_dust = .true.; call optNZ(param%nwave, lumtmp, lumtmp_d, LAE%NZ)
       call extcont(LAE%NZ, LAE%LS03, LAE%LS03d)
    ELSE
       logic_dust = .false.; lumtmp_d(:) = lumtmp(:); LAE%LS03d(:) = LAE%LS03(:)
    ENDIF
    call CalGalMag(logic_dust, lumtmp, lumtmp_d, LAE%mag, LAE%mag_d, param%nwave)
    DO iband = 1, param%nwave
       IF(LAE%flag_mag(iband) == 1) THEN
          ! LAE%mag(:), LAE%mag_d(:)
          ! = log10[NLyC/(photons/s/h^2) or log10[Llam/(erg/s/A/h^2)]
          LAE%mag(iband)   = -0.4d0 * LAE%mag(iband)
          LAE%mag_d(iband) = -0.4d0 * LAE%mag_d(iband)
       ENDIF
    ENDDO
  END SUBROUTINE CalStarburst
! -----------------------------------------------------------------------
  SUBROUTINE SubroutineSetForLyaBrightGalaxies
    call CalMagS03(logic_dust)
    call CalLyaLFDist(LAE%b_or_q, bin, mordLAE, LAE%weight, LAE%Ngas, &
                      LAE%Zc, LAE%NZ, LAE%Mstar)
    call TstarDist(LAE%b_or_q, mordLAE, LAE%taust_yr, LAE%weight)
    call WritePhysPropOfEachLAE(iLAE%numLya, iforest, morLAE, mordLAE)
  END SUBROUTINE SubroutineSetForLyaBrightGalaxies
! -----------------------------------------------------------------------
  SUBROUTINE SubroutineSetForUVBrightGalaxies
    call CalMagS03(logic_dust)
    call CalLFS03(LAE%sftype, morLAE, mordLAE, LAE%Mgas, LAE%weightAll)
    call WritePhysPropOfEachGal(iLAE%numUV, iforest, morLAE, mordLAE)
    call CalLF(LAE%b_or_q, morLAE,  LAE%mag,   param%nwave, LAE%weightAll, &
               LF_LAE(nz))
    call CalLF(LAE%b_or_q, mordLAE, LAE%mag_d, param%nwave, LAE%weightAll, &
               LFd_LAE(nz))
    call CalDists(LAE%b_or_q, mordLAE, nz, LAE%weightAll, LAE%SFR, LAE%reff, &
                  LAE%Ngas, LAE%Zc, LAE%NZ, LAE%Mstar, LAE%Mhost)
    ! --- for Makiya related (2011/01/31)
    ! call Output_Makiya(mordLAE, LAE%weightAll, weight)
  END SUBROUTINE SubroutineSetForUVBrightGalaxies
! -----------------------------------------------------------------------
END SUBROUTINE MainCalcForLAE
! ========================================================================
SUBROUTINE CalLF(b_or_q, mor, LAEmag, nwave, weight, LF_nz)
  use LAErelated; implicit none
  INTEGER, INTENT(IN) :: b_or_q, mor, nwave
  DOUBLE PRECISION, INTENT(IN) :: weight
  DOUBLE PRECISION, INTENT(IN) :: LAEmag(nwave)
  TYPE(DistributionFunction_LAE), INTENT(INOUT) :: LF_nz ! LF_LAE(nz) or LFd_LAE(nz)
  INTEGER :: iband, imag, CalBin

  DO iband = 1, nwave
     imag = CalBin(LAEmag(iband), 2.d0, -40.d0)
     IF(imag <= LF_nz%Nbin .and. imag >= 1) THEN
        LF_nz%n(b_or_q, mor, iband, imag) = &
             LF_nz%n(b_or_q, mor, iband, imag) + weight
     ENDIF
  ENDDO
END SUBROUTINE CalLF
! ========================================================================
DOUBLE PRECISION FUNCTION CalAbsMag(band, lum) RESULT(abs_mag)
  use LAErelated
  implicit none
  INTEGER, INTENT(IN) :: band
  DOUBLE PRECISION, INTENT(IN) :: lum ! [erg/s/h^2]
  DOUBLE PRECISION, PARAMETER :: EPS = 1.d-10

  abs_mag = 100.d0
  IF(lum > EPS) abs_mag = -2.5d0 * log10(lum) + corrS03(band)
END FUNCTION CalAbsMag
! ========================================================================
SUBROUTINE CalPreLS03(tel, metal, taueff, SFR)
  use LAErelated
  implicit none
  DOUBLE PRECISION, INTENT(IN) :: tel ! [yr]
  DOUBLE PRECISION, INTENT(IN) :: metal ! [Zsun]
  DOUBLE PRECISION, INTENT(IN) :: taueff ! [yr]
  DOUBLE PRECISION, INTENT(IN) :: SFR ! [Msun/yr]
  INTEGER :: i, istype, imax
  INTEGER :: CalIsType, CalImax ! functions
  DOUBLE PRECISION :: temp
  DOUBLE PRECISION :: LS03_constSF ! function

  istype = CalIsType(tel)
  imax   = CalImax(istype, tel)
  DO i = 1, NWAVE_S03
     LAE%LS03pre(i) = LS03_constSF(istype, i, tel, metal, imax)
  ENDDO

  temp = tel / taueff
  temp = SFR * temp * constLAE%convH2 / (1.d0 - exp(-temp))
  ! sfr(1,i) = mean SFR in the time step --> SFR at start of time step
  LAE%LS03pre(:) = LAE%LS03pre(:) * temp ![erg/s/h^2]
END SUBROUTINE CalPreLS03
! ========================================================================
DOUBLE PRECISION FUNCTION CalTbrstRatio(num, taueff1, taueff2, tstep, id, &
     flag_cent, MDM, weight) RESULT(tratio)
  implicit none
  INTEGER, INTENT(IN) :: num, id, flag_cent
  DOUBLE PRECISION, INTENT(IN) :: taueff1, taueff2, tstep
  DOUBLE PRECISION, INTENT(IN) :: MDM, weight
  CHARACTER*20 :: ctype(2) = (/CHARACTER*20 ::'satellite', 'central'/)

  tratio = taueff1 / tstep
  IF(num >= 2) THEN
     print '(A, I3, A, I10, 2(A, G11.3))', &
          '# (LAE)burst occured more than 2 times (', num, ' times) at i=', &
          id, ' ('//trim(ctype(flag_cent+1))//'), M_DM =', MDM, ', sum =', weight
     print '(2(A, 1X, E11.3), A)', '# (LAE)t_eff1 =', taueff1,  &
          '[yr] t_eff2 =', taueff2, '[yr]'
  ENDIF
END FUNCTION CalTbrstRatio
! ========================================================================
DOUBLE PRECISION FUNCTION CalLyaPhase(Rtbrst, tstep) RESULT(lphase)
  ! this term represents the time [yr] in which the Lya emission
  !  of a galaxy is bright enough to be detected
  ! if 'paramLAE%fphaseLya * t_eff^burst' is shorter than the lifetime
  !  of O stars (10^7 yr), 'lphase' is estimated by '1.d+7'
  ! the lifetime of O stars is taken from Maeder & Meynet
  !  (1989)
  use LAErelated
  implicit none
  DOUBLE PRECISION, INTENT(IN) :: Rtbrst, tstep

  IF(LAE%b_or_q == 2) THEN ! quiescent
     lphase = tstep
  ELSEIF(LAE%b_or_q == 1) THEN ! starburst
     lphase = max(2.d+7, paramLAE%fphaseLya * Rtbrst * tstep)
     IF(lphase > tstep) lphase = tstep
     ! in the case that the 'lphase' value of all bursting galaxies is set
     !  to be equal to 'tstep', there is no difference from the case that
     !  the 'lphase' value is estimated above method
     ! this is because the number weight of a sample galaxy is reduced
     !  by a factor of 'lphase' / 'tstep'
     ! so, this assumption in which 'paramLAE%fphaseLya' is equal to 4.0
     !  is adepuate that bursting galaxy can be detected via its Lya
     !  emission only in 'lphase' [yr] after the commencement of
     !  the bursting star formation
  ENDIF
END FUNCTION CalLyaPhase
!============================================================================
SUBROUTINE TstarDist(b_or_q, mor, Tstar, weight)
  use LAErelated
  implicit none
  INTEGER, INTENT(IN) :: b_or_q ! b_or_q = 1:starburst, 2:quiescent
  INTEGER, INTENT(IN) :: mor
  DOUBLE PRECISION, INTENT(IN) :: Tstar ! [yr]
  DOUBLE PRECISION, INTENT(IN) :: weight ! [h^3/Mpc^3]
  INTEGER :: bin, CalBin
  DOUBLE PRECISION, PARAMETER :: EPS = 1.d-10

  IF(Tstar > EPS) THEN
     bin = CalBin(log10(Tstar), distLAE(n_distLAE)%invstep, &
                  distLAE(n_distLAE)%base)
!!$         bin = int(aint(10.d0 * (log10(Tstar) - 1.d0)))
     IF(bin < 1) THEN
        bin = 1
     ELSEIF(bin > NLF1) THEN
        bin = NLF1
     ENDIF
  ELSE
     bin = 1
  ENDIF
  distLAE(n_distLAE)%count(b_or_q, mor, bin) &
       = distLAE(n_distLAE)%count(b_or_q, mor, bin) + weight
  distLAE(n_distLAE)%wcount(b_or_q, mor, bin) &
       = distLAE(n_distLAE)%wcount(b_or_q, mor, bin) + Tstar * weight
END SUBROUTINE TstarDist
! ========================================================================
SUBROUTINE CalLS03(sftype, tel, tend, taustinv, SFRini, metal, Lpre, LLya, NLyC)
  use LAErelated
  implicit none
  INTEGER, INTENT(IN) :: sftype ! 1:quiescent, 2:single major merger, 3:multiple-merger
  DOUBLE PRECISION, INTENT(IN) :: tel, tend, taustinv, SFRini
                                  ! tend and taustinv are used only for
                                  !  sftype = 2 or 3 and paramLAE%type_b = 2
  DOUBLE PRECISION, INTENT(IN) :: metal ! [Zsun] (used only for b_or_q = 1)
  DOUBLE PRECISION, INTENT(IN) :: Lpre(NWAVE_S03) ! used only for sftype = 3
  DOUBLE PRECISION, INTENT(INOUT) :: LLya, NLyC
  INTEGER :: i, istype, imax, imin
  INTEGER :: CalIsType, CalImax, CalImin ! functions
  DOUBLE PRECISION :: tmp
  DOUBLE PRECISION :: LS03_expSF, LS03_constSF ! functions

  IF(sftype == 2 .or. sftype == 3) THEN ! starburst
     IF(paramLAE%type_b == 1) THEN ! in the case of instantaneous burst
        istype = CalIsType(tel); imax = CalImax(istype, tel)
        tmp = SFRini * constLAE%convH2
        DO i = 1, NWAVE_S03
           LAE%LS03(i) = S03(istype)%lum(i,imax) * tmp ![erg/s/h^2]
        ENDDO
     ELSEIF(paramLAE%type_b == 2) THEN ! in the case of exponential burst
        istype = CalIsType(tel)
        imax   = CalImax(istype, tel); imin = CalImin(istype, tel, tend)
        IF(imax == imin) imax = imin + 1
        tmp = S03base(istype)%tstep * SFRini * constLAE%convH2
        DO i = 1, NWAVE_S03
           LAE%LS03(i) = LS03_expSF(i, istype, imax, imin, &
                                    taustinv, tel) * tmp ![erg/s/h^2]
        ENDDO
     ENDIF
     IF(sftype == 3) LAE%LS03(:) = LAE%LS03(:) + Lpre(:)
  ELSEIF(sftype == 1) THEN ! quiescent
     istype = CalIsType(tel); imax = CalImax(istype, tel)
     tmp = SFRini * constLAE%convH2
     DO i = 1, NWAVE_S03
        LAE%LS03(i) = LS03_constSF(istype, i, tel, metal, imax) * tmp ![erg/s/h^2]
     ENDDO
  ENDIF
  LLya = LAE%LS03(1); NLyC = LAE%LS03(1) / constLAE%fLya
  ! LLya [erg/s/h^2], NLyC [/s/h^2]
END SUBROUTINE CalLS03
! ========================================================================
INTEGER FUNCTION CalIsType(time) RESULT(istype)
  use LAErelated
  implicit none
  DOUBLE PRECISION, INTENT(IN) :: time

  IF(time > thLAE%time) THEN
     istype = 2
  ELSE
     istype = 1
  ENDIF
END FUNCTION CalIsType
! ========================================================================
INTEGER FUNCTION CalImax(istype, tel) RESULT(imax)
  use LAErelated
  implicit none
  INTEGER, INTENT(IN) :: istype
  DOUBLE PRECISION, INTENT(IN) :: tel
  INTEGER :: CalIS03

  IF(istype == 1 .or. istype == 2) THEN
     imax = min(NT_S03-1, CalIS03(istype, tel))
  ELSE
     print '(A, I5)', '# (LAE)Error in CalImax: called with '// &
          'undefined variable!', istype; stop
  ENDIF
END FUNCTION CalImax
! ========================================================================
INTEGER FUNCTION CalImin(istype, tel, tend) RESULT(imin)
  implicit none
  INTEGER, INTENT(IN) :: istype
  DOUBLE PRECISION, INTENT(IN) :: tel, tend
  INTEGER :: CalIS03

  IF(istype == 1 .or. istype == 2) THEN
     imin = max(1, CalIS03(istype, tel - tend))
  ELSE
     print '(A, I5)', '# (LAE)Error in CalImin: called with '// &
          'undefined variable!', istype; stop
  ENDIF
END FUNCTION CalImin
! ========================================================================
INTEGER FUNCTION CalIS03(istype, time) RESULT(iS03)
  use LAErelated
  implicit none
  INTEGER, INTENT(IN) :: istype
  DOUBLE PRECISION, INTENT(IN) :: time

  iS03 = int(aint(time * S03base(istype)%invtstep + S03base(istype)%base_iS03))
  ! anint: -4.d4 + 1.d5*i <= tel < 6.d4+1.d5*i ==> imax = i+1
  ! aint:   1.d4 + 1.d5*i = S03base(1)%time(i+1) <= tel
  !         < 1.1d5+1.d5*i = S03base(1)%time(i+2) ==> imax = i+1
  ! aint should be used: S03base(1)%time(imax) <= tel
  !                      < S03base(1)%time(imax+1)
END FUNCTION CalIS03
! ========================================================================
SUBROUTINE CalBaseVarUsingMetal(metal)
  use LAErelated
  implicit none
  DOUBLE PRECISION, INTENT(IN) :: metal ! [Zsun]
  INTEGER :: i, j

  DO i = 1, 2
     DO j = 1, NWAVE_S03
        call CalQH(i, j, metal)
     ENDDO
  ENDDO
CONTAINS
  SUBROUTINE CalQH(istype, num, metal)
    INTEGER, INTENT(IN) :: istype, num
    DOUBLE PRECISION, INTENT(IN) :: metal ! [Zsun]
    INTEGER :: i, iZ, iZm1, NZ
    DOUBLE PRECISION :: Z, x1, x2, y1, y2, tmp
    DOUBLE PRECISION, PARAMETER :: EPS = 1.d-9
    DOUBLE PRECISION :: LinearInterp

    Z  = metal * constLAE%Zsun ! metallicity (not in [Zsun])
    NZ = NZ_S03(istype)
    IF(Z > S03base(istype)%metal(NZ)) THEN
       S03(istype)%lum(num,:) = S03base(istype)%Lis(num,NZ,:)
    ELSE IF(Z < EPS) THEN
       S03(istype)%lum(num,:) = S03base(istype)%Lis(num, 1,:)
    ELSE
       i = 1
       DO WHILE(Z > S03base(istype)%metal(i))
          i = i + 1
       ENDDO
       iZ = i; iZm1 = i - 1
       ! S03base(istype)%metal(iZm1) < Z < S03base(istype)%metal(iZ)

       x1 = S03base(istype)%metal(iZ); x2 = S03base(istype)%metal(iZm1)
       tmp = 1.d0 / (x1 - x2)

       DO i = 1, NT_S03
          y1 = S03base(istype)%Lis(num,iZ,  i)
          y2 = S03base(istype)%Lis(num,iZm1,i)
          S03(istype)%lum(num, i) = LinearInterp(Z, x1, x2, y1, y2, tmp)
       ENDDO
    ENDIF
  END SUBROUTINE CalQH
END SUBROUTINE CalBaseVarUsingMetal
! ========================================================================
DOUBLE PRECISION FUNCTION CalElapsedTime(idum, tmax) RESULT(tel)
  ! --- this function returns the elapsed time in yr since
  !      the onset of starburst at the output time
  implicit none
  INTEGER, INTENT(INOUT) :: idum
  DOUBLE PRECISION, INTENT(IN) :: tmax
  REAL :: ran1

  tel = dble(ran1(idum)) * tmax
  DO WHILE(tel < 1.d+6)
     tel = dble(ran1(idum)) * tmax
  ENDDO
END FUNCTION CalElapsedTime
! ========================================================================
DOUBLE PRECISION FUNCTION LS03_constSF(istype, band, tel, metal_eff, iT) &
     RESULT(lum)
  use LAErelated
  implicit none
  INTEGER, INTENT(IN) :: istype, band, iT
  DOUBLE PRECISION, INTENT(IN) :: tel
  DOUBLE PRECISION, INTENT(IN) :: metal_eff ! [Zsun]
  INTEGER :: i, iTp1, iZ, iZm1, NZ
  DOUBLE PRECISION :: metal, x1, x2, y1, y2, base1, base2, tmp
  DOUBLE PRECISION, PARAMETER :: EPS = 1.d-9
  DOUBLE PRECISION :: LinearInterp

  metal = metal_eff * constLAE%Zsun ! metallicity (not in [Zsun])
  iTp1  = iT + 1; NZ = NZ_S03(istype)
  IF(metal > S03base(istype)%metal(NZ)) THEN
     x1 = S03base(istype)%time(iT);   y1 = S03base(istype)%Lcs(band,NZ,iT)
     x2 = S03base(istype)%time(iTp1); y2 = S03base(istype)%Lcs(band,NZ,iTp1)
     lum = LinearInterp(tel, x1, x2, y1, y2, 1.d0/(x1-x2))
  ELSE IF(metal < EPS) THEN
     x1 = S03base(istype)%time(iT);   y1 = S03base(istype)%Lcs(band,1,iT)
     x2 = S03base(istype)%time(iTp1); y2 = S03base(istype)%Lcs(band,1,iTp1)
     lum = LinearInterp(tel, x1, x2, y1, y2, 1.d0/(x1-x2))
  ELSE
     i = 1
     DO WHILE(metal > S03base(istype)%metal(i))
        i = i + 1
     ENDDO
     iZ = i; iZm1 = i - 1
     ! S03base(istype)%metal(iZm1) < metal <= S03base(istype)%metal(iZ)

     x1 = S03base(istype)%metal(iZ); x2 = S03base(istype)%metal(iZm1)
     tmp = 1.d0 / (x1 - x2)

     y1 = S03base(istype)%Lcs(band, iZ,   iT)
     y2 = S03base(istype)%Lcs(band, iZm1, iT)
     base1 = LinearInterp(metal, x1, x2, y1, y2, tmp)

     y1 = S03base(istype)%Lcs(band, iZ,   iTp1)
     y2 = S03base(istype)%Lcs(band, iZm1, iTp1)
     base2 = LinearInterp(metal, x1, x2, y1, y2, tmp)

     x1  = S03base(istype)%time(iT); x2 = S03base(istype)%time(iTp1)
     lum = LinearInterp(tel, x1, x2, base1, base2, 1.d0/(x1-x2))
  ENDIF
END FUNCTION LS03_constSF
! ========================================================================
DOUBLE PRECISION FUNCTION LS03_expSF(band, istype, imax, imin, taustinv, &
                                     tstart) RESULT(lum)
  use LAErelated
  implicit none
  INTEGER, INTENT(IN) :: band, istype, imax, imin
  DOUBLE PRECISION, INTENT(IN) :: taustinv, tstart
  INTEGER :: i, imaxp1
  DOUBLE PRECISION :: x1, x2, y1, y2, tmp
  DOUBLE PRECISION :: trapzd, LinearInterp

  lum = 0.d0
  DO i = imin, imax
     IF(i == imin .or. i == imax) THEN
        trapzd = 0.5d0
     ELSE
        trapzd = 1.d0
     ENDIF
     tmp = (tstart - S03base(istype)%time(i)) * taustinv
     lum = lum + trapzd * exp(-tmp) * S03(istype)%lum(band, i)
  ENDDO
  imaxp1 = imax + 1
  x1 = S03base(istype)%time(imax); x2 = S03base(istype)%time(imaxp1)
  y1 = S03(istype)%lum(band,imax); y2 = S03(istype)%lum(band,imaxp1)
  lum = lum + LinearInterp(tstart, x1, x2, y1, y2, 1.d0/(x1-x2))
END FUNCTION LS03_expSF
! ========================================================================
DOUBLE PRECISION FUNCTION CalTmass(Tmass0, tel, twind, taust, Mc0, abratio) &
     RESULT(Tmass) ! mass-weighted age[yr]
  DOUBLE PRECISION, INTENT(IN) :: Tmass0, tel, twind, taust, Mc0, abratio
  DOUBLE PRECISION :: x, xwind

  IF(tel >= twind) THEN
     x = tel / taust; xwind = twind / taust
     Tmass = Mc0 * taust * exp(-x) * (abratio * (x - xwind + 1.e0) - xwind)
  ELSE
     x = tel / taust
     Tmass = Mc0 * taust * (1.e0 - (x + 1.e0) * exp(-x))
  ENDIF
  Tmass = Tmass + Tmass0
END FUNCTION CalTmass
! ========================================================================
SUBROUTINE extcont(NZ, LS03, LS03d)
  use LAErelated
  implicit none
  DOUBLE PRECISION, INTENT(IN) :: NZ, LS03(NWAVE_S03)
  DOUBLE PRECISION, INTENT(INOUT) :: LS03d(NWAVE_S03)
  INTEGER :: iband
  DOUBLE PRECISION :: tauV, taud, fesc
  DOUBLE PRECISION :: f_scr, f_slab ! function

  tauV = constLAE%OptV * NZ
  DO iband = 2, NWAVE_S03
     taud = tauV * xi_S03(iband)
     IF(paramLAE%exttype == 1) THEN ! screen dust geometry
        fesc = f_scr(taud)
     ELSEIF(paramLAE%exttype == 2) THEN ! slab dust geometry
        fesc = f_slab(taud)
     ENDIF
     LS03d(iband) = fesc * LS03(iband)
  ENDDO
END SUBROUTINE extcont
! ========================================================================
SUBROUTINE CalGalMag(logic_dust, lum, lumd, LAEmag, LAEmag_d, nwave)
  implicit none
  LOGICAL, INTENT(IN) :: logic_dust ! true: w/ dust, false: w/o dust
  INTEGER, INTENT(IN) :: nwave
  DOUBLE PRECISION, INTENT(IN) :: lum(2*nwave), lumd(2*nwave)
  DOUBLE PRECISION, INTENT(INOUT) :: LAEmag(nwave), LAEmag_d(nwave)
  INTEGER :: iband, iband2
  DOUBLE PRECISION :: CalGalAllMag ! function

  IF(logic_dust) THEN
     DO iband = 1, nwave
        iband2 = nwave + iband
        LAEmag(iband)   = CalGalAllMag(lum(iband),  lum(iband2))
        LAEmag_d(iband) = CalGalAllMag(lumd(iband), lumd(iband2))
     ENDDO
  ELSE
     DO iband = 1, nwave
        iband2 = nwave + iband
        LAEmag(iband) = CalGalAllMag(lum(iband), lum(iband2))
     ENDDO
     LAEmag_d(:) = LAEmag(:)
  ENDIF
END SUBROUTINE CalGalMag
! ========================================================================
DOUBLE PRECISION FUNCTION CalGalAllMag(lumb, lumd) RESULT(mag_gal)
  use global_var
  implicit none
  DOUBLE PRECISION, INTENT(IN) :: lumb, lumd

  mag_gal = const%corr -2.5d0 * log10(lumb + lumd)
END FUNCTION CalGalAllMag
! ========================================================================
DOUBLE PRECISION FUNCTION CalNgas(Mcold, r_eff)
  ! --- calculate the column density 'column_denst' of the galaxy
  !      according to its mass of the cold gas 'Mcold' [M_sun]
  !      and the effective radius 'r_eff' [kpc/h]
  use LAErelated
  implicit none
  DOUBLE PRECISION, INTENT(IN) :: Mcold, r_eff
  DOUBLE PRECISION :: Square

  IF(Mcold > 0.d0) THEN
     ! --- mean column density N_0 [cm^-2] is estimated
     CalNgas = constLAE%Ngas * Mcold / Square(r_eff)
  ELSE
     CalNgas = 0.d0
  ENDIF
END FUNCTION CalNgas
! ========================================================================
DOUBLE PRECISION FUNCTION Mstar_brst(Mcold0, Mstar0, Mcold)
  ! --- The stellar mass of a burst galaxy at t after a burst
  ! --- d(Mstar)/dt = SFR(t) = Mc(t) / t_brst
  !                 = Mc(t=0) * exp[-t/t_brst] / t_brst
  !      ==> Mstar(t) = Mstar(t=0) + Mc(t=0) * [1-exp(-t/t_brst)]
  ! --- The relationship between 'Mcold0' and 'Mcold' is the following:
  !       Mcold = Mcold0 * exp[-t/t_brst]
  ! --- 'Mcold0' [Msun], 'Mstar0' [Msun], 'Mcold' [Msun], 'Mstar_brst' [Msun]
  DOUBLE PRECISION, INTENT(IN) :: Mcold0, Mstar0, Mcold

  Mstar_brst = Mstar0 + (Mcold0 - Mcold) ! [Msun]
END FUNCTION Mstar_brst
! ========================================================================
DOUBLE PRECISION FUNCTION Zc_brst(Zcold0, t, t_eff, b, b_or_q)
  ! --- The metallicity in cold gas phase of a burst galaxy remaind at t
  !      after a burst
  ! --- The metal mass in cold phase increases according to the chemical yield
  !      by the star formation whose typical timescale is t_brst and decreases
  !      by the reheating by SN feedback mechanism
  ! --- d(McZc)/dt = [p-(alpha+beta)]*SFR(t)
  !                = p*Mc(t=0)*exp(-(alpha+beta)*t/t_brst)/t_brst
  !                   - ((alpha+beta)/t_brst)*Mc(t)*Zc(t)
  !                = p*[Mc(t=0)/t_brst]*exp(-t/t_eff) - Mc(t)*Zc(t)/t_eff
  !      ==> Mc(t)Zc(t) = Mc(t=0) * (Zc(t=0) + p * (t/t_brst)) * exp(-t/t_eff)
  !      ==> Zc(t) = Zc(t=0) + p * (t/t_brst)
  ! --- p = y * alpha
  ! --- 'Zcold0' [Zsun]
  use global_var
  implicit none
  INTEGER, INTENT(IN) :: b_or_q ! 1:starubrst, 2:quiescent
  DOUBLE PRECISION, INTENT(IN) :: Zcold0, t, t_eff, b

  Zc_brst = min(1.d0, Zcold0 * const%Zsun &
                      + ssp%p(b_or_q) * t / (t_eff * (ssp%alp(b_or_q) + b)))
END FUNCTION Zc_brst
! ========================================================================
SUBROUTINE CalLyaLFDist(b_or_q, bin, mor, weight, Ngas, Zc, NgasZc, Mstar)
  use LAErelated
  implicit none
  INTEGER, INTENT(IN) :: b_or_q, bin, mor
  DOUBLE PRECISION, INTENT(IN) :: weight, Ngas, Zc, NgasZc, Mstar
  INTEGER :: i
  DOUBLE PRECISION :: x

  call CalLyaLF
  DO i = 1, n_distLya
     IF(i == 1) THEN
        x = Ngas
     ELSEIF(i == 2) THEN
        x = Zc
     ELSEIF(i == 3) THEN
        x = Mstar
     ELSEIF(i == 4) THEN
        x = NgasZc
     ENDIF
     call CalDistLya(i, bin, x, distLya(i)%invstep, distLya(i)%base)
  ENDDO
CONTAINS
! -----------------------------------------------------------------------
  SUBROUTINE CalLyaLF
    LyaLF%count(b_or_q, mor, bin) = LyaLF%count(b_or_q, mor, bin) + weight
  END SUBROUTINE CalLyaLF
! -----------------------------------------------------------------------
  SUBROUTINE CalDistLya(num, binLya, x, invstep, base)
    INTEGER, INTENT(IN) :: num, binLya
    DOUBLE PRECISION, INTENT(IN) :: x, invstep, base
    INTEGER :: bin, CalBin
    DOUBLE PRECISION, PARAMETER :: EPS = 1.d-10

    IF(x > EPS) THEN
       bin = CalBin(log10(x), invstep, base)
       IF(bin <= NLF2 .and. bin >= 1) &
            distLya(num)%count(b_or_q, bin, binLya) &
            = distLya(num)%count(b_or_q, bin, binLya) + weight
    ENDIF
  END SUBROUTINE CalDistLya
! -----------------------------------------------------------------------
END SUBROUTINE CalLyaLFDist
! ========================================================================
SUBROUTINE CalMagS03(logic_dust)
  use LAErelated
  implicit none
  LOGICAL, INTENT(IN) :: logic_dust ! true: w/ dust, false: w/o dust
  INTEGER :: iband
  DOUBLE PRECISION :: CalAbsMag ! function

  IF(logic_dust) THEN ! w/ dust
     DO iband = 2, NWAVE_S03
        LAE%magS03(iband)  = CalAbsMag(iband, LAE%LS03(iband))
        LAE%magS03d(iband) = CalAbsMag(iband, LAE%LS03d(iband))
     ENDDO
  ELSE
     DO iband = 2, NWAVE_S03
        LAE%magS03(iband) = CalAbsMag(iband, LAE%LS03(iband))
     ENDDO
     LAE%magS03d(:) = LAE%magS03(:)
  ENDIF
END SUBROUTINE CalMagS03
! ========================================================================
SUBROUTINE CalDists(b_or_q, mor, nz, weight, SFR, r_eff, Ngas, Zc, &
                    NgasZc, Mstar, Mhost)
  use LAErelated
  implicit none
  INTEGER, INTENT(IN) :: b_or_q, mor, nz
  DOUBLE PRECISION, INTENT(IN) :: weight, SFR, r_eff, Ngas, Zc, NgasZc
  DOUBLE PRECISION, INTENT(IN) :: Mstar, Mhost
  INTEGER :: i, iend
  DOUBLE PRECISION :: x

  iend = n_distLAE - 1
  DO i = 1, iend
     IF(i == 1) THEN
        x = Ngas ! [cm^-2]
     ELSEIF(i == 2) THEN
        x = Zc ! [Zsun]
     ELSEIF(i == 3) THEN
        x = NgasZc / constLAE%Zsun ! [Zsun/cm^2]
     ELSEIF(i == 4) THEN
        x = Mstar ! [Msun]
     ELSEIF(i == 5) THEN
        x = Mhost ! [Msun]
     ELSEIF(i == 6) THEN
        x = r_eff ! [kpc/h]
     ELSEIF(i == 7) THEN
        x = SFR ! [Msun/yr]
     ENDIF
     call CalDistLAE(i, x, distLAE(i)%invstep, distLAE(i)%base)
  ENDDO
  CSFRD(nz, b_or_q, mor) = CSFRD(nz, b_or_q, mor) + SFR * weight
CONTAINS
! -----------------------------------------------------------------------
  SUBROUTINE CalDistLAE(num, x, invstep, base)
    INTEGER, INTENT(IN) :: num
    DOUBLE PRECISION, INTENT(IN) :: x, invstep, base
    INTEGER :: bin, CalBin
    DOUBLE PRECISION, PARAMETER :: EPS = 1.d-10
    CHARACTER(LEN=100), PARAMETER :: name(n_distLAE) &
         = (/CHARACTER(LEN=100):: 'Ngas', 'Zc', 'NcZc', 'Mstar', 'Mhost', &
                                  'reff', 'SFR', 'taueff'/)

    IF(x > EPS) THEN
       bin = CalBin(log10(x), invstep, base)
       IF(bin < 1) THEN
          bin = 1
       ELSEIF(bin > NLF1) THEN
          bin = NLF1
       ENDIF
    ELSE
       bin = 1
       IF(num /= 1 .and. num /= 2 .and. num /= 3 .and. num /= 7) &
            print '(A, G10.2)', '# (LAE)CalDistLAE: x < EPS for '//&
            trim(name(num))//': x=', x
    ENDIF
    IF(bin <= NLF1 .and. bin >= 1) THEN
       distLAE(num)%count(b_or_q, mor,bin) &
            = distLAE(num)%count(b_or_q, mor, bin) + weight
       distLAE(num)%wcount(b_or_q,mor,bin) &
            = distLAE(num)%wcount(b_or_q, mor, bin) + x * weight
    ENDIF
  END SUBROUTINE CalDistLAE
! -----------------------------------------------------------------------
END SUBROUTINE CalDists
! ========================================================================
DOUBLE PRECISION FUNCTION CalFdecline(sb_type, time, taust) RESULT(fdecline)
  implicit none
  INTEGER, INTENT(IN) :: sb_type ! 1:instantaneous, 2:exp-decaying
  DOUBLE PRECISION, INTENT(IN) :: time, taust ! [yr]

  IF(sb_type == 1) THEN
     fdecline = 0.d0
  ELSEIF(sb_type == 2) THEN
     fdecline = - time / taust; fdecline = exp(fdecline)
  ELSE
     print '(A, I5)', &
          '# (LAE)CalFdecline is called with wrong argument!:', sb_type; stop
  ENDIF
END FUNCTION CalFdecline
! ========================================================================
DOUBLE PRECISION FUNCTION CalSFRburst(time, tend, SFRini, fdec) RESULT(SFRb)
  implicit none
  DOUBLE PRECISION, INTENT(IN) :: time, tend ! [yr]
  DOUBLE PRECISION, INTENT(IN) :: SFRini ! [Msun/yr]
  DOUBLE PRECISION, INTENT(IN) :: fdec

  SFRb = 0.d0
  IF(time < tend) SFRb = SFRini * fdec ! [Msun/yr]
END FUNCTION CalSFRburst
! ========================================================================
SUBROUTINE CalLFS03(sftype, mor, mord, Mgas, weight)
  use LAErelated
  implicit none
  INTEGER, INTENT(IN) :: sftype, mor, mord
  DOUBLE PRECISION, INTENT(IN) :: Mgas, weight
  INTEGER :: i, bin, CalBin

  DO i = 2, NWAVE_S03
     bin = CalBin(LAE%magS03(i), 5.d0, -30.d0)
     IF(bin <= NLF1 .and. bin >= 1) THEN
        LF_S03(i)%count(sftype, mor, bin) &
             = LF_S03(i)%count(sftype, mor, bin) + weight
     ENDIF

     IF(sftype == 1 .or. Mgas > 0.d0) THEN
        bin = CalBin(LAE%magS03d(i), 5.d0, -30.d0)
     ENDIF
     IF(bin <= NLF1 .and. bin >= 1 .and. mord /= 0) THEN
        LF_S03(i)%countd(sftype, mord, bin) &
             = LF_S03(i)%countd(sftype, mord, bin) + weight
     ENDIF
  ENDDO
END SUBROUTINE CalLFS03
! ========================================================================
SUBROUTINE WritePhysPropOfEachLAE(num, iforest, morLAE, mordLAE)
  ! --- write physical properties of the bursting galaxies whose L(Lya) are
  !      larger than thLAE%LLya [erg/s/h^2])
  use global_var;  use LAErelated
  implicit none
  INTEGER, INTENT(IN) :: iforest, morLAE, mordLAE
  INTEGER, INTENT(INOUT) :: num
  INTEGER :: iband, iout
  DOUBLE PRECISION :: CalMeanSFR, EBV, CalEB_V

  IF(LAE%LLya >= thLAE%LLya .and. LAE%Mstar >= thLAE%Mstar) THEN
     EBV = CalEB_V(LAE%NZ); num = num + 1
     call MagCheck(LAE%mag, LAE%mag_d); LAE%mSFR = LAE%SFR
     IF(LAE%b_or_q == 1) THEN ! starburst
        LAE%mSFR = CalMeanSFR(LAE%tel, LAE%twind, LAE%taust_yr, &
                              LAE%Mgas_pre, LAE%Ms_b, LAE%beta, LAE%ab)
        IF(LAE%sftype == 3 .and. gal(LAE%id)%Mcoold > 0.d0) & ! multiple merger
             LAE%mSFR = LAE%mSFR + allgal(LAE%id)%SFR(2)
     ENDIF

     iout = iLAE%eachLAE
     IF(param%run_type == 1) iout = paramLAE%base + paramLAE%nfile + inode + 1

     ! --- integers
     write(iout, '(I8, 4I2, $)') LAE%id, morLAE, mordLAE, &
          LAE%b_or_q, LAE%sftype
     write(iout, '(I2, I8, 4I16, $)') LAE%flag_c, iforest, &
          gal(LAE%id)%mpi, mrgt%mpi(gal(LAE%id)%IDhost), gal(LAE%id)%IDhost, &
          gal(LAE%id)%id_cgal
     ! --- double precisions
     write(iout, '(10X, 2(10E13.4, 10X), 3E13.4, $)') LAE%NZ, &
          LAE%weight, LAE%LLya, LAE%Mstar, LAE%reff, LAE%Zc, LAE%Mhost, &
          LAE%beta, LAE%taust_yr, LAE%Mgas, LAE%Ms_pre, LAE%Mprog, LAE%Vc, &
          LAE%tel, LAE%Vcc, LAE%Vel, LAE%SFR, LAE%twind, EBV, LAE%NLyC, &
          LAE%Mgas_pre, LAE%mSFR, LAE%Tmass
     DO iband = 2, NWAVE_S03
        write(iout, '(2(1X, F9.4), 2X, $)') &
             LAE%magS03(iband), LAE%magS03d(iband)
     ENDDO
!!$     mag_q(1) = 100.d0
!!$    IF(gal(LAE%id)%lumq(1) > 0.d0) &
!!$         mag_q(1) = - 2.5d0 * log10(gal(LAE%id)%lumq(1) / const%nub) &
!!$           - 65.35d0 + const%corr
!!$    write(iout, '(10X, F9.4, 2E13.4, $)') mag_q(1), &
!!$         gal(LAE%id)%lumq(2), gal(LAE%id)%lumq(3)
     write(iout, '(10X, 2(1X, F9.4), 2X, $)') LAE%mag(1), LAE%mag_d(1)
     DO iband = 2, param%nwave
        write(iout, '(2(1X, F9.4), 2X, $)') &
             LAE%mag(iband), LAE%mag_d(iband)
     ENDDO
     write(iout, *)
  ENDIF
CONTAINS
! -----------------------------------------------------------------------
  SUBROUTINE MagCheck(LAEmag, LAEmag_d)
    DOUBLE PRECISION, INTENT(INOUT) :: LAEmag(param%nwave), LAEmag_d(param%nwave)
    INTEGER :: i

    DO i = 1, param%nwave
       IF(LAEmag(i)   > 200.d0) LAEmag(i)   = 200.d0
       IF(LAEmag_d(i) > 200.d0) LAEmag_d(i) = 200.d0
    ENDDO
  END SUBROUTINE MagCheck
! -----------------------------------------------------------------------
END SUBROUTINE WritePhysPropOfEachLAE
! ========================================================================
SUBROUTINE WritePhysPropOfEachGal(num, iforest, morLAE, mordLAE)
  use global_var; use LAErelated
  implicit none
  INTEGER, INTENT(IN) :: iforest, morLAE, mordLAE
  INTEGER, INTENT(INOUT) :: num
  INTEGER :: iband, iout
  DOUBLE PRECISION :: CalMeanSFR, EBV, CalEB_V ! functions

!!$  IF(LAE%magS03(2) <= thLAE%magUV .or. LAE%flag_c == 1) THEN
  IF(LAE%magS03(2) <= thLAE%magUV) THEN
     EBV = CalEB_V(LAE%NZ); num = num + 1; LAE%mSFR = LAE%SFR
     IF(LAE%b_or_q == 1) THEN ! starburst
        LAE%mSFR = CalMeanSFR(LAE%tel, LAE%twind, LAE%taust_yr, &
                              LAE%Mgas_pre, LAE%Ms_b, LAE%beta, LAE%ab)
        IF(LAE%sftype == 3 .and. gal(LAE%id)%Mcoold > 0.d0) & ! multiple merger
             LAE%mSFR = LAE%mSFR + allgal(LAE%id)%SFR(2)
     ENDIF

     iout = iLAE%each
     IF(param%run_type == 1) &
          iout = paramLAE%base + paramLAE%nfile + nnode + inode + 1

     ! --- integers
     write(iout, '(I8, 4I2, $)') LAE%id, morLAE, mordLAE, &
          LAE%b_or_q, LAE%sftype
     write(iout, '(I2, I8, 4I16, $)') LAE%flag_c, iforest, gal(LAE%id)%mpi, &
          mrgt%mpi(gal(LAE%id)%IDhost), gal(LAE%id)%IDhost, gal(LAE%id)%id_cgal
     ! --- double precisions
     write(iout, '(10X, 2(10E13.4, 10X), 3E13.4, 10X, $)') LAE%NZ, &
          LAE%weightAll, LAE%LLya, LAE%Mstar, LAE%reff, LAE%Zc, LAE%Mhost, &
          LAE%beta, LAE%taust_yr, LAE%Mgas, LAE%Ms_pre, LAE%Mprog, LAE%Vc, &
          LAE%tel, LAE%Vcc, LAE%Vel, LAE%SFR, LAE%twind, EBV, LAE%NLyC, &
          LAE%Mgas_pre, LAE%mSFR, LAE%Tmass
     DO iband = 2, NWAVE_S03
        write(iout, '(2(1X, F9.4), 2X, $)') &
             LAE%magS03(iband), LAE%magS03d(iband)
     ENDDO
!!$     mag_q(1) = 100.d0
!!$     IF(gal(LAE%id)%lumq(1) > 0.d0) &
!!$          mag_q(1) = - 2.5d0 * log10(gal(LAE%id)%lumq(1) / const%nub) &
!!$          -65.35d0 + const%corr
!!$     write(iout, '(10X, F9.4, 2E13.4, 10X, $)') mag_q(1), &
!!$          gal(LAE%id)%lumq(2), gal(LAE%id)%lumq(3)
     write(iout, '(10X, 2(1X, F9.4), 2X, $)') LAE%mag(1), LAE%mag_d(1)
     DO iband = 2, param%nwave
        write(iout, '(2(1X, F9.4), 2X, $)') &
             LAE%mag(iband), LAE%mag_d(iband)
     ENDDO
     write(iout, *)
  ENDIF
END SUBROUTINE WritePhysPropOfEachGal
! ========================================================================
DOUBLE PRECISION FUNCTION CalMeanSFR(tel, twind, taust, Mc0, Msb, beta, ab) &
     RESULT(mSFR)
  ! --- mean SFR past t0 [yr] for starburst galaxies
  DOUBLE PRECISION, INTENT(IN) :: tel, twind, taust ! [yr]
  DOUBLE PRECISION, INTENT(IN) :: Mc0, Msb ! [Msun]
  DOUBLE PRECISION, INTENT(IN) :: beta, ab
  DOUBLE PRECISION, PARAMETER :: t0 = 1.d+7 ! [yr]
  DOUBLE PRECISION :: tmp1, tmp2

  IF(tel - twind > t0) THEN
     mSFR = 0.d0
  ELSE
     IF(tel < t0) THEN
        IF(tel < twind) THEN
           tmp1 = tel / taust
           mSFR = Mc0 * (1.d0 - exp(-tmp1)) / tel
        ELSE
           mSFR = Msb / tel
        ENDIF
     ELSE
        IF(tel < twind) THEN
           tmp1 = (tel - t0) / taust; tmp2 = tel / taust
           mSFR = Mc0 * (exp(-tmp1) - exp(-tmp2)) / t0
        ELSE
           tmp1 = (tel - t0) / taust
           mSFR = Mc0 * (exp(-tmp1) - beta / ab) / t0
        ENDIF
     ENDIF
  ENDIF
  IF(mSFR < 0.d0) THEN
     print '(A)', '# (LAE)Mean SFR is negative!!'
     print '(3(A, G11.4))', '  --- tel = ', tel, ', twind = ', twind, &
          ', taust = ', taust
     print '(4(A, G11.4))', '  --- Mc0 = ', Mc0, ', Msb = ', Msb, &
          ', beta = ', beta, ', ab = ', ab
  ENDIF
END FUNCTION CalMeanSFR
! ========================================================================
SUBROUTINE WriteLyaLF
  use LAErelated
  implicit none
  INTEGER :: i, j, k
  DOUBLE PRECISION, PARAMETER :: step = 0.1d0
  DOUBLE PRECISION, PARAMETER :: base = 30.d0
  DOUBLE PRECISION :: x, invstep, cumn(2,3), tmp(2,3)
  DOUBLE PRECISION, PARAMETER :: EPS = 1.d-15

  invstep = 1.d0 / step
  cumn(:,:) = 0.d0; tmp(:,:) = 0.d0
  DO i = NLF1, 1, -1
     x = dble(i) * step + base; cumn(:,:) = cumn(:,:) + LyaLF%count(:,:,i)
     DO k = 1, 3
        IF(LyaLF%count(1,k,i) < EPS) LyaLF%count(1,k,i) = EPS
        IF(LyaLF%count(2,k,i) < EPS) LyaLF%count(2,k,i) = EPS
        DO j = 1, 2
           IF(cumn(j,k) < EPS) THEN
              tmp(j,k) = EPS
           ELSE
              tmp(j,k) = cumn(j,k)
           ENDIF
        ENDDO
     ENDDO
     IF(x <= 45.d0 .and. x >= 35.d0) THEN
        write(iLAE%LyaLF, '(13(F9.5, 1X))')  x, &
             (log10(invstep * LyaLF%count(2,k,i)), k=1,3), & ! quiescent
             (log10(invstep * LyaLF%count(1,k,i)), k=1,3), & ! starburst
             (log10(tmp(2,k)), k=1,3), & ! quiescent
             (log10(tmp(1,k)), k=1,3)    ! burst
     ENDIF
  ENDDO
END SUBROUTINE WriteLyaLF
! ========================================================================
SUBROUTINE WriteS03LFs
  use LAErelated
  implicit none
  INTEGER :: i, j, k
  CHARACTER(LEN=10) :: lam(NWAVE_S03) &
       = (/CHARACTER(LEN=10):: '***','1500','2800','6563','4861','1216'/)
  DOUBLE PRECISION, PARAMETER :: step = 0.2d0
  DOUBLE PRECISION, PARAMETER :: base = -30.d0
  DOUBLE PRECISION :: x, invstep

  invstep = 1.d0 / step ! for mag^-1
  DO j = 2, NWAVE_S03
     write(iLAE%contLF, '(A, I2, A)') '# index ', j-2, ': '//lam(j)//'A LF'
     DO i = 1, NLF1
        x = dble(i) * step + base ! -29.8~10.0
        write(iLAE%contLF, '(F9.5, 18(2X, E13.6))') x,&
             (invstep * LF_S03(j)%count(1, k,i), k=1,3),&
             (invstep * LF_S03(j)%count(2, k,i), k=1,3),&
             (invstep * LF_S03(j)%count(3, k,i), k=1,3),&
             (invstep * LF_S03(j)%countd(1,k,i), k=1,3),&
             (invstep * LF_S03(j)%countd(2,k,i), k=1,3),&
             (invstep * LF_S03(j)%countd(3,k,i), k=1,3)
     ENDDO
     write(iLAE%contLF, *); write(iLAE%contLF, *)
  ENDDO
END SUBROUTINE WriteS03LFs
! ========================================================================
SUBROUTINE WriteDistributionsForLAE
  use LAErelated
  implicit none
  INTEGER :: i, j, k
  DOUBLE PRECISION :: x, cumn(2,3)

  DO i = 1, n_distLAE
     cumn(:,:) = 0.d0
     DO j = NLF1, 1, -1
        x = dble(j) * distLAE(i)%step + distLAE(i)%base
        cumn(:,:) = cumn(:,:) + distLAE(i)%count(:,:,j)
        write(iLAE%dist(i), '(19(G12.4, 1X))') x, &
             (distLAE(i)%invstep * distLAE(i)%wcount(2,k,j), k=1,3), & ! quiescent
             (distLAE(i)%invstep * distLAE(i)%wcount(1,k,j), k=1,3), & ! starburst
             (distLAE(i)%invstep * distLAE(i)%count(2, k,j), k=1,3), & ! quiescent
             (distLAE(i)%invstep * distLAE(i)%count(1, k,j), k=1,3), & ! starburst
             (cumn(2,k), k=1,3), & ! quiescent
             (cumn(1,k), k=1,3)    ! starburst
     ENDDO
  ENDDO
END SUBROUTINE WriteDistributionsForLAE
! ========================================================================
SUBROUTINE Write2DxLyaDistribution
  use LAErelated
  implicit none
  INTEGER :: i, ii, j, k
  DOUBLE PRECISION, PARAMETER :: step = 0.1d0
  DOUBLE PRECISION, PARAMETER :: base = 30.d0
  DOUBLE PRECISION :: x(NLF1), invstep, cumn(2,NLF2), tmp(2,NLF2)
  DOUBLE PRECISION, PARAMETER :: EPS = 1.d-15

  DO i = NLF1, 1, -1
     x(i) = dble(i) * step + base ! log10[L(Lya)/(erg/s/h^2)]
  ENDDO

  invstep = 1.d0 / step
  DO ii = 1, n_distLya
     cumn(:,:) = 0.d0; tmp(:,:) = 0.d0
     DO i = NLF1, 1, -1
        cumn(:,:) = cumn(:,:) + distLya(ii)%count(:,:,i)
        DO k = 1, NLF2
           IF(distLya(ii)%count(1,k,i) < EPS) distLya(ii)%count(1,k,i) = EPS
           IF(distLya(ii)%count(2,k,i) < EPS) distLya(ii)%count(2,k,i) = EPS
           DO j = 1, 2
              IF(cumn(j,k) < EPS) THEN
                 tmp(j,k) = 1.d-10
              ELSE
                 tmp(j,k) = cumn(j,k)
              ENDIF
           ENDDO
        ENDDO
        IF(x(i) <= 45.d0 .and. x(i) >= 35.d0) THEN
           write(iLAE%distLya(ii), '(41(F9.5, 1X))') x(i),&
                (log10(invstep * distLya(ii)%count(1,k,i)), k=1,NLF2),&
                (log10(invstep * distLya(ii)%count(2,k,i)), k=1,NLF2),&
                (log10(tmp(1,k)), k=1,NLF2), (log10(tmp(2,k)), k=1,NLF2)
        ENDIF
     ENDDO
  ENDDO
END SUBROUTINE Write2DxLyaDistribution
! ========================================================================
SUBROUTINE WriteCaptionsForLAE(ilog, run_redshift)
  use global_var; use LAErelated
  implicit none
  INTEGER, INTENT(IN) :: ilog, run_redshift
  INTEGER :: i, ii, j, idist, idistLya, iband, iband1, iband2, ier
  INTEGER :: ilength(paramLAE%nfile), ilength_max
  CHARACTER(LEN=500) :: ref(paramLAE%nfile), cspace
  CHARACTER(LEN=50)  :: cmor = '# (3i-1)E (3i)S0 (3i+1)S'
  CHARACTER(LEN=100) :: cdist_title(n_distLAE) = (/CHARACTER(LEN=100):: &
       'column density Ngas [cm^-2]', 'metallicity Zc [Zsun]', &
       'the metal column density NgasZc [Zsun cm^-2]', &
       'stellar mass Mstar [Msun]', 'host DM halo mass Mhost [Msun]', &
       'effective radius r_eff [kpc/h]', 'SFR [Msun/yr]', &
       'effective star formation timescale tau_eff'/)
  CHARACTER(LEN=50)  :: cdist_x(n_distLAE) = (/CHARACTER(LEN=50):: &
       'Ngas/cm^-2', 'Zc/Zsun', 'NgasZc/(Zsun cm^-2)', 'Mstar/Msun', &
       'Mhost/Msun', 'r_eff/(kpc/h)', 'SFR/(Msun/yr)', 'tau_eff/yr'/)
  CHARACTER(LEN=20)  :: cdist_y1(n_distLAE) = (/CHARACTER(LEN=20)::  &
       'Ngas', '(Zc/Z_sun)', '(NgasZc/Zsun)', 'Mstar', 'Mhost', &
       'r_eff', 'SFR', 'tau_eff'/)
  CHARACTER(LEN=20)  :: cdist_y2(n_distLAE) = (/CHARACTER(LEN=20)::  &
       'cm^-2', 'Zsun', 'Zsun cm^-2', 'Msun', 'Msun', 'kpc/h', 'M_sun/yr', 'yr'/)
  CHARACTER(LEN=100) :: cdist_y3 = &
       ' h^3 Mpc^-3 dex^-1] (8-13)dn/dx[h^3 Mpc^-3 dex^-1] (14-19)n(>x) [h^3 Mpc^-3]'
  CHARACTER(LEN=100) :: cdist_type = &
       '# (2-4,8-10,14-16)quiescent (5-7,11-13,17-19)starburst'
  CHARACTER(LEN=100) :: cdistLya_title = &
       '# Lya LFs divided in 10 groups according to the'
  CHARACTER(LEN=50)  :: cdistLya_title_y(n_distLya) = (/CHARACTER(LEN=50)::  &
       'column density Ngas[cm^-2]', 'metallicity Z [Zsun]', &
       'stellar mass Mstar [Msun]', 'metal column density NgasZc [cm^-2]'/)
  CHARACTER(LEN=50)  :: cdistLya_y1(n_distLya) = (/CHARACTER(LEN=50)::  &
       '[Ngas/cm^-2](=17.75-22.75)', '[Z/Zsun](=-4.25~+0.75)', &
       '[Mstar/Msun](=2.5-12.5)', '[NgasZc/cm^-2](=13.5-23.5)'/)
  CHARACTER(LEN=10)  :: cdistLya_y2(n_distLya) = (/'0.5', '0.5', '1.0', '1.0'/)
  CHARACTER(LEN=500) :: cdistLya_y3(n_distLya) = (/CHARACTER(LEN=500)::  &
       '17.75-18.25 (3+10i)18.25-18.75 (4+10i)18.75-19.25 '//&
       '(5+10i)19.25-19.75 (6+10i)19.75-20.25 (7+10i)20.25-20.75 '//&
       '(8+10i)20.75-21.25 (9+10i)21.25-21.75 (10+10i)21.75-22.25 '//&
       '(11+10i)22.25-22.75 ', &
       '-4.25~-3.75 (3+10i)-3.75~-3.25 (4+10i)-3.25~-2.75 '//&
       '(5+10i)-2.75~-2.25 (6+10i)-2.25~-1.75 (7+10i)-1.75~-1.25 '//&
       '(8+10i)-1.25~-0.75 (9+10i)-0.75~-0.25 (10+10i)-0.25~+0.25 '//&
       '(11+10i)0.25-0.75', &
       '2.5-3.5 (3+10i)3.5-4.5 (4+10i)4.5-5.5(5+10i)5.5-6.5 (6+10i)6.5-7.5 '//&
       '(7+10i)7.5-8.5 (8+10i)8.5-9.5 (9+10i)9.5-10.5 (10+10i)10.5-11.5 '//&
       '(11+10i)11.5-12.5', &
       '13.5-14.5 (3+10i)14.5-15.5 (4+10i)15.5-16.5 (5+10i)16.5-17.5 '//&
       '(6+10i)17.5-18.5 (7+10i)18.5-19.5 (8+10i)19.5-20.5 (9+10i)20.5-21.5 '//&
       '(10+10i)21.5-22.5 (11+10i)22.5-23.5'/)
  CHARACTER(LEN=100) :: cdistLya = '# (1)x=log10[L_Lya / (ergs/s)] '//&
       '(2-21)log10[dn/dx / (h^3 Mpc^-3)] (22-41)log10[n(>x) / (h^3 Mpc^-3)]'
  CHARACTER(LEN=100) :: cdistLya_type = &
       '# (2-11,22-31)quiescent (12-21,32-41)starburst'
  DOUBLE PRECISION :: redshift
  DOUBLE PRECISION :: z2t, Square ! functions


  tlife = (z2t(param%zsp1) - z2t(mrgp%zp1ar(param%izout-1))) * param%th_yr
          ! the final timestep [yr]
  redshift = param%zsp1 - 1.d0
  DO i = 1, paramLAE%nfile
     ii = i + paramLAE%base
     IF(i == 15) THEN
        ii = paramLAE%base + paramLAE%nfile + 1
     ELSEIF(i == 16) THEN
        ii = paramLAE%base + paramLAE%nfile + nnode + 1
     ENDIF
     write(ii, '(A, F7.4)')     '# output z = ', redshift
     write(ii, '(A, F9.4, A)')  '# t_hubble(z) = ', z2t(param%zsp1)*param%th, ' [Gyr]'
     write(ii, '(A, G13.6, A)') '# tlife =', tlife, ' [yr]'
     write(ii, '(A, I5)')       '# paramLAE%nloop = ', paramLAE%nloop
     write(ii, '(A, F8.3)')     '# param%h = ', param%h
  ENDDO

  i = iLAE%LyaLF ! LyaLF
  write(i, '(A)') '# Lya LFs'
  write(i, '(A)') '# (1)x=log10(L_Lya/[ergs/s]) '//&
       '(2-7)log10[dn/dx / (h^3/Mpc^3/dex)] (8-13)log10[n(>x) / (h^3/Mpc^3)]'
  write(i, '(A)') '# (2-4,8-10)quiescent (5-7,11-13)starburst'
  write(i, '(A)') trim(cmor)//' (i=1-4)'

  i = iLAE%contLF ! continuum LFs
  write(i, '(A)') '# Continuum LFs calculated by using the Schaerer '//&
       'population synthesis model'
  write(i, '(A)') '# (1)M-5logh[ABmag] (2-10,11-19)dn/dM w/o,w/ dust[h^3/Mpc^3/mag]'
  write(i, '(A)') '# (2-4,11-13)quiescent (5-7,14-16)starburst '//&
       '(8-10,17-19)multiple starburst'
  write(i, '(A)') trim(cmor)//' (i=1-6)'
  write(i, '(A, I2, A)') '# index 0-', NWAVE_S03-2, ': 1500A, 2800A, 6563A, '//&
       '4861A, 1216A LFs'

  DO idist = 1, n_distLAE ! distribution functions for physical quantities
     i = iLAE%dist(idist)
     write(i, '(A)') '# '//trim(cdist_title(idist))//' distribution'
     write(i, '(A)') '# (1)x=log10['//trim(cdist_x(idist))//'] '//&
          '(2-7)'//trim(cdist_y1(idist))//'*dn/dx['//&
          trim(cdist_y2(idist))//trim(cdist_y3)
     write(i, '(A)') trim(cdist_type)
     write(i, '(A)') trim(cmor)//' (i=1-6)'
  ENDDO

  DO idistLya = 1, n_distLya ! Lya LFs separated into 10 groups
     i = iLAE%distLya(idistLya)
     write(i, '(A)') trim(cdistLya_title)//' '//trim(cdistLya_title_y(idistLya))
     write(i, '(A)') '# y=log10'//trim(cdistLya_y1(idistLya))//', dy='//&
          trim(cdistLya_y2(idistLya))
     write(i, '(A)') '# (2+10i)'//trim(cdistLya_y3(idistLya))//' (i=0-3)'
     write(i, '(A)') trim(cdistLya)
     write(i, '(A)') trim(cdistLya_type)
  ENDDO

  ! i = iLAE%eachLAE
  i = paramLAE%base + paramLAE%nfile + 1
  write(i, '(2(A, G10.3), A)') '# physical properties of the galaxies '//&
       'with L(Lya)^int >=', thLAE%LLya, ' [erg/s/h^2] and Mstar >=', &
       thLAE%Mstar, ' [Msun]'
  write(i, '(A)') '# '//&
       ! --- integers
       '(1)ID '//&
       '(2,3)mor,mor_d[1:E,2:S0,3:S] '//&
       '(4)1:starburst,2:quiescent '//&
       '(5)1:quiescent,2:starburst,3:multiple-merger '//&
       '(6)0:satellite,1:central '//&
       '(7)iforest '//&
       '(8)mpi(gal) '//&
       '(9)mpi(host) '//&
       '(10)ID(host) '//&
       '(11)ID(central) '//&
       ' --- space --- '
  write(i, '(A)') '# '//&
       ! --- double precisions
       '(12)NgasZc[cm^-2] '//&
       '(13)sum(quiescent) or sum*[max(1.d+7, fphaseLya*{t_brst/(alpha+'//&
       'beta)})/tlife](starburst)[h^3/Mpc^3] '//&
       '(14)L(Lya)^int[erg/s/h^2] '//&
       '(15)Mstar[Msun] '//&
       '(16)reff[kpc/h] '//&
       '(17)Zc/Zsun '//&
       '(18)Mhost[Msun] '//&
       '(19)beta '//&
       '(20)tau^eff[yr] for quiescent or tau_star[yr] for starburst '//&
       '(21)Mgas[Msun] '//&
       ' --- space --- '
  write(i, '(A)') '# '//&
       '(22)Mstar^pre[Msun] '//&
       '(23)Mprog[Msun] '//&
       '(24)Vc[km/s] '//&
       '(25)time after the onset of SF[yr] '//&
       '(26)Vc(host halo)[km/s] '//&
       '(27)Vbulge[km/s] for E/S0 or Vdisk[km/s] for S '//&
       '(28)SFR[Msun/yr] '//&
       '(29)tGW[yr] '//&
       '(30)E(B-V) '//&
       '(31)N_LyC[1/s/h^2] '//&
       ' --- space --- '
  write(i, '(A)') '# '//&
       '(32)Mgas prior to SF[Msun] '//&
       '(33)mean SFR past 10Myr[Msun/yr] '//&
       '(34)mass-weighted age[yr] '//&
       '(35,36)M(1500A)-5logh[ABmag] w/o, w/ dust '//&
       '(37,38)M(2800A)-5logh[ABmag] w/o, w/ dust '//&
       '(39,40)M(6563A)-5logh[ABmag] w/o, w/ dust '//&
       '(41,42)M(4861A)-5logh[ABmag] w/o, w/ dust '//&
       '(43,44)M(1216A)-5logh[ABmag] w/o, w/ dust '//&
       ' --- space --- '
!!$  write(i, '(A)') '# '//&
!!$       '(45)L_Br^AGN[erg/s] '//&
!!$       '(46)Eddington ratio in B-band '//&
!!$       '(47)accretion rate to SMBH[Msun/yr]'//&
!!$       ' --- space --- '
!!$  ii = 48 ! the last row + 1
  ii = 45 ! the last row + 1
  write(i, '(A, $)') '#'
  DO iband = 1, param%nwave
     iband1 = 2*(iband-1) + ii; iband2 = param%iwave(iband)
     IF(iband1 >= 100) THEN
        write(i, '(2(A, I3), $)') ' (', iband1, ',', iband1+1
     ELSE
        write(i, '(2(A, I2), $)') ' (', iband1, ',', iband1+1
     ENDIF
     IF((iband2 >= 50 .and. iband2 <= 58) .or. &
          (iband2 >= 50+param%nwave .and. iband2 <= 58+param%nwave)) THEN
        write(i, '(A, $)') ')log10['//trim(ssp%bandname(iband2))//'/'
        IF(iband2 == 50 .or. iband2 == 50+param%nwave) THEN
           write(i, '(A, $)') '(photons/s/h^2)] w/o, w/ dust'
        ELSE
           write(i, '(A, $)') '(erg/s/A/h^2)] w/o, w/ dust'
        ENDIF
     ELSE
        write(i, '(A, $)') ')M('//trim(ssp%bandname(iband2))//&
             ')-5logh[ABmag] w/o, w/ dust'
     ENDIF
  ENDDO
  write(i, *)

  ! i = iLAE%each
  i = paramLAE%base + paramLAE%nfile + nnode + 1
  write(i, '(A, F7.2, A)') '# physical properties of the '//&
       'galaxies with M(UVint)-5logh <= ', thLAE%magUV, ' [ABmag]'
  write(i, '(A)') '# '//&
       ! --- integers
       '(1)ID '//&
       '(2,3)mor,mor_d[1:E,2:S0,3:S] '//&
       '(4)1:starburst,2:quiescent '//&
       '(5)1:quiescent,2:starburst,3:multiple-merger '//&
       '(6)0:satellite,1:central '//&
       '(7)iforest '//&
       '(8)mpi(gal) '//&
       '(9)mpi(host) '//&
       '(10)ID(host) '//&
       '(11)ID(central) '//&
       ' --- space --- '
  write(i, '(A)') '# '//&
       ! --- double precisions
       '(12)NgasZc[cm^-2] '//&
       '(13)sum(quiescent) or sum*[max(1.d+7, fphaseLya*{t_brst/(alpha+'//&
       'beta)})/tlife](starburst)[h^3/Mpc^3] '//&
       '(14)L(Lya)^int[erg/s/h^2] '//&
       '(15)Mstar[Msun] '//&
       '(16)reff[kpc/h] '//&
       '(17)Zc/Zsun '//&
       '(18)Mhost[Msun] '//&
       '(19)beta '//&
       '(20)tau^eff[yr] for quiescent or tau_star[yr] for starburst '//&
       '(21)Mgas[Msun] '//&
       ' --- space --- '
  write(i, '(A)') '# '//&
       '(22)Mstar^pre[Msun] '//&
       '(23)Mprog[Msun] '//&
       '(24)Vc[km/s] '//&
       '(25)time after the onset of SF[yr] '//&
       '(26)Vc(host halo)[km/s] '//&
       '(27)Vbulge[km/s] for E/S0 or Vdisk[km/s] for S '//&
       '(28)SFR[Msun/yr] '//&
       '(29)tGW[yr] '//&
       '(30)E(B-V) '//&
       '(31)N_LyC[1/s/h^2] '//&
       ' --- space --- '
  write(i, '(A)') '# '//&
       '(32)Mgas prior to SF[Msun] '//&
       '(33)mean SFR past 10Myr[Msun/yr] '//&
       '(34)mass-weighted age[yr] '//&
       '(35,36)M(1500A)-5logh[ABmag] w/o, w/ dust '//&
       '(37,38)M(2800A)-5logh[ABmag] w/o, w/ dust '//&
       '(39,40)M(6563A)-5logh[ABmag] w/o, w/ dust '//&
       '(41,42)M(4861A)-5logh[ABmag] w/o, w/ dust '//&
       '(43,44)M(1216A)-5logh[ABmag] w/o, w/ dust '//&
       ' --- space --- '
!!$  write(i, '(A)') '# '//&
!!$       '(45)L_Br^AGN[erg/s] '//&
!!$       '(46)Eddington ratio in B-band '//&
!!$       '(47)accretion rate to SMBH[Msun/yr]'//&
!!$       ' --- space --- '
!!$  ii = 48 ! the last row + 1
  ii = 45 ! the last row + 1
  write(i, '(A, $)') '#'
  DO iband = 1, param%nwave
     iband1 = 2*(iband-1) + ii; iband2 = param%iwave(iband)
     IF(iband1 >= 100) THEN
        write(i, '(2(A, I3), $)') ' (', iband1, ',', iband1+1
     ELSE
        write(i, '(2(A, I2), $)') ' (', iband1, ',', iband1+1
     ENDIF
     IF((iband2 >= 50 .and. iband2 <= 58) .or. &
          (iband2 >= 50+param%nwave .and. iband2 <= 58+param%nwave)) THEN
        write(i, '(A, $)') ')log10['//trim(ssp%bandname(iband2))//'/'
        IF(iband2 == 50 .or. iband2 == 58+param%nwave) THEN
           write(i, '(A, $)') '(photons/s/h^2)] w/o, w/ dust'
        ELSE
           write(i, '(A, $)') '(erg/s/A/h^2)] w/o, w/ dust'
        ENDIF
     ELSE
        write(i, '(A, $)') ')M('//trim(ssp%bandname(iband2))//&
             ')-5logh[ABmag] w/o, w/ dust'
     ENDIF
  ENDDO
  write(i, *)

  write(ilog, *); write(ilog, *)
  write(ilog, '(A)') '# +--------------------------------------------+'
  write(ilog, '(A)') '# |          LAE related parameters            |'
  write(ilog, '(A)') '# +--------------------------------------------+'
  write(ilog, '(A, I5)') '# paramLAE%nloop:     ',   paramLAE%nloop
  write(ilog, '(A, I1)') '# paramLAE%type_b:    ',   paramLAE%type_b
  write(ilog, '(A, F6.2)') '# paramLAE%Rtburst:   ', paramLAE%Rtburst
  write(ilog, '(A, F6.2)') '# paramLAE%fphaseLya: ', paramLAE%fphaseLya
  write(ilog, '(A)') '# Thresholds used in the LAE related calculation'
  write(ilog, '(A, G9.2, A)') '# --- thLAE%time:   ', thLAE%time,   ' [yr]'
  write(ilog, '(A, G9.2, A)') '# --- thLAE%LLya:   ', thLAE%LLya,   ' [erg/s/h^2]'
  write(ilog, '(A, G9.2, A)') '# --- thLAE%Mstar:  ', thLAE%Mstar,  ' [Msun]'
  write(ilog, '(A, G9.2, A)') '# --- thLAE%magUV:  ', thLAE%magUV,  ' [ABmag]'
  write(ilog, '(A)') '# filenames used in the LAE related calculation:'
  ref(1)  = ' = intrinsic (i.e., w/o dust extinction) Lya LF calculated '//&
       'in the Mitaka model'
  ref(2)  = ' = continuum LFs w/o and w/ dust extinction'
  ref(3)  = ' = distribution of cold gas column density'
  ref(4)  = ' = distribution of metallicity'
  ref(5)  = ' = distribution of metal column density in cold phase'
  ref(6)  = ' = distribution of stellar mass'
  ref(7)  = ' = distribution of host halo mass'
  ref(8)  = ' = distribution of effective radius'
  ref(9)  = ' = distribution of SFR'
  ref(10) = ' = distribution of effective timescale of SF'
  ref(11) = ' = distribution in Lya-Ngas plane'
  ref(12) = ' = distribution in Lya-Zc plane'
  ref(13) = ' = distribution in Lya-Mstar plane'
  ref(14) = ' = distribution in Lya-NcZc plane'
  ref(15) = ' = physical properties of each Lya-bright galaxies'
  ref(16) = ' = physical properties of each UV-bright galaxies'
  ilength_max = 0
  DO i = 1, paramLAE%nfile
     ilength(i) = len_trim(fnameLAE(i))
     IF(ilength(i) > ilength_max) ilength_max = ilength(i)
  ENDDO
  DO i = 1, paramLAE%nfile
     cspace = ' '
     DO ii = 1, ilength_max - ilength(i)
        cspace = cspace//' '
     ENDDO
     ii = ilength_max - ilength(i)
     write(ilog, '(A, I2, A)') '# --- ', i, ': '//trim(fnameLAE(i))//&
          cspace(1:ii)//trim(ref(i))
  ENDDO
END SUBROUTINE WriteCaptionsForLAE
! ========================================================================
SUBROUTINE CloseFilesForLAE(run_redshift)
  use LAErelated
  implicit none
  INTEGER, INTENT(IN) :: run_redshift ! 1:run at a redshift, 2:run at all redshift
  INTEGER :: i, ier
  CHARACTER(LEN=50) :: cerr = '# (LAE)CloseFilesForLAE: fail to deallocate LAE%'

  deallocate(LAE%flag_mag,stat=ier); call CheckIerr(ier, trim(cerr)//'flag_mag')
  deallocate(LAE%lumg,    stat=ier); call CheckIerr(ier, trim(cerr)//'lumg')
  deallocate(LAE%lumg_d,  stat=ier); call CheckIerr(ier, trim(cerr)//'lumg_d')
  deallocate(LAE%mag,     stat=ier); call CheckIerr(ier, trim(cerr)//'mag')
  deallocate(LAE%mag_d,   stat=ier); call CheckIerr(ier, trim(cerr)//'mag_d')
  IF(run_redshift == 1) THEN
     DO i = 1, paramLAE%nfile
        close(i+paramLAE%base)
     ENDDO
  ENDIF
END SUBROUTINE CloseFilesForLAE
! ========================================================================
INTEGER FUNCTION CalBin(x, invstep, base) RESULT(bin)
  implicit none
  DOUBLE PRECISION, INTENT(IN) :: x, invstep, base

  bin = int(anint(invstep * (x - base)))
        ! base+(i+0.5)/invstep <= x < base+(i+1.5)/invstep ==> bin = i+1
        ! the center of bin=i+1: x = base+(i+1)/invstep
END FUNCTION CalBin
! ========================================================================
DOUBLE PRECISION FUNCTION CalT_GW(tau, alp, beta)
! --- the end time of star formation after the onset of the starburst
! --- this is defined by the following equation:
!      \int_0^t_GW SFR(t')dt' = [alpha/(alpha+beta)]Mc0
!      where SFR(t') = [Mc0/taust] exp[-t'/taust]
! --- t_GW = taust * log(1 + alpha/beta)
  implicit none
  DOUBLE PRECISION, INTENT(IN) :: tau, alp, beta
  DOUBLE PRECISION :: tmp

  tmp = 1.d0 + alp / beta
  CalT_GW = tau * log(tmp)
END FUNCTION CalT_GW
! ========================================================================
SUBROUTINE lum_burst_LAE(zplus1pc, zplus1, taustar, lumtmp2, mz, mt, tel, &
                         t_GW, metal0)
  use global_var
  implicit none
  DOUBLE PRECISION, INTENT(IN)    :: zplus1pc, zplus1, taustar
  DOUBLE PRECISION, INTENT(INOUT) :: lumtmp2(param%nwave)
  DOUBLE PRECISION, INTENT(INOUT) :: mz, mt
  DOUBLE PRECISION, INTENT(INOUT) :: tel, t_GW
  DOUBLE PRECISION, INTENT(IN)    :: metal0
  INTEGER, PARAMETER :: b_or_q = 1
  INTEGER :: k, i, tt, ttp1, zli, zlip1
  DOUBLE PRECISION :: tstart, tfinal, tsinv, tfinv, tauinv
  DOUBLE PRECISION :: ttm1, ttm2, ttm3
  DOUBLE PRECISION :: tstep, tbranch
  DOUBLE PRECISION :: sfrt, sfrtt
  DOUBLE PRECISION :: fxyi0, fxyi1, fxy(param%nwave)
  DOUBLE PRECISION :: zt
  DOUBLE PRECISION, PARAMETER :: EPS  = 1.d-8
  DOUBLE PRECISION, PARAMETER :: EPS2 = 1.d-14
  DOUBLE PRECISION :: tmp1, tmp2
  INTEGER :: IntTimeSSP, IntMetalSSP ! functions
  DOUBLE PRECISION :: z2t ! function

  lumtmp2(:) = 0.d0; mz = 0.d0; mt = 0.d0
  IF(zplus1pc == zplus1)  return

  tauinv  = 1.d0 / taustar

  tstart = z2t(zplus1) + tel
  tfinal = min(z2t(zplus1pc), tstart + t_GW)
  tsinv  = const%tout - tstart; tfinv = const%tout - tfinal
  tbranch = tfinal - tstart

  zt = metal0; zli = IntMetalSSP(metal0)
  zlip1 = zli + 1
  sfrtt = 0.d0; ttm1 = 0.d0; ttm2 = 0.d0
  loop1: DO WHILE (abs(tbranch - ttm2) > EPS)
     tt = IntTimeSSP(tsinv-ttm1)
     tstep = ssp%time(tt+1) - ssp%time(tt)
     IF((ttm2 + tstep) > tbranch) tstep = tbranch - ttm2

     ttm1 = ttm2
     ttm2 = ttm2 + tstep

     tmp1 = -ttm1  * tauinv; tmp2 = -tstep * tauinv
     sfrt = exp(tmp1) * (1.d0 - exp(tmp2))
     IF(sfrt <= EPS2) goto 21
     sfrtt = sfrtt + sfrt

     ttm3 = tsinv - ttm1 - 0.5d0 * tstep
     tt   = IntTimeSSP(ttm3)
     ttp1 = tt + 1

     tmp1 = (ttm3 - ssp%time(tt)) / (ssp%time(ttp1) - ssp%time(tt))
     tmp2 = (zt - ssp%chem(zli)) / (ssp%chem(zlip1) - ssp%chem(zli))
     loop2: DO k = 1, param%nwave
        i = param%iwave(k)
        fxyi0 = (ssp%lumi(b_or_q, i, ttp1, zli)   - ssp%lumi(b_or_q, i, tt, zli))   &
                * tmp1 + ssp%lumi(b_or_q, i, tt, zli)
        fxyi1 = (ssp%lumi(b_or_q, i, ttp1, zlip1) - ssp%lumi(b_or_q, i, tt, zlip1)) &
                * tmp1 + ssp%lumi(b_or_q, i, tt, zlip1)
        IF(zt > ssp%chem(NCHEM)) THEN
           fxy(k) = fxyi1
        ELSEIF(zt < ssp%chem(1)) THEN
           fxy(k) = fxyi0
        ELSE
           fxy(k) = (fxyi1 - fxyi0) * tmp2 + fxyi0
        ENDIF
        lumtmp2(k) = lumtmp2(k) + sfrt * fxy(k)
     ENDDO loop2

     tmp1 = sfrt * fxy(3)
     mz = mz + tmp1 * zt; mt = mt + tmp1 * (tstart + ttm1)
  ENDDO loop1

21 continue
END SUBROUTINE lum_burst_LAE
! ========================================================================
SUBROUTINE optNZ(nwave, lumtmp, lumtmp_d, NZ)
! --- dust optical depth for specified pass-bands
! --- dust geometry: screen (paramLAE%exttype=1), slab (paramLAE%exttype=2)
! --- this subroutine runs in each galaxy calculation
  use global_var; use LAErelated
  implicit none
  INTEGER, INTENT(IN) :: nwave
  DOUBLE PRECISION, INTENT(INOUT) :: lumtmp(param%tnw), lumtmp_d(param%tnw)
  DOUBLE PRECISION, INTENT(IN)    :: NZ ! metal column density [cm^-2]
  INTEGER :: i, j
  DOUBLE PRECISION :: tauV, taud, fescV
  DOUBLE PRECISION :: f_scr, f_slab ! function

  tauV = constLAE%OptV * NZ
  DO i = 1, nwave
     taud = tauV * xi(i)
     IF(paramLAE%exttype == 1) THEN ! screen dust geometry
        fescV = f_scr(taud)
     ELSEIF(paramLAE%exttype == 2) THEN ! slab dust geometry
        fescV = f_slab(taud)
     ENDIF
     j = i + nwave
     lumtmp_d(i) = fescV * lumtmp(i) ! spheroid
     lumtmp_d(j) = fescV * lumtmp(j) ! disk
  ENDDO
END SUBROUTINE optNZ
! ========================================================================
DOUBLE PRECISION FUNCTION CalEB_V(NZ)
! --- this subroutine runs in each galaxy calculation 'NZ' [cm^-2]
  use LAErelated
  implicit none
  DOUBLE PRECISION, INTENT(IN) :: NZ ! [cm^-2] (not in [Zsun cm^-2])
  DOUBLE PRECISION :: tauV, tauB, fescV, fescB, AV, AB
  DOUBLE PRECISION :: f_slab ! function

  tauV = constLAE%OptV * NZ; tauB = constLAE%OptB * NZ
  IF(tauV < 1.d-12) THEN
     CalEB_V = 0.d0
  ELSE
     IF(paramLAE%exttype == 1) THEN ! screen
        AV = constLAE%AV0 * tauV; AB = constLAE%AB0 * tauV
     ELSEIF(paramLAE%exttype == 2) THEN ! slab
        fescV = f_slab(tauV); AV = -2.5d0 * log10(fescV)
        fescB = f_slab(tauB); AB = -2.5d0 * log10(fescB)
     ENDIF
     CalEB_V = AB - AV
  ENDIF
END FUNCTION CalEB_V
! ========================================================================
DOUBLE PRECISION FUNCTION f_scr(tau)
  implicit none
  DOUBLE PRECISION, INTENT(IN) :: tau

  IF(tau < 1.d-5) THEN
     f_scr = 1.d0 - tau
  ELSE
     f_scr = exp(-tau)
  ENDIF
END FUNCTION f_scr
! ========================================================================
DOUBLE PRECISION FUNCTION f_slab(tau)
  implicit none
  DOUBLE PRECISION, INTENT(IN) :: tau

  IF(tau > 20.d0) THEN
     f_slab = 1.d0 / tau
  ELSEIF(tau < 1.d-8) THEN
     f_slab = 1.d0 - (0.5d0 - tau / 6.d0) * tau
  ELSE
     f_slab = (1.d0 - exp(-tau)) / tau
  ENDIF
END FUNCTION f_slab
! ========================================================================
DOUBLE PRECISION FUNCTION f_sand(tau)
   implicit none
   DOUBLE PRECISION, INTENT(IN) :: tau
   DOUBLE PRECISION, PARAMETER :: delta = 0.5d0

   IF(tau > 20.d0) THEN
      f_sand = (0.5d0 * (1-delta) + delta / tau)
   ELSE IF(tau < 1.d-8) THEN
      f_sand = 1.d0 - tau * (0.5d0 +  0.25 * tau * (1.d0 - delta / 3.d0))
   ELSE
      f_sand = (0.5d0 * (1-delta)) * (1.d0 + exp(-tau)) &
               + (delta / tau) * (1.d0 - exp(-tau))
   ENDIF
   END FUNCTION f_sand
! ========================================================================
DOUBLE PRECISION FUNCTION LinearInterp(x, x1, x2, y1, y2, inv) RESULT(y)
  implicit none
  DOUBLE PRECISION, INTENT(IN) :: x, x1, x2, y1, y2
  DOUBLE PRECISION, INTENT(IN) :: inv ! inv = 1/(x1-x2)

  y = ((y1 - y2) * x + x1*y2 - x2*y1) * inv
END FUNCTION LinearInterp
! ========================================================================
SUBROUTINE GatherResultsLAE(run_redshift, inode)
  use LAErelated
  implicit none

  include 'mpif.h'

  INTEGER, INTENT(IN) :: run_redshift ! 1:run at a redshift, 2:run at all redshift
  INTEGER, INTENT(IN) :: inode
  INTEGER :: i, ier, istat
  CHARACTER(LEN=5000) :: command
  TYPE(DistributionForLAE) :: tmp_dist1 ! for LyaLF, distLAE(1:n_distLAE)
  TYPE(LuminosityFunction) :: tmp_dist2 ! for LF_S03(1:NWAVE_S03)
  TYPE(Lya_x_Distribution) :: tmp_dist3 ! for distLya(1:n_distLya)

  IF(run_redshift == 1 .and. inode == 0) THEN
     ! --- LyA-bright galaxy catalog
     write(fnameLAE_catalog_LyA(len_trim(fnameLAE_catalog_LyA)-5:&
          len_trim(fnameLAE_catalog_LyA)-4), '(A2)') '??'
     command = 'cat '//trim(fnameLAE_catalog_LyA)//' > '//&
          trim(fnameLAE(15))
     call system(trim(command))
     print '(A)', '# (LAE)GatherResultsLAE executing the following '// &
          'command: '//trim(command)

     ! --- UV-bright galaxy catalog
     write(fnameLAE_catalog_UV(len_trim(fnameLAE_catalog_UV)-5:&
          len_trim(fnameLAE_catalog_UV)-4), '(A2)') '??'
     command = 'cat '//trim(fnameLAE_catalog_UV)//' > '//&
          trim(fnameLAE(16))
     call system(trim(command))
     print '(A)', '# (LAE)GatherResultsLAE: executing the following '// &
          'command: '//trim(command)

     deallocate(fnameLAE, stat=ier); call CheckIerr(ier, &
          '# (LAE)GatherResultsLAE: fail to deallocate: fnameLAE')
  ENDIF

  ! --- Lya LF
  tmp_dist1%count(:,:,:) = 0.d0
  call MPI_ALLREDUCE(LyaLF%count, tmp_dist1%count, 2*3*NLF1, &
                     MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD, istat)
  LyaLF%count = tmp_dist1%count

  ! --- S03 LFs
  DO i = 1, NWAVE_S03
     tmp_dist2%count(:,:,:) = 0.d0; tmp_dist2%countd(:,:,:) = 0.d0
     call MPI_ALLREDUCE(LF_S03(i)%count,  tmp_dist2%count,  3*3*NLF1, &
                        MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD, istat)
     call MPI_ALLREDUCE(LF_S03(i)%countd, tmp_dist2%countd, 3*3*NLF1, &
                        MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD, istat)
     LF_S03(i)%count = tmp_dist2%count; LF_S03(i)%countd = tmp_dist2%countd
  ENDDO

  ! --- distribution functions
  DO i = 1, n_distLAE
     tmp_dist1%count(:,:,:) = 0.d0; tmp_dist1%wcount(:,:,:) = 0.d0
     call MPI_ALLREDUCE(distLAE(i)%count,  tmp_dist1%count,  2*3*NLF1, &
                        MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD, istat)
     call MPI_ALLREDUCE(distLAE(i)%wcount, tmp_dist1%wcount, 2*3*NLF1, &
                        MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD, istat)
     distLAE(i)%count = tmp_dist1%count; distLAE(i)%wcount = tmp_dist1%wcount
  ENDDO

  ! --- 2D distribution functions
  DO i = 1, n_distLya
     tmp_dist3%count(:,:,:) = 0.d0
     call MPI_ALLREDUCE(distLya(i)%count, tmp_dist3%count, 2*NLF2*NLF1, &
                        MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD, istat)
     distLya(i)%count = tmp_dist3%count
  ENDDO
END SUBROUTINE GatherResultsLAE
! ========================================================================
