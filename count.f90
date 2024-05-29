!!$ === ifc version 7 is required ===
!!$ === for HDF/SDF ===
!!$ === including the effects of intergalactic absorption (03/12/14)
!!$ count050403.f based on count040811jhk.f
!!$ see Astronomical Quantities p.150
!!$ 1 Jy = 1e-26 W/m^2/Hz
!!$ band lam[micron] Flam [W/m^2/micron] Fnu for Vega [Jy] lam*Flam[W/m^2]
!!$ J    1.215       3.31e-9             1630              4.02165e-9
!!$ H    1.654       1.15e-9             1050              1.90210e-9
!!$ K    2.179       4.14e-10             655              9.02106e-10
!!$ The units of the EBL should be nW/m^2/sr .
  ! +-----------------------------------------------+
  ! |  Usage: (ex) ./count.out 140hMpc_2048_Planck2 |
  ! +-----------------------------------------------+
PROGRAM count
  implicit none

  DOUBLE PRECISION :: PI
  INTEGER :: i, j, k, l, m
  DOUBLE PRECISION :: x, y
  INTEGER :: ix, iy

  INTEGER, PARAMETER :: Nwave = 9 ! 1:2MASS_Ks, 2:IRACch1, 3:IRACch2, 4:CIBER_I, 5:CIBER_H,
                                   !  6--8:AKARI_N{2,3,4}, 9--13:Suprime_{B,V,R,i',z'}
  INTEGER, PARAMETER :: Ngalmax = 1000000
  INTEGER :: izslice, nzslice, iz, iz0, izp1
  INTEGER :: ihalo, igal, iforest, iforest0, irow, iwave
  INTEGER :: ngal, ngal_tot, nhalo_tot
  INTEGER :: iread = 20, iwrite = 30, ier, ier0
  INTEGER, PARAMETER :: N_COLUMN_I = 8, N_COLUMN_I8 = 2, N_COLUMN_DP = 43
  INTEGER :: ic, ixin(N_COLUMN_I)
  INTEGER(KIND=8) :: ixin8(N_COLUMN_I8)
  DOUBLE PRECISION :: xin(N_COLUMN_DP)
  DOUBLE PRECISION :: rad2min, rad2sec
  CHARACTER(LEN=3)   :: ci, ci2
  CHARACTER(LEN=100) :: c_af  = '# Allocation Fault'
  CHARACTER(LEN=100) :: c_daf = '# Deallocation Fault'
  CHARACTER(LEN=100) :: c_fop = '# Fail to open file = '
  ! CHARACTER(LEN=100) :: fbase = '/home/kobayasi/nugc/data/NIRB/'
  ! CHARACTER(LEN=100) :: fbase_r = '/work2/kobayasi/data/NIRB/2015May/'
  ! CHARACTER(LEN=100) :: fbase_w = '/work2/kobayasi/data/NIRB_map/2015May/'
  CHARACTER(LEN=100) :: fbase_r = 'data/'
  CHARACTER(LEN=100) :: fbase_w = 'data/'

  INTEGER, PARAMETER :: N_FOREST = 512
  INTEGER, PARAMETER :: iz_start = 1
  INTEGER :: i_nbody
  INTEGER :: iz_end(2) = (/57, 49/)
  CHARACTER(LEN=100) :: fbase_nbody(2) = &
       (/'/work/140hMpc_2048_Planck_MRGT2/', '/work/280hMpc_2048_Planck_MRGT/'/)
  CHARACTER(LEN=100) :: fbase_lc(2)    = &
       (/'/work/140hMpc_2048_Planck_Lightcone/', &
         '/work/280hMpc_2048_Planck_Lightcone_modified/'/)
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: zslice
  CHARACTER(LEN=100) :: cn
  CHARACTER(LEN=500) :: fname_in, fname_out, fname_wind, fname_wpos
  INTEGER :: irep_gal, irep_forest, irep_slice
  CHARACTER(LEN=1000) :: buf

  ! --- Nbody lightcone data
  TYPE data_wprm
     ! Nhalo(i): # of mpi in the redshift slice of "i"
     ! mpi0(i) : the total # of mpi in the redshift slice of < "i"
     INTEGER :: Nhalo(100)
     INTEGER :: mpi0(100)
  END type data_wprm
  TYPE(data_wprm) :: wprm

  ! --- cosmology related quantities
  TYPE data_cosmology
     INTEGER :: nz
     DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: z
     DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: drad
     DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: dlum
     DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: dmod
     ! z(1:nz)   : redshift
     ! drad(1:nz): comoving radial distance [Mpc]
     ! dlum(1:nz): luminosity distance [Mpc/h]
     ! dmod(1:nz): distance modulus for M-5logh
  END type data_cosmology
  TYPE(data_cosmology) :: cos
  DOUBLE PRECISION, PARAMETER :: H0 = 0.68d0 ! present hubble parameter

  ! --- positional information from Nbody
  TYPE position
     INTEGER :: nhalo
     INTEGER(KIND=8), DIMENSION(:), ALLOCATABLE :: mpi ! mpi(1:nhalo)
     REAL, DIMENSION(:,:), ALLOCATABLE :: x ! position, x(1:nhalo, 1:3)
                           ! x(:,1):   comoving radial distance [Mpc]
                           ! x(:,2:3): comoving tangential coordinates [Mpc]
  END type position
  TYPE(position) :: p

  DOUBLE PRECISION :: drad, theta0, theta_size0
  TYPE galaxy
     INTEGER :: ngal
     INTEGER, DIMENSION(:), ALLOCATABLE :: mord
     INTEGER(KIND=8),  DIMENSION(:), ALLOCATABLE :: mpi 
     DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: rdomi, SFR, Mstar, Mhalo
     DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: magd
     ! mord(1:ngal) : morphology[1:E,2:S0,3:S]
     ! mpi(1:ngal)  : marker particle index
     ! rdomi(1:ngal): effective radius of dominant comp.
     ! magd(1:ngal,1:Nwave)
  END type galaxy
  TYPE(galaxy), ALLOCATABLE :: gal(:) ! gal(1:N_FOREST)
  ! --- tentative arrays for the elements of the structure galaxy
  INTEGER :: mord(Ngalmax)
  INTEGER(KIND=8)  :: mpi(Ngalmax)
  DOUBLE PRECISION :: magd(Ngalmax,Nwave), rdomi(Ngalmax), SFR(Ngalmax), Mstar(Ngalmax), Mhalo(NgalMax)
  DOUBLE PRECISION :: z, dmod, theta_x, theta_y, theta_size, tmp
  DOUBLE PRECISION :: LinearInterp ! function


  PI = acos(-1.d0)
  rad2min = 180.d0 * 60.d0 / PI ! conversion factor from radian to arcmin
  rad2sec = rad2min * 60.d0     ! conversion factor from radian to arcsec

  ! --- reading the # of halos for lightcone written in each N-body data slice
  call getarg(1, cn)
  IF(index(cn, "140hMpc_2048_Planck2") /= 0) THEN ! 140 Mpc/h w/ 2048^3
     i_nbody = 1
     open(iread, file = trim(fbase_lc(i_nbody))//'140hMpc_2048_Planck.wprm', status='old')
  ELSEIF(index(cn, "280hMpc_2048_Planck") /= 0) THEN ! 280 Mpc/h w/ 2048^3
     i_nbody = 2
     open(iread, file = trim(fbase_lc(i_nbody))//'280hMpc_2048_Planck.wprm', status='old')
  ENDIF
  ier = 0; wprm%mpi0(:) = 0; wprm%Nhalo(:) = 0
  DO WHILE(ier == 0)
     read(iread, '(A)', iostat=ier) buf
     IF(ier == 0) THEN
        read(buf, *) iz, wprm%mpi0(iz), wprm%Nhalo(iz)
     ENDIF
  ENDDO
  close(iread)


  ! --- reading the redshifts of N-body data slice
  open(iread, file = trim(fbase_nbody(i_nbody))//'sspfiles/redshifts.dat', status='old')
  ier = 0; irow = 1
  allocate(zslice(iz_end(i_nbody)), stat=ier); call CheckIOE(ier, &
       trim(c_af)//': zslice'); zslice(:) = 0.d0
  DO WHILE(ier == 0)
     read(iread, '(A)', iostat=ier) buf
     IF(irow <= 2) THEN
        irow = irow + 1
     ELSE
        IF(ier == 0) THEN
           read(buf, *) iz, zslice(iz)
        ENDIF
     ENDIF
  ENDDO
  close(iread)


  ! --- reading the distance scales
  ! --- first, reading the number of redshift in the data file and allocate
  !       the relevant arrays
  open(iread, file = 'distance_Planck.dat', status='old')
  ier = 0; irow = 1; cos%nz = 0
  DO WHILE(ier == 0)
     read(iread, '(A)', iostat=ier) buf
     IF(irow <= 3) THEN
        irow = irow + 1
     ELSE
        IF(ier == 0) cos%nz = cos%nz + 1
     ENDIF
  ENDDO
  close(iread)
  allocate(cos%z(cos%nz),    stat=ier); call CheckIOE(ier, trim(c_af)//': cos%z')
  allocate(cos%drad(cos%nz), stat=ier); call CheckIOE(ier, trim(c_af)//': cos%drad')
  allocate(cos%dlum(cos%nz), stat=ier); call CheckIOE(ier, trim(c_af)//': cos%dlum')
  allocate(cos%dmod(cos%nz), stat=ier); call CheckIOE(ier, trim(c_af)//': cos%dmod')
  cos%z(:) = 0.d0; cos%drad(:) = 0.d0; cos%dlum(:) = 0.d0; cos%dmod(:) = 0.d0

  ! --- then, reading the data of distance
  open(iread, file = 'distance_Planck.dat', status='old')
  ier = 0; irow = 1; iz = 1
  DO WHILE(ier == 0)
     read(iread, '(A)', iostat=ier) buf
     IF(irow <= 3) THEN
        irow = irow + 1
     ELSE
        IF(ier == 0) THEN
           read(buf, *) (xin(ic), ic=1,14)
           cos%z(iz)    = xin(1) ! redshift
           cos%drad(iz) = xin(4) * 1.d+3 / H0
                          ! comoving radial distance [Mpc] (not in [Mpc/h])
           cos%dlum(iz) = xin(5) ! luminosity distance [Gpc/h]
           cos%dmod(iz) = 5.d0 * log10(cos%dlum(iz) * 1.d+8)
                          ! distance modulus for M-5logh
           iz = iz + 1
        ENDIF
     ENDIF
  ENDDO
  close(iread)


  IF(i_nbody == 1) THEN ! 140hMpc_2048_Planck2
     fname_in  = trim(fbase_r)//'best_forLC_1_001.dat'
     fname_out = trim(fbase_w)//'best_forLC_1_LC.dat'
  ELSEIF(i_nbody == 2) THEN ! 280hMpc_2048_Planck
     fname_in  = trim(fbase_r)//'best_280hMpc_2048_Planck_1_001.dat'
     fname_out = trim(fbase_w)//'best_280hMpc_2048_Planck_1_background.dat'
  ENDIF
  open(iwrite, file = trim(fname_out), iostat=ier); call &
       CheckIOE(ier, trim(c_fop)//trim(fname_out))
  irep_gal = index(fname_in, '001.dat')
  write(ci, '(I3.3)') iz_start; write(ci2, '(I3.3)') iz_end(i_nbody)
  write(iwrite, '(A)') '# Data for Lightcone Observation calculated from '//&
       fname_in(1:irep_gal-1)//'***.dat (*** = '//ci//'--'//ci2//')'
  !write(iwrite, '(A)') '# (1)izslice (2)mpi(gal) [space] (3)mor_d[1:E,2:S0,3:S] '//&
  !     '(4)z (5)comoving radial distance [Mpc] [space] (6--18)Mag-5logh(6:2MASS_Ks,'//&
  !     '7:IRACch1,8:IRACch2,9:CIBER_I,10:CIBER_H,11:AKARI_N2,12:AKARI_N3,13:AKARI_N4,'//&
  !     '14:Suprime_B,15:Suprime_V,16:Suprime_R,17:Suprime_ip,18:Suprime_zp-band) w/ '//&
  !     'dust[AB mag] (19)distance modulus from Mag-5logh to mag [space] '//&
  !     '(20,21)theta_x,y[arcmin] (22)theta_size[arcsec]'
  write(iwrite, '(A)') '# (1)morphology[1:E,2:S0,3:S] '//&
       '(2)redshift (3)comoving radial distance [Mpc] (4) halo mass [Msun] (5) stellar mass [Msun] '//&
       '(6) SFR [Msun/yr] (7--12)Mag-5logh (7:SDSS_up,'//&
       '8:SDSS_gp, 9:SDSS_rp, 10:SDSS_ip, 11:SDSS_zp, 12:2MASS_Ks-band'//&
       '[AB mag] (13)distance modulus from Mag-5logh to mag [space] '//&
       '(14,15)theta_x,y[arcmin] (16)theta_size[arcsec]'


  IF(i_nbody == 1) THEN ! 140hMpc_2048_Planck2
     fname_wind = trim(fbase_lc(i_nbody))//'140hMpc_2048_Planck.00000000.001'
  ELSEIF(i_nbody == 2) THEN ! 280hMpc_2048_Planck
     fname_wind = trim(fbase_lc(i_nbody))//'280hMpc_2048_Planck.00000000.001'
  ENDIF
  fname_wpos = trim(fname_wind)//'.wpos'
  fname_wind = trim(fname_wind)//'.wind'
  irep_slice  = index(fname_wind, '001')
  irep_forest = irep_slice - 4
  nzslice = iz_end(i_nbody) - iz_start + 1
  DO_zslice: DO izslice = iz_start, iz_end(i_nbody)
     print '(2(A, I2), A, I12, A, F10.3)', &
          '(', izslice-iz_start+1, '/', nzslice, ') Nhalo = ', &
          wprm%Nhalo(izslice), ' at z = ', zslice(izslice)
     IF(wprm%Nhalo(izslice) > 0) THEN
        ! --- reading the data of galaxy
        write(ci, '(I3.3)') izslice
        fname_in(irep_gal:irep_gal+2) = ci
        open(iread, file=trim(fname_in), status='old', iostat=ier); call &
             CheckIOE(ier, trim(c_fop)//trim(fname_in))
        print '(A)', ' --- reading the data of galaxy from '//trim(fname_in)
        allocate(gal(N_FOREST), stat=ier); call CheckIOE(ier, trim(c_af)//': gal')
        gal(:)%ngal = 0
        ier0 = 0; irow = 1; iforest0 = 0; ngal_tot = 0; mord(:) = 0; mpi(:) = 0
        rdomi(:) = 0.d0; magd(:,:) = 0.d0
        DO WHILE(ier0 == 0)
           read(iread, '(A)', iostat=ier0) buf
           ! IF(izslice == 7) print '(A)', trim(buf)
           IF(irow <= 2) THEN
              irow = irow + 1
           ELSE
              IF(ier0 == 0) THEN
                 read(buf, *) (ixin(ic), ic=1,N_COLUMN_I), &
                      (ixin8(ic), ic=1,N_COLUMN_I8), &
                      (xin(ic), ic=1,N_COLUMN_DP)
                 iforest = ixin(7) + 1 ! 0-offset --> 1-offset
                 IF(iforest /= iforest0) THEN
                    IF(iforest0 /= 0) THEN
                       ngal = gal(iforest0)%ngal
                       ! print '(A, 2I4, I10)', 'A: ', ier0, iforest0, ngal
                       allocate(gal(iforest0)%mord(ngal),       stat=ier); call  &
                            CheckIOE(ier, trim(c_af)//': gal(:)%mord')
                       allocate(gal(iforest0)%mpi(ngal),        stat=ier); call &
                            CheckIOE(ier, trim(c_af)//': gal(:)%mpi')
                       allocate(gal(iforest0)%rdomi(ngal),      stat=ier); call &
                            CheckIOE(ier, trim(c_af)//': gal(:)%rdomi')
                       allocate(gal(iforest0)%magd(ngal,Nwave), stat=ier); call &
                            CheckIOE(ier, trim(c_af)//': gal(:)%magd')
                       gal(iforest0)%mord(1:ngal)   = mord(1:ngal)
                       gal(iforest0)%mpi(1:ngal)    = mpi(1:ngal)
                       gal(iforest0)%rdomi(1:ngal)  = rdomi(1:ngal)
                       gal(iforest0)%magd(1:ngal,:) = magd(1:ngal,:)
                       mord(:) = 0; mpi(:) = 0; rdomi(:) = 0.d0; magd(:,:) = 0.d0
                    ENDIF
                    iforest0 = iforest; gal(iforest)%ngal = 0
                 ENDIF
                 gal(iforest)%ngal = gal(iforest)%ngal + 1
                 ngal     = gal(iforest)%ngal ! # of galaxies in the forest
                 ngal_tot = ngal_tot + 1      ! # of galaxies in all forests
                 IF(ngal > Ngalmax) &
                      print '(A)', '# Ngalmax is insufficient!!'
                 ! --- mord, rdomi, mpi, mags are tentatively substituted into
                 !       the arrays not in the arrays in the structure of gal(:)
                 mord(ngal) = ixin(4)  ! morphology from (B/T)dust
                 mpi(ngal)  = ixin8(1) ! mpi(gal) (not mpi(host))
                 IF(mord(ngal) == 1 & .or. mord(ngal) == 2) THEN
                    rdomi(ngal) = xin(12) ! rb0 [kpc/h]
                 ELSE
                    rdomi(ngal) = xin(14) ! rd0 [kpc/h]
                 ENDIF
                 Mhalo(ngal) = xin(5) 
                 Mstar(ngal) = xin(2)+xin(3)
                 SFR(ngal) = xin(21)
                 ! IF(ixin(4) /= 3) print '(2I2, 3G10.3)', &
                 !      ixin(4), mord(ngal), xin(11), xin(13), rdomi(ngal)
                 DO iwave = 1, Nwave
                    magd(ngal, iwave) = xin(25+iwave*2)
                 END DO
                 !magd(ngal, 1) = xin(30) ! 2MASS_Ks w/ dust
                 !magd(ngal, 2) = xin(32) ! IRACch1  w/ dust
                 !magd(ngal, 3) = xin(34) ! IRACch2  w/ dust
                 !magd(ngal, 4) = xin(36) ! CIBER_I  w/ dust
                 !magd(ngal, 5) = xin(38) ! CIBER_H  w/ dust
                 !magd(ngal, 6) = xin(40) ! AKARI_N2 w/ dust
                 !magd(ngal, 7) = xin(42) ! AKARI_N3 w/ dust
                 !magd(ngal, 8) = xin(44) ! AKARI_N4 w/ dust
                 !magd(ngal, 9) = xin(46) ! Suprime_B w/ dust
                 !magd(ngal,10) = xin(48) ! Suprime_V w/ dust
                 !magd(ngal,11) = xin(50) ! Suprime_R w/ dust
                 !magd(ngal,12) = xin(52) ! Suprime_ip w/ dust
                 !magd(ngal,13) = xin(54) ! Suprime_zp w/ dust
               ELSE
                 ngal = gal(iforest0)%ngal
                 ! print '(A, 2I4, I10)', 'B: ', ier0, iforest0, ngal
                 allocate(gal(iforest0)%mord(ngal),       stat=ier); call &
                      CheckIOE(ier, trim(c_af)//': gal(:)%mord')
                 allocate(gal(iforest0)%mpi(ngal),        stat=ier); call &
                      CheckIOE(ier, trim(c_af)//': gal(:)%mpi')
                 allocate(gal(iforest0)%rdomi(ngal),      stat=ier); call &
                      CheckIOE(ier, trim(c_af)//': gal(:)%rdomi')
                 allocate(gal(iforest0)%magd(ngal,Nwave), stat=ier); call &
                      CheckIOE(ier, trim(c_af)//': gal(:)%magd')
                 gal(iforest0)%mord(1:ngal)   = mord(1:ngal)
                 gal(iforest0)%mpi(1:ngal)    = mpi(1:ngal)
                 gal(iforest0)%rdomi(1:ngal)  = rdomi(1:ngal)
                 gal(iforest0)%magd(1:ngal,:) = magd(1:ngal,:)
              ENDIF
           ENDIF
        ENDDO
        close(iread); print '(A, I10)', '     ngal_tot = ', ngal_tot

        ! --- reading wind+wpos data and searching corresponding galaxy data
        print '(A)', ' --- reading wind+wpos data & searching corresponding galaxy data'
        write(ci, '(I3.3)') izslice
        fname_wind(irep_slice:irep_slice+2) = ci
        fname_wpos(irep_slice:irep_slice+2) = ci
        nhalo_tot = 0; iz0 = 1; izp1 = iz0 + 1
        DO WHILE(cos%z(izp1) < zslice(izslice+1) .and. izp1 < cos%nz)
           iz0 = iz0 + 1; izp1 = iz0 + 1
        ENDDO ! --> cos%z(iz0) < zslice(izslice+1) <= cos%z(iz0+1)
        DO_forest: DO iforest = 1, N_FOREST
           IF(gal(iforest)%ngal > 0) THEN
              write(ci, '(I3.3)') iforest - 1
              fname_wind(irep_forest:irep_forest+2) = ci
              fname_wpos(irep_forest:irep_forest+2) = ci
              open(iread, file = trim(fname_wind), form='binary', &
                   iostat=ier); call CheckIOE(ier, trim(c_fop)//trim(fname_wind))
              read(iread, iostat=ier) p%nhalo
              ! print '(A, X, I3, I8)', trim(fname_wind), iforest, p%nhalo
              IF(p%nhalo == 0) THEN
                 close(iread)
              ELSE
                 nhalo_tot = nhalo_tot + p%nhalo

                 ! --- reading wind data
                 allocate(p%mpi(p%nhalo), stat=ier); call CheckIOE(ier, &
                      trim(c_af)//': p%mpi'); p%mpi(:) = 0
                 read(iread, iostat=ier) (p%mpi(ihalo), ihalo=1,p%nhalo)
                 close(iread)

                 ! --- reading wpos data
                 allocate(p%x(p%nhalo,3), stat=ier); call CheckIOE(ier, &
                      trim(c_af)//': p%x'); p%x(:,:) = 0.d0
                 open(iread, file = trim(fname_wpos), form='binary', &
                      iostat=ier); call CheckIOE(ier, trim(c_fop)//trim(fname_wpos))
                 read(iread, iostat=ier) ix
                 read(iread, iostat=ier) ((p%x(ihalo,ix), ix=1,3), ihalo=1,p%nhalo)
                 close(iread)

                 ! --- searching for the corresponding galaxy data
                 DO ihalo = 1, p%nhalo
                    drad   = dble(p%x(ihalo, 1)) ! comoving radial distance [Mpc]
                    theta0 = rad2min / drad
                    iz = iz0; izp1 = iz + 1
                    DO WHILE(cos%drad(izp1) < drad .and. izp1 < cos%nz)
                       iz = iz + 1; izp1 = iz + 1
                    ENDDO ! --> cos%drad(iz) < drad <= cos%drad(iz+1)
                    z    = LinearInterp(drad, cos%drad(iz), cos%drad(izp1), &
                                        cos%z(iz), cos%z(izp1))
                           ! redshift of the halo 'ihalo'
                    dmod = LinearInterp(drad, cos%drad(iz), cos%drad(izp1), &
                                        cos%dmod(iz), cos%dmod(izp1))
                    theta_size0 = rad2sec * 1.d-3 / H0 / (drad / (1.d0+z))
                    ! print '(5G15.6)', cos%drad(iz), drad, cos%drad(izp1), z, dmod
                    DO igal = 1, gal(iforest)%ngal
                       IF(gal(iforest)%mpi(igal) == p%mpi(ihalo)) THEN
                          theta_x    = theta0 * dble(p%x(ihalo,2)) ! [arcmin]
                          theta_y    = theta0 * dble(p%x(ihalo,3)) ! [arcmin]
                          theta_size = theta_size0 * gal(iforest)%rdomi(igal) ! [arcsec]
                          !write(iwrite, '(I3, I15, 5X, I2, 2F12.5, 5X, \)') &
                          !     izslice, gal(iforest)%mpi(igal), &
                          !     gal(iforest)%mord(igal), z, drad
                          write(iwrite, '(I2, 2F10.3, 5X, 2E11.2e2, F10.3\)') &
                               gal(iforest)%mord(igal), z, drad, &
                               Mhalo(igal), Mstar(igal), SFR(igal)
                          DO iwave = 4, Nwave
                             write(iwrite, '(F10.3, \)') gal(iforest)%magd(igal,iwave)
                          ENDDO
                          write(iwrite, '(F12.5, 5X, 3F10.3)') &
                               dmod, theta_x, theta_y, theta_size
                       ENDIF
                    ENDDO ! igal = 1, gal(iforest)%ngal
                 ENDDO ! ihalo = 1, p%nhalo
                 deallocate(p%mpi, stat=ier); call CheckIOE(ier, trim(c_daf)//': p%mpi')
                 deallocate(p%x,   stat=ier); call CheckIOE(ier, trim(c_daf)//': p%x')
              ENDIF
           ENDIF
        ENDDO DO_forest
        print '(A, 2I10)', '     nhalo_tot = ', nhalo_tot

        DO iforest = 1, N_FOREST
           IF(gal(iforest)%ngal > 0) THEN
              deallocate(gal(iforest)%mord,  stat=ier); call CheckIOE(ier, &
                   trim(c_daf)//': gal(:)%mord')
              deallocate(gal(iforest)%mpi,   stat=ier); call CheckIOE(ier, &
                   trim(c_daf)//': gal(:)%mpi')
              deallocate(gal(iforest)%magd,   stat=ier); call CheckIOE(ier, &
                   trim(c_daf)//': gal(:)%magd')
              deallocate(gal(iforest)%rdomi, stat=ier); call CheckIOE(ier, &
                   trim(c_daf)//': gal(:)%rdomi')
           ENDIF
        ENDDO
        deallocate(gal, stat=ier); call CheckIOE(ier, trim(c_daf)//': gal')
     ENDIF
  ENDDO DO_zslice
  close(iwrite)

  deallocate(cos%z,    stat=ier); call CheckIOE(ier, trim(c_daf)//': cos%z')
  deallocate(cos%drad, stat=ier); call CheckIOE(ier, trim(c_daf)//': cos%drad')
  deallocate(cos%dlum, stat=ier); call CheckIOE(ier, trim(c_daf)//': cos%dlum')
  deallocate(cos%dmod, stat=ier); call CheckIOE(ier, trim(c_daf)//': cos%dmod')
END PROGRAM count
!!$============================================================================
SUBROUTINE CheckIOE(ier, cerr)
  implicit none
  INTEGER, INTENT(IN) :: ier
  CHARACTER(LEN=*), INTENT(IN) :: cerr

  IF(ier /= 0) THEN
     print '(A, I4)', trim(cerr)//', stat = ', ier
     stop
  ENDIF
END SUBROUTINE CheckIOE
!!$============================================================================
DOUBLE PRECISION FUNCTION LinearInterp(x, x1, x2, y1, y2) RESULT(y)
  implicit none
  DOUBLE PRECISION, INTENT(IN) :: x, x1, x2, y1, y2

  y = (y1 - y2) / (x1 - x2) * (x - x1) + y1
END FUNCTION LinearInterp
!!$============================================================================
