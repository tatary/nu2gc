! ifort -o makeHaloMF.out makeHaloMF.f90
PROGRAM makeHaloMF
  implicit none
  INTEGER, PARAMETER:: N = 25000000, NN = 200, NNN = 200000
  INTEGER :: n_forest, iforest
  CHARACTER(LEN=50)  :: fbase_nbody = '/work/'
!!$    === mgrt ===
  INTEGER :: num_tot
  INTEGER :: f_des(0:N), n_des(0:N), f_prg(0:N), numb(0:N)
  INTEGER :: hori(0:N)
  INTEGER(KIND=8) :: mpi(0:N)
  DOUBLE PRECISION :: mhalo(0:N)
!!$    === mgrp ===
  INTEGER :: num_now(NN), st_halo(NN), num_step
  DOUBLE PRECISION :: times(NN), zplus1ar(NN)

  INTEGER, PARAMETER :: Nbin = 100
  TYPE DistributionFunction
     INTEGER :: ibin_max, ibin_min, flag
     INTEGER :: nprog_min(Nbin), nprog_max(Nbin)
     DOUBLE PRECISION :: z ! redshift
     DOUBLE PRECISION :: t_univ ! [Myr]
     DOUBLE PRECISION :: base, step, invstep
     DOUBLE PRECISION :: x(Nbin)
     DOUBLE PRECISION :: n(Nbin), ncum(Nbin), nprog(Nbin)
  END type DistributionFunction
  TYPE(DistributionFunction), ALLOCATABLE :: mf(:) ! mf(1:num_step)
  INTEGER :: itime, ibin, ihalo, me, itime2
  DOUBLE PRECISION :: zmin, zmax, tmin, tmax
  DOUBLE PRECISION, PARAMETER :: EPS = 1.d-10

  INTEGER   :: i, j, ii, iread, iwrite, ier, iprint
  CHARACTER :: ci*5, cz*3
  CHARACTER(LEN=100) :: fname_nbody, fname_mrgt, fname_mrgp, fname_mf
  CHARACTER(LEN=1000) :: buf
  DOUBLE PRECISION :: x, h0, OMEGA_M, OMEGA_b, OMEGA_L, sum
  DOUBLE PRECISION :: Cube ! function


  call getarg(1, fname_nbody)
  IF(index(fname_nbody, '1120') == 0) THEN
     fbase_nbody = trim(fbase_nbody)//trim(fname_nbody)//'/'
  ELSE
     fbase_nbody = trim(fbase_nbody(1:len_trim(fbase_nbody)-1))//'2/'//&
          trim(fname_nbody)//'/'
  ENDIF
  iread = 1
  open(iread, file = trim(fbase_nbody)//'param', status = 'old', iostat=ier); call &
       CheckIerr(ier, '# fail to open file: '//trim(fbase_nbody)//'param')
  read(iread, *) n_forest
  read(iread, *) x ! simulation box size [Mpc/h]
  read(iread, *) OMEGA_M
  read(iread, *) h0
  read(iread, *) OMEGA_b
  close(iread)
  sum = 1.d0 / Cube(x) ! [h^3/Mpc^3]
  OMEGA_L = 1.d0 - OMEGA_M ! Omega_L


  ii = index(fname_nbody, '_MRGT')
  fname_mrgt = trim(fbase_nbody)//trim(fname_nbody(1:ii-1))//'.00000000.mrgt'
  IF(n_forest >= 10000) THEN
     ii = index(fname_mrgt, '00000.mrgt')
  ELSEIF(n_forest >= 1000) THEN
     ii = index(fname_mrgt, '0000.mrgt')
  ELSE
     ii = index(fname_mrgt, '000.mrgt')
  ENDIF
  fname_mrgp = fname_mrgt(1:len_trim(fname_mrgt)-1)//'p'; iprint = 0
  DO iforest = 0, n_forest-1
     IF(n_forest >= 10000) THEN
        write(ci, '(I5.5)') iforest
     ELSEIF(n_forest >= 1000) THEN
        write(ci, '(I4.4)') iforest
     ELSE
        write(ci, '(I3.3)') iforest
     ENDIF

     ! normalize
     f_des(:) = 0; n_des(:) = 0; f_prg(:) = 0; numb(:) = 0; hori(:) = 0; mpi(:) = 0
     mhalo(:) = 0.d0; times(:) = 0.d0; num_now(:) = 0; st_halo(:) = 0
     zplus1ar(:) = 0.d0

     IF(n_forest >= 10000) THEN
        fname_mrgt(ii:ii+4) = ci; fname_mrgp(ii:ii+4) = ci
     ELSEIF(n_forest >= 1000) THEN
        fname_mrgt(ii:ii+3) = ci; fname_mrgp(ii:ii+3) = ci
     ELSE
        fname_mrgt(ii:ii+2) = ci; fname_mrgp(ii:ii+2) = ci
     ENDIF
     call PrintMessage(iprint, n_forest, iforest, fname_mrgt)


     ! --- reading the data from mrgp file
     iread = 1; j = 1
     open(iread, file = trim(fname_mrgp), status = 'old', iostat = ier); call &
         CheckIerr(ier, '# fail to open file: '//trim(fname_mrgp))
     ier = 0; x = 0.d0
     DO WHILE(ier == 0 .and. x <= 1.d0)
        read(iread, '(A)', iostat=ier) buf
        IF(ier == 0) THEN
           read(buf, *) times(j), x, num_now(j), st_halo(j)
                        ! times(:) : [Myr]
           zplus1ar(j) = 1.d0 / x
           j = j + 1
        ENDIF
     ENDDO
     close(iread)
     num_step = j - 1
     ! times(:) = times(:) * 1.d+3 ! [Gyr] --> [Myr]


     ! --- allocate the structure for mass function
     IF(iforest == 0) THEN
        allocate(mf(num_step), stat=ier); call CheckIerr(ier, '# Allocation Fault: mf')
        zmin = 100.d0; zmax = -100.d0
        DO itime = 1, num_step
           mf(itime)%flag = 0

           mf(itime)%nprog_min(1:Nbin) = 10000; mf(itime)%nprog_max(1:Nbin) = 1
           mf(itime)%n(1:Nbin) = 0.d0; mf(itime)%ncum(1:Nbin) = 0.d0
           mf(itime)%nprog(1:Nbin) = 0.d0

           mf(itime)%z      = zplus1ar(itime) - 1.d0
           mf(itime)%t_univ = times(itime) ! [Myr]
           IF(zmin > mf(itime)%z) THEN
              zmin = mf(itime)%z; tmax = mf(itime)%t_univ
           ENDIF
           IF(zmax < mf(itime)%z) THEN
              zmax = mf(itime)%z; tmin = mf(itime)%t_univ
           ENDIF

           mf(itime)%base = 0.d0; mf(itime)%step = 0.2d0
           mf(itime)%invstep = 1.d0 / mf(itime)%step
           DO ibin = 1, Nbin
              mf(itime)%x(ibin) = mf(itime)%base + mf(itime)%step * dble(ibin)
           ENDDO
        ENDDO
     ENDIF


     ! --- reading the data from mrgt file
     iread = 1
     open(iread, file = trim(fname_mrgt), status = 'old', iostat = ier); call &
         CheckIerr(ier, '# fail to open file: '//trim(fname_mrgt))
     ier = 0
     DO WHILE(ier == 0)
        read(iread, '(A)', iostat=ier) buf
        IF(ier == 0) THEN
           read(buf, *) i, f_des(i), n_des(i), f_prg(i), numb(i), &
                hori(i), mpi(i), mhalo(i)
           mhalo(i) = mhalo(i) * h0 ! [Msun] --> [Msun/h]
        ENDIF
     ENDDO
     close(iread)
     num_tot = i


     ! --- calculating the halo mass function
     DO itime = 1, num_step
        DO ihalo = 1, num_now(itime)
           me = st_halo(itime) + ihalo - 1 ! (ihalo,itime) --> me
           ibin = int(anint((log10(mhalo(me)) - mf(itime)%base) * mf(itime)%invstep))
           IF(ibin >= 1 .and. ibin <= Nbin) THEN
              IF(mf(itime)%flag == 0) mf(itime)%flag = 1
              mf(itime)%n(ibin)     = mf(itime)%n(ibin) + sum
              mf(itime)%nprog(ibin) = mf(itime)%nprog(ibin) + sum * numb(me)
              IF(mf(itime)%nprog_min(ibin) > numb(me)) &
                   mf(itime)%nprog_min(ibin) = numb(me)
              IF(mf(itime)%nprog_max(ibin) < numb(me)) THEN
                 mf(itime)%nprog_max(ibin) = numb(me)
              ENDIF
           ENDIF
        ENDDO
     ENDDO
  ENDDO


  ! --- calculating the cumulative number density of halo mass function
  DO itime = 1, num_step
     i = 1
     DO ibin = Nbin, 1, -1
        IF(ibin == Nbin) THEN
           mf(itime)%ncum(ibin) = mf(itime)%n(ibin)
        ELSE
           IF(mf(itime)%n(ibin) > 0.d0) THEN
              mf(itime)%ncum(ibin) = mf(itime)%ncum(ibin+1) + mf(itime)%n(ibin)
              i = 2
           ELSE
              IF(i == 1) THEN
                 mf(itime)%ibin_max = ibin
              ELSEIF(i == 2) THEN
                 mf(itime)%ibin_min = ibin; i = 3
              ENDIF
           ENDIF
        ENDIF
     ENDDO
     mf(itime)%ibin_min = mf(itime)%ibin_min + 1
     mf(itime)%ibin_max = mf(itime)%ibin_max - 1
  ENDDO


  ! --- writing the halo mass function
  ii = index(fname_nbody, '_MRGT')
  fname_mf = trim(fbase_nbody)//trim(fname_nbody(1:ii-1))//'_HaloMF.dat'
  iwrite = 1
  open(iwrite, file=trim(fname_mf), status='unknown', iostat=ier); call &
       CheckIerr(ier, '# fail to open file: '//trim(fname_mf))
  write(iwrite, '(4(A, F10.5), A)') '# Halo Mass Functions for '//&
       trim(fname_nbody)//' N-body run at z = ', zmin, '--', zmax, &
       ' (t_univ = ', tmin, '--', tmax, ' Gyr)'
  write(iwrite, '(4(A, F8.3))') '# used cosmological parameters: h = ', &
       h0, ', Omega_b = ', OMEGA_b, ', Omega_M = ', OMEGA_M, ', Omega_L = ', OMEGA_L
  itime2 = 1
  DO itime = 1, num_step ! writing the information of each index
     IF(mf(itime)%flag == 1) THEN
        write(iwrite, '(A, I2, 2(A, F10.5), A)') '# index ', itime2-1, &
             ': z = ', mf(itime)%z, ' (t_univ = ', mf(itime)%t_univ, ' Gyr)'
        itime2 = itime2 + 1
     ENDIF
  ENDDO
  write(iwrite, '(A)') '# ----------------------------------------------------'
  itime2 = 1
  DO itime = 1, num_step ! writing the information of each index
     IF(mf(itime)%flag == 1) THEN
        write(iwrite, '(A, I2, 2(A, F10.5), A)') '# index ', itime2-1, &
          ': z = ', mf(itime)%z, ' (t_univ = ', mf(itime)%t_univ, ' Gyr)'
        write(iwrite, '(A)') '# (1)log10[Mhalo/(Msun/h)] (2)dn/dlog10Mhalo[h^3/'// &
             'Mpc^3/dex] (3)n(>Mhalo)[h^3/Mpc^3] (4)dn_prog/dlog10Mhalo[h^3/Mpc^3'//&
             '/dex] (5)Nhalo (6,7)Nprog_{min,max} per halo'
        ! DO ibin = 1, Nbin
        DO ibin = mf(itime)%ibin_min, mf(itime)%ibin_max
           IF(mf(itime)%n(ibin) < EPS) THEN
              mf(itime)%nprog_min(ibin) = 0; mf(itime)%nprog_max(ibin) = 0
           ENDIF
           IF(int(mf(itime)%n(ibin) / sum) == 1) &
                mf(itime)%nprog_min(ibin) = mf(itime)%nprog_max(ibin)
           IF(mf(itime)%nprog_min(ibin) > mf(itime)%nprog_max(ibin)) &
                mf(itime)%nprog_min(ibin) = mf(itime)%nprog_max(ibin)
           write(iwrite, '(F7.2, 3G13.5, 5X, 3I10)') mf(itime)%x(ibin), &
                mf(itime)%n(ibin) * mf(itime)%invstep, mf(itime)%ncum(ibin), &
                mf(itime)%nprog(ibin) * mf(itime)%invstep, &
                int(mf(itime)%n(ibin) / sum), &
                mf(itime)%nprog_min(ibin), mf(itime)%nprog_max(ibin)
        ENDDO
        write(iwrite, *); write(iwrite, *)
        itime2 = itime2 + 1
     ENDIF
  ENDDO
  close(iwrite)
END PROGRAM makeHaloMF
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
DOUBLE PRECISION FUNCTION Cube(x) RESULT(x3)
  implicit none
  DOUBLE PRECISION, INTENT(IN) :: x

  x3 = x * x * x
END FUNCTION Cube
!!$============================================================================
SUBROUTINE PrintMessage(iprint, n_forest, iforest, fname_mrgt)
  implicit none
  INTEGER, INTENT(INOUT) :: iprint
  INTEGER, INTENT(IN) :: n_forest, iforest
  CHARACTER(LEN=*), INTENT(IN) :: fname_mrgt


  IF(iprint == 0 &
       .or. (iprint == 1 .and. iforest > int(dble(n_forest) * 0.1d0)) &
       .or. (iprint == 2 .and. iforest > int(dble(n_forest) * 0.2d0)) &
       .or. (iprint == 3 .and. iforest > int(dble(n_forest) * 0.3d0)) &
       .or. (iprint == 4 .and. iforest > int(dble(n_forest) * 0.4d0)) &
       .or. (iprint == 5 .and. iforest > int(dble(n_forest) * 0.5d0)) &
       .or. (iprint == 6 .and. iforest > int(dble(n_forest) * 0.6d0)) &
       .or. (iprint == 7 .and. iforest > int(dble(n_forest) * 0.7d0)) &
       .or. (iprint == 8 .and. iforest > int(dble(n_forest) * 0.8d0)) &
       .or. (iprint == 9 .and. iforest > int(dble(n_forest) * 0.9d0))) THEN
     iprint = iprint + 1

     IF(n_forest >= 10000) THEN
        print '(2(A, I5), \)', '(', iforest+1, '/', n_forest
     ELSEIF(n_forest >= 1000) THEN
        print '(2(A, I4), \)', '(', iforest+1, '/', n_forest
     ELSEIF(n_forest >= 100) THEN
        print '(2(A, I3), \)', '(', iforest+1, '/', n_forest
     ELSEIF(n_forest >= 10) THEN
        print '(2(A, I2), \)', '(', iforest+1, '/', n_forest
     ELSE
        print '(2(A, I1), \)', '(', iforest+1, '/', n_forest
     ENDIF
     print '(A, I3, A)', ') reading '//fname_mrgt(1:len_trim(fname_mrgt)-1)//&
          '{t,p} (', int(100.d0 * dble(iforest) / dble(n_forest)), &
          '% finished)'
  ENDIF
END SUBROUTINE PrintMessage
!!$============================================================================
