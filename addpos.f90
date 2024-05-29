!!$ ifort -o addpos.out addpos.f90
  PROGRAM addpos
    implicit none

    INTEGER, PARAMETER :: TNWAVE = 18 !30   ! NWAVE*2
!   INTEGER, PARAMETER :: N_FILE = 2
    INTEGER, PARAMETER :: N_FILE = 1 !Makiya
    INTEGER :: ifile

    INTEGER :: i_in, i_out, ier
    CHARACTER(LEN=10) :: cz
    CHARACTER(LEN=50) :: fbase = '/home/makiya/nugc/data/'
!   CHARACTER(LEN=50) :: fend(N_FILE) = (/'400Mpc_LyAbright', '400Mpc_UVbright'/)
    CHARACTER(LEN=50) :: fend(N_FILE) = (/'_UVbright'/) !Makiya
    CHARACTER(LEN=500) :: fname_in, fname_out

    INTEGER, PARAMETER :: N_FOREST = 512 !8 for 100Mpc, 512 for 400Mpc box !Makiya
    INTEGER :: iforest, iforest0, irep, iread = 20, ihalo, ix, nhalo_max
    CHARACTER(LEN=3)  :: ci
    CHARACTER(LEN=10) :: czout
!!$    CHARACTER(LEN=100) :: fbase_pos = '/home/nugc/ngc_400Mpc_2048/position/'
!!$    CHARACTER(LEN=100) :: fbase_pos = '/home/nugc/position/400Mpc_2048_MPOS/'
!    CHARACTER(LEN=100) :: fbase_pos = '/work/100Mpc_512_MPOS2/'
    CHARACTER(LEN=100) :: fbase_pos = '/work/400Mpc_2048_MPOS/'
    CHARACTER(LEN=1000) :: fname_mind, fname_mpos
    TYPE position
       INTEGER :: nhalo
       INTEGER(KIND=8), DIMENSION(:), ALLOCATABLE :: mpi ! mpi(1:nhalo)
!!$       DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: x ! position, x(1:nhalo, 1:3)
       REAL, DIMENSION(:,:), ALLOCATABLE :: x ! position, x(1:nhalo, 1:3)
    END type position
    TYPE(position), ALLOCATABLE :: p(:) ! p(1:N_FOREST)
    INTEGER :: ihalo_cum, nhalo_tot

!   INTEGER, PARAMETER :: N_COLUMN = 53 ! # of total column - 8 (integer)
    INTEGER, PARAMETER :: N_COLUMN = 54 !53 ! # of total column - 8 (integer) !Makiya
    INTEGER, PARAMETER :: N_COLUMN_I = 8 ! # of integer column
    INTEGER :: ic, ic0, ic1, iloop, flag, igal_tot, igal_fail
    INTEGER(KIND=8) :: mpi, mpi_pre, mpi_pre_f
    CHARACTER(LEN=10000) :: buf
    INTEGER(KIND=8), DIMENSION(:) :: ixin(N_COLUMN_I)
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: xin ! xin(1:N_COLUMN)

    call getarg(1, cz)
    IF(cz == '0') THEN
       czout = '018'; cz = 'z0p0' !Makiya
!   IF(cz == '3') THEN
    ELSEIF(cz == '3') THEN
       czout = 'z2.947'; cz = 'z2p9'
    ELSEIF(cz == '4') THEN
       czout = 'z4.040'; cz = 'z4p0'
    ELSEIF(cz == '5') THEN
       czout = 'z4.867'; cz = 'z4p8'
    ELSEIF(cz == '6') THEN
       czout = 'z5.966'; cz = 'z6p0'
    ELSEIF(cz == '7') THEN
       czout = 'z6.976'; cz = 'z7p0'
    ELSE
       print '(A)', '# Select an integer from 3 to 7!!!'; stop
    ENDIF
    allocate(p(N_FOREST), stat=ier); call CheckIerror(ier, &
         '# addpos: allocation fault --- p')
!    fname_mind = trim(fbase_pos)//trim(czout)//'.00000000'
    fname_mind = trim(fbase_pos)//'400Mpc_2048'//'.00000000.'//trim(czout)  !Makiya
    fname_mpos = trim(fname_mind)//'.mpos'
    fname_mind = trim(fname_mind)//'.mind'
!   irep = index(fname_mind, '.mind') - 3; nhalo_tot = 0; nhalo_max = 0
    irep = index(fname_mind, '.mind') - 7; nhalo_tot = 0; nhalo_max = 0 !Makiya

    !全ての iforest の mpos and mind を読み込む
    DO iforest = 1, N_FOREST
       write(ci, '(I3.3)') iforest - 1   !substitution iforest-1 for ci 
       fname_mind(irep:irep+2) = ci; fname_mpos(irep:irep+2) = ci

       ! --- reading MPI from "fname_mind"
       open(iread, file=trim(fname_mind), form='binary', iostat=ier)
       call CheckIerror(ier, '# ReadHaloNumb: fail to open file = '//trim(fname_mind))
       read(iread, iostat=ier) p(iforest)%nhalo
!!$       print '(A, I3, A, I16)', 'p(',iforest,')%nhalo = ', p(iforest)%nhalo
       IF(p(iforest)%nhalo > nhalo_max) nhalo_max = p(iforest)%nhalo
       allocate(p(iforest)%mpi(p(iforest)%nhalo), stat=ier); call &
            CheckIerror(ier, '# (addpos)allocation fault: p(iforest)%mpi')
       p(iforest)%mpi(:) = 0
       read(iread, iostat=ier) (p(iforest)%mpi(ihalo), ihalo=1,p(iforest)%nhalo)
!!$       DO ihalo = 1, p(iforest)%nhalo
!!$          print '(A, I3, A, I6, A, I16)', &
!!$               '--- p(',iforest,')%mpi(',ihalo,')= ', p(iforest)%mpi(ihalo)
!!$       ENDDO
       close(iread)

       ! --- reading position from "fname_mpos"
       allocate(p(iforest)%x(3,p(iforest)%nhalo), stat=ier); call &
            CheckIerror(ier, '# (addpos)allocation fault: p(iforest)%x')
       p(iforest)%x(:,:) = 0.d0
       open(iread, file=trim(fname_mpos), form='binary', iostat=ier); call &
            CheckIerror(ier, '# (addpos)cannot open file='//trim(fname_mpos))
       read(iread, iostat=ier) ix
       read(iread, iostat=ier) ((p(iforest)%x(ix, ihalo), ix=1,3), &
                                ihalo=1,p(iforest)%nhalo)
       close(iread)

       nhalo_tot = nhalo_tot + p(iforest)%nhalo
    ENDDO
    print '(2(A,I10))', '# nhalo_tot = ', nhalo_tot, ', nhalo_max = ', nhalo_max


    i_in = 20; i_out = 21
    DO ifile = 1, N_FILE
!      fname_in = trim(fbase)//trim(cz)//'LAE'//trim(fend(ifile))
       fname_in = trim(fbase)//'kobayashi'//trim(fend(ifile))  !Makiya
!       fname_in = trim(fbase)//'z2p9LAE400Mpc'//trim(fend(ifile))  !Makiya
       fname_out = trim(fname_in)//'_pos.dat'
       fname_in  = trim(fname_in)//'.dat'
       open(i_in, file=trim(fname_in), status='old', iostat=ier); call &
            CheckIerror(ier, '# (addpos)fail to open file='//trim(fname_in))
       open(i_out, file=trim(fname_out), iostat=ier); call &
            CheckIerror(ier, '# (addpos)fail to open file='//trim(fname_out))

       allocate(xin(N_COLUMN), stat=ier); call CheckIerror(ier, &
            '# (addpos)allocation fault: xin')

       ier = 0; iloop = 1; ic0 = N_COLUMN - TNWAVE - 3; mpi_pre = 0; mpi_pre_f = 0
       igal_tot = 0; igal_fail = 0; iforest = 1
       DO WHILE(ier == 0)
          read(i_in, '(A)', iostat=ier) buf
!         IF(iloop <= 5) THEN   
          IF(iloop <= 6) THEN   !writing a header !Makiya
!            IF(iloop == 5) THEN
             IF(iloop == 6) THEN !Makiya
                write(i_out, '(2(A,I2), A)') trim(buf)//' (', &
                     N_COLUMN+N_COLUMN_I+1, '-', N_COLUMN+N_COLUMN_I+3, ')position(x,y,z)'
             ELSE
                write(i_out, '(A)') trim(buf)
             ENDIF
             iloop = iloop + 1
          ELSE
             IF(ier == 0) THEN
                read(buf, *) (ixin(ic), ic=1,2), &
                    (xin(ic), ic=1,ic0), (ixin(ic), ic=3,N_COLUMN_I), &
                    (xin(ic), ic=ic0+1,N_COLUMN)
!                     (xin(ic), ic=1,20), (ixin(ic), ic=3,N_COLUMN_I), & !Makiya
!                     (xin(ic), ic=20+1,N_COLUMN)
!               mpi = ixin(4) ! mpi(gal)
                mpi = ixin(5) ! mpi(gal) !Makiya
                iforest = ixin(4)+1!Makiya
                IF(mpi /= mpi_pre) THEN
                   igal_tot = igal_tot + 1; mpi_pre = mpi
                   flag = 0

!                   DO_iforest: DO WHILE(flag == 0)
!                      DO ihalo = 1, p(iforest)%nhalo
!                         IF(mpi == p(iforest)%mpi(ihalo)) THEN
!                            flag = 1
!                            write(i_out, '(A, 3E13.4)') trim(buf), &
!                                 (p(iforest)%x(ix,ihalo), ix=1,3)
!                            goto 100
!                         ENDIF
!                      ENDDO
!                      IF(iforest < N_FOREST) iforest = iforest + 1
!                   ENDDO DO_iforest

                   DO ihalo = 1, p(iforest)%nhalo
                      IF(mpi == p(iforest)%mpi(ihalo)) THEN
                         flag = 1
                         write(i_out, '(A, 3E13.4)') trim(buf), &
                              (p(iforest)%x(ix,ihalo), ix=1,3)
                      ENDIF
                   ENDDO
                   goto 100

!!$                   print '(A, 2I10)', '1:', iforest, ihalo
                ELSE
!!$                   print '(A, 2I10)', '2:', iforest, ihalo
!                   write(i_out, '(A, I5, 3E13.4)') trim(buf), &
!                        iforest-1, (p(iforest)%x(ix,ihalo), ix=1,3) !comment out by Makiya
                   write(i_out, '(A, I5, 3E13.4)') trim(buf), &
                        (p(iforest)%x(ix,ihalo), ix=1,3) 
                ENDIF
100             IF(flag == 0 .and. mpi /= mpi_pre_f) THEN
                   print '(I10, I16)', ixin(1), mpi
                   mpi_pre_f = mpi; igal_fail = igal_fail + 1
                ENDIF
             ENDIF
         ENDIF
       ENDDO
       print '(3I10)', ixin(1), ixin(2), ixin(3)
!      print '(A, G8.3, A)', '#('//trim(cz)//') fraction of the failed galaxies:', & !Makiya
!            100.d0 * dble(igal_fail)/dble(igal_tot), '[%]'
       print '(A, G10.3, A)', '#('//trim(cz)//') fraction of the failed galaxies:', &
            100.d0 * dble(igal_fail)/dble(igal_tot), '[%]'
       close(i_in); close(i_out)
       deallocate(xin, stat=ier); call CheckIerror(ier, &
            '# (addpos)dallocation fault: xin')
    ENDDO
!!$============================================================================
  CONTAINS
!!$============================================================================
    SUBROUTINE CheckIerror(ier, cerr)
      implicit none
      INTEGER, INTENT(IN) :: ier
      CHARACTER(LEN=*), INTENT(IN) :: cerr

      IF(ier /= 0) THEN
         print '(A, I4)', trim(cerr)//', stat = ', ier; stop
      ENDIF
    END SUBROUTINE CheckIerror
!!$============================================================================
  END PROGRAM addpos
