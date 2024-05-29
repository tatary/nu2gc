  PROGRAM addpos
    implicit none

    INTEGER, PARAMETER :: N_WAVE = 2

    INTEGER :: iwave

    INTEGER :: i_in, i_out, ier, iz
    CHARACTER(LEN=10) :: cz
    CHARACTER(LEN=50) :: fbase = 'data/LAE/2013Sep/'
    CHARACTER(LEN=50) :: fend(N_WAVE) = (/'400Mpc_LyAbright', '400Mpc_UVbright'/)
    CHARACTER(LEN=500) :: fname_in, fname_out
    CHARACTER(LEN=1000) :: header, buf

    INTEGER, PARAMETER :: N_FOREST = 512
    INTEGER :: iforest, irep
    CHARACTER(LEN=3)  :: ci
    CHARACTER(LEN=10) :: czout
    CHARACTER(LEN=100) :: fbase_pos = 'position/400Mpc_2048_MPOS/'
    CHARACTER(LEN=1000) :: fname_mind, fname_mpos
    TYPE position
       INTEGER :: nhalo
       INTEGER, DIMENSION(:), ALLOCATABLE :: mpi ! mpi(1:nhalo)
       DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: x ! position, x(1:nhalo, 1:3)
    END type position
    TYPE(position), ALLOCATABLE :: pin(:) ! pin(1:N_FOREST)

    ! --- functions
    INTEGER :: ReadHaloNumb


    call getarg(1, cz)
    IF(cz == '3') THEN
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
    allocate(pin(N_FOREST), stat=ier); call CheckIerror(ier, &
         '# addpos: allocation fault --- pin')
    fname_mind = trim(fbase_pos)//trim(czout)//'.00000000'
    fname_mpos = trim(fname_mind)//'.mpos'
    fname_mind = trim(fname_mind)//'.mind'
    irep = index(fname_mind, '.mind') - 3
    DO iforest = 1, N_FOREST
       write(ci, '(I3.3)') iforest - 1
       fname_mind(irep:irep+2) = ci; fname_mpos(irep:irep+2) = ci
       pin(iforest)%nhalo = ReadHaloNumb(fname_mind)
       print '(A, I3, A, I10)', 'pin(',iforest,')%nhalo = ', pin(iforest)%nhalo
       allocate(pin(iforest)%mpi(pin(iforest)%nhalo), stat=ier); call &
            CheckIerror(ier, '# addpos: allocation fault --- pin(iforest)%mpi')
!!$       allocate(pin(iforest)%x(pin(iforest)%nhalo, 3), stat=ier); call &
!!$            CheckIerror(ier, '# addpos: allocation fault --- pin(iforest)%x')
       allocate(pin(iforest)%x(3,pin(iforest)%nhalo), stat=ier); call &
            CheckIerror(ier, '# addpos: allocation fault --- pin(iforest)%x')
       call ReadMPIandPosition(pin(iforest)%nhalo, pin(iforest)%mpi, pin(iforest)%x, &
                               fname_mind, fname_mpos)
    ENDDO

    i_in = 20; i_out = 21
    DO iwave = 1, N_WAVE
       fname_in = trim(fbase)//trim(cz)//'LAE'//trim(fend(iwave))
       fname_out = trim(fname_in)//'_pos.dat'
       fname_in  = trim(fname_in)//'.dat'
       open(i_in, file=trim(fname_in), status='old', iostat=ier)
       call CheckIerror(ier, '# addpos: fail to open file='//trim(fname_in))
       open(i_out, file=trim(fname_out), iostat=ier)
       call CheckIerror(ier, '# addpos: fail to open file='//trim(fname_out))

       call SearchAndWritePosition(iwave)

       close(i_in); close(i_out)
    ENDDO

    DO iforest = 1, N_FOREST
       deallocate(pin(iforest)%mpi, stat=ier); call &
            CheckIerror(ier, '# addpos: deallocation fault --- pin(iforest)%mpi')
       deallocate(pin(iforest)%x, stat=ier); call &
            CheckIerror(ier, '# addpos: deallocation fault --- pin(iforest)%x')
    ENDDO
    deallocate(pin, stat=ier); call CheckIerror(ier, &
         '# addpos: deallocation fault --- pin')
  CONTAINS
!!$============================================================================
    SUBROUTINE SearchAndWritePosition(iwave)
      INTEGER, INTENT(IN) :: iwave
      INTEGER, PARAMETER :: nc(N_WAVE) = (/56, 58/)
      INTEGER :: ier, ic, iforest, ihalo, ix
      CHARACTER(LEN=5000) :: buf
      CHARACTER(LEN=100) :: cerr(2) = (/'# (SearchAndWritePosition)Allocation Fault:', &
           '# (SearchAndWritePosition)Deallocation Fault:'/)
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: xin ! xin(1:nc)
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: px  ! px(1:3)
      DOUBLE PRECISION :: mpi

      allocate(xin(nc(iwave)), stat=ier); call CheckIerror(ier, trim(cerr(1))//' xin')
      allocate(px(3), stat=ier); call CheckIerror(ier, trim(cerr(1))//' px')

      ier = 0
      DO WHILE(ier == 0)
         read(i_in, *, iostat=ier) buf
         IF(ier == 0) THEN
            IF(buf(1:1) == '#') THEN ! header
               write(i_out, '(A)') trim(buf)
            ELSE ! data
               read(buf, *) (xin(ic), ic=1, nc(iwave))
               IF(iwave == 1) THEN
                  mpi = xin(35)
               ELSEIF(iwave == 2) THEN
                  mpi = xin(37)
               ENDIF
               DO iforest = 1, N_FOREST
                  DO ihalo = 1, pin(iforest)%nhalo
                     IF(mpi == pin(iforest)%mpi(ihalo)) THEN
                        write(i_out, '(A, 3G12.5)') &
                             trim(buf), (pin(iforest)%x(ihalo,ix), ix=1,3)
                        exit
                     ENDIF
                  ENDDO
               ENDDO
            ENDIF
         ENDIF
      ENDDO

      deallocate(xin, stat=ier); call CheckIerror(ier, trim(cerr(2))//' xin')
      deallocate(px, stat=ier); call CheckIerror(ier, trim(cerr(2))//' px')
    END SUBROUTINE SearchAndWritePosition
!!$============================================================================
  END PROGRAM addpos
!!$============================================================================
!!$============================================================================
  INTEGER FUNCTION ReadHaloNumb(fname_mind) RESULT(nhalo)
    implicit none
    CHARACTER(LEN=*), INTENT(IN) :: fname_mind
    INTEGER :: i_in = 20, ier

    open(i_in, file=trim(fname_mind), status='old', form='binary', iostat=ier)
    call CheckIerror(ier, '# ReadHaloNumb: fail to open file = '//trim(fname_mind))
    read(i_in) nhalo
    close(i_in)
  END FUNCTION ReadHaloNumb
!!$============================================================================
  SUBROUTINE ReadMPIandPosition(nhalo, mpi, x, fname_mind, fname_mpos)
    implicit none
    INTEGER, INTENT(IN) :: nhalo
    CHARACTER(LEN=*), INTENT(IN) :: fname_mind, fname_mpos
    INTEGER, INTENT(INOUT) :: mpi(nhalo)
!!$    DOUBLE PRECISION, INTENT(INOUT) :: x(nhalo,3)
    DOUBLE PRECISION, INTENT(INOUT) :: x(3,nhalo)
    INTEGER :: iread = 20, ier, i_mind, i_mpos, ihalo, ix
    CHARACTER(LEN=100) :: cerr = '# ReadMPIandPosition: fail to open file = '

    open(iread, file=trim(fname_mind), form='binary', iostat=ier); call &
         CheckIerror(ier, trim(cerr)//trim(fname_mind))
    read(iread) i_mind
    IF(i_mind /= nhalo) &
         print '(2(A,I10))', '# Different!: N_HALO=', nhalo, 'N_MIND=', i_mind
    read(iread) (mpi(ihalo), ihalo=1,i_mind)
    close(iread)

    open(iread, file=trim(fname_mpos), form='binary', iostat=ier); call &
         CheckIerror(ier, trim(cerr)//trim(fname_mpos))
    read(iread) i_mpos
    IF(i_mpos /= nhalo) &
         print '(2(A,I10))', '# Different!: N_HALO=', nhalo, ', N_MPOS=', i_mpos
!!$    read(iread) ((x(ihalo,ix), ix=1,3), ihalo=1,i_mpos)
    read(iread) ((x(ix,ihalo), ix=1,3), ihalo=1,i_mpos)
    close(iread)

    return
  END SUBROUTINE ReadMPIandPosition
!!$============================================================================
  SUBROUTINE CheckIerror(ier, cerr)
    implicit none
    INTEGER, INTENT(IN) :: ier
    CHARACTER(LEN=*) :: cerr

    IF(ier /= 0) THEN
       print '(A, I4)', trim(cerr)//', stat = ', ier; stop
    ENDIF
  END SUBROUTINE CheckIerror
!!$============================================================================

