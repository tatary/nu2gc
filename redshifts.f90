PROGRAM main

  implicit none

  INTEGER :: num_nbody
!  INTEGER, INTENT(INOUT) :: run, num_nbody, type_ssp
  DOUBLE PRECISION :: z_nbody(999)

  INTEGER, PARAMETER :: NUM_NBODY_MAX = 999
  INTEGER :: i, j
  CHARACTER :: num*50
  CHARACTER :: base*100 = '/home/nugc/'
  DOUBLE PRECISION :: scale_factor(NUM_NBODY_MAX), x

!  open(1, file = '/home/nugc/'//'ngc_100Mpc_512'//'/'&
!       //'ngc_100Mpc_512'//'.000.mrgp', status = 'old')
!  open(1, file = '/work/100Mpc_512_MRGT2/100Mpc_512.00000000.mrgp', status = 'old')
  open(1, file = '/home/nugc/ngc_400Mpc_2048/400Mpc_2048_MRGT/400Mpc_2048.00000000.mrgp', status = 'old')
  DO i = 1, NUM_NBODY_MAX
     read(1, *, end = 1001) x, scale_factor(i), j, j
  ENDDO
1001 close(1)

  num_nbody = i - 1
  DO i = 1, num_nbody
     z_nbody(i) = 1.d0 / scale_factor(i) - 1.d0
  ENDDO
  
  DO i = 1, num_nbody
     print *, i, z_nbody(i)
  ENDDO

END PROGRAM main

