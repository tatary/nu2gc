! ifort -o ascii2bin.out ascii2bin.f90
! (e.g.) ./ascii2bin.out 70hMpc_512_Planck_MRGT
PROGRAM ascii2bin
  implicit none
  INTEGER, PARAMETER:: N = 25000000, NN = 200, NNN = 200000
  INTEGER :: iforest, iforest_max
!!$    === mgrt ===
  INTEGER :: num_tot
  INTEGER :: f_des(0:N), n_des(0:N), f_prg(0:N), numb(0:N)
  INTEGER :: hori(0:N)
  INTEGER(KIND=8) :: mpi(0:N)
  DOUBLE PRECISION :: mhalo(0:N)
!!$    === mgrp ===
  INTEGER :: num_now(NN), st_halo(NN), num_step
  DOUBLE PRECISION :: times(NN), zplus1ar(NN)

  INTEGER :: i, j, ii, iii, iz, iread, iwrite, ier
  CHARACTER(LEN=100) :: base = '/work/'
  CHARACTER(LEN=100) :: c_nbody
  CHARACTER(LEN=100) :: fname_mrgt, fname_mrgp, fname_bin, fname_z
  CHARACTER(LEN=500) :: buf
  CHARACTER :: ci*5, cz*3
  DOUBLE PRECISION :: x


  call getarg(1, c_nbody)
  IF(index(c_nbody, '1120') == 0) THEN
     base = trim(base)//trim(c_nbody)//'/' ! (e.g.,) base = /work/70hMpc_512_Planck_MRGT/
  ELSE ! for 1120 Mpc/h w/ 8192^3
     base = '/work2/'//trim(c_nbody)//'/' ! base = /work2/1120hMpc_8192_Planck_MRGT/
  ENDIF
  iread = 1
  open(iread, file = trim(base)//'param', status = 'old', iostat=ier); call &
       CheckIerr(ier, '# fail to open file: '//trim(base)//'param')
  read(iread, *) iforest_max
  read(iread, *) x ! simulation box size [Mpc/h]
  read(iread, *) x ! Omega_M
  read(iread, *) x ! h0
  read(iread, *) x ! Omega_b
  close(iread)

  i = index(c_nbody, '_MRGT')
  fname_mrgt = trim(base)//trim(c_nbody(1:i-1))//'.00000000.mrgt'

  IF(iforest_max >= 10000) THEN
     ii = index(fname_mrgt, '00000.mrgt')
  ELSEIF(iforest_max >= 1000) THEN
     ii = index(fname_mrgt, '0000.mrgt')
  ELSE
     ii = index(fname_mrgt, '000.mrgt')
  ENDIF
  fname_mrgp = fname_mrgt(1:len_trim(fname_mrgt)-1)//'p'

  IF(iforest_max >= 10000) THEN
     fname_bin = trim(base)//trim(c_nbody)//'.00000.mrgb'
     iii = index(fname_bin, '00000.mrgb')
  ELSEIF(iforest_max >= 1000) THEN
     fname_bin = trim(base)//trim(c_nbody)//'.0000.mrgb'
     iii = index(fname_bin, '0000.mrgb')
  ELSE
     fname_bin = trim(base)//trim(c_nbody)//'.000.mrgb'
     iii = index(fname_bin, '000.mrgb')
  ENDIF

  DO iforest = 0, iforest_max-1
     IF(iforest_max >= 10000) THEN
        write(ci, '(I5.5)') iforest
     ELSEIF(iforest_max >= 1000) THEN
        write(ci, '(I4.4)') iforest
     ELSE
        write(ci, '(I3.3)') iforest
     ENDIF

     ! normalize
     f_des(:) = 0; n_des(:) = 0; f_prg(:) = 0; numb(:) = 0; hori(:) = 0; mpi(:) = 0
     mhalo(:) = 0.d0; times(:) = 0.d0; num_now(:) = 0; st_halo(:) = 0
     zplus1ar(:) = 0.d0

     IF(iforest_max >= 10000) THEN
        fname_mrgt(ii:ii+4) = ci; fname_mrgp(ii:ii+4) = ci; fname_bin(iii:iii+4) = ci
        print '(2(A, I5), A)', '(', iforest+1, '/', iforest_max, &
             ')'//trim(fname_mrgt)//', '//trim(fname_mrgp)
     ELSEIF(iforest_max >= 1000) THEN
        fname_mrgt(ii:ii+3) = ci; fname_mrgp(ii:ii+3) = ci; fname_bin(iii:iii+3) = ci
        print '(2(A, I4), A)', '(', iforest+1, '/', iforest_max, &
             ')'//trim(fname_mrgt)//', '//trim(fname_mrgp)
     ELSEIF(iforest_max >= 100) THEN
        fname_mrgt(ii:ii+2) = ci; fname_mrgp(ii:ii+2) = ci; fname_bin(iii:iii+3) = ci
        print '(2(A, I3), A)', '(', iforest+1, '/', iforest_max, &
             ')'//trim(fname_mrgt)//', '//trim(fname_mrgp)
     ELSEIF(iforest_max >= 10) THEN
        fname_mrgt(ii:ii+2) = ci; fname_mrgp(ii:ii+2) = ci; fname_bin(iii:iii+3) = ci
        print '(2(A, I2), A)', '(', iforest+1, '/', iforest_max, &
             ')'//trim(fname_mrgt)//', '//trim(fname_mrgp)
     ELSE
        fname_mrgt(ii:ii+2) = ci; fname_mrgp(ii:ii+2) = ci; fname_bin(iii:iii+2) = ci
        print '(2(A, I1), A)', '(', iforest+1, '/', iforest_max, &
             ') '//trim(fname_mrgt)//' + '//trim(fname_mrgp)
     ENDIF


     ! --- reading mrgp file
!!$     j = 1
!!$     open(12, file = file_mrgp, status = 'old')
!!$30   read(12, *, end = 40) times(j), x, num_now(j), st_halo(j)
!!$     ! times(:) : [Myr]
!!$     if(x > 1.d0) goto 40
!!$     zplus1ar(j) = 1.d0 / x
!!$     j = j + 1
!!$     goto 30
!!$40   continue
!!$     close(12)
     iread = 1
     open(iread, file = trim(fname_mrgp), status = 'old', iostat = ier); call &
          CheckIerr(ier, '# fail to open file: '//trim(fname_mrgp))
     ier = 0; x = 0.d0; j = 1
     DO WHILE(ier == 0 .and. x <= 1.d0)
        read(iread, '(A)', iostat=ier) buf
        IF(ier == 0) THEN
           read(buf, *) times(j), x, num_now(j), st_halo(j)
           IF(x <= 1.d0) THEN
              zplus1ar(j) = 1.d0 / x
              j = j + 1
           ENDIF
        ENDIF
     ENDDO
     close(iread)
     num_step = j - 1
     times(:) = times(:) * 1.d+3 ! [Gyr] --> [Myr]


     ! --- reading mrgt file
!!$     open(11, file = fname_mrgt, status = 'old')
!!$10   read(11, *, end=20) i, f_des(i), n_des(i), f_prg(i), numb(i),&
!!$          hori(i), mpi(i), mhalo(i)
!!$     goto 10
!!$20   continue
!!$     close(11)
     iread = 1
     open(iread, file = trim(fname_mrgt), status = 'old', iostat = ier); call &
         CheckIerr(ier, '# fail to open file: '//trim(fname_mrgt))
     ier = 0
     DO WHILE(ier == 0)
        read(iread, '(A)', iostat=ier) buf
        IF(ier == 0) THEN
           read(buf, *) i, f_des(i), n_des(i), f_prg(i), numb(i), &
                hori(i), mpi(i), mhalo(i)
        ENDIF
     ENDDO
     close(iread)
     num_tot = i


     ! --- writing the N-body data into mrgb file
     print '(A)', '     --> fname_bin = '//trim(fname_bin)
     iwrite = 1
     open(iwrite, file = trim(fname_bin), form = 'unformatted', iostat = ier); call &
         CheckIerr(ier, '# fail to open file: '//trim(fname_bin))
     write(iwrite) num_step, num_tot
     print '(A, I3, A, I8)', '     num_step = ', num_step, ', num_tot = ', num_tot

     ! from mrgp
     write(iwrite) (times(i),    i = 1, num_step)
     write(iwrite) (zplus1ar(i), i = 1, num_step)
     write(iwrite) (num_now(i),  i = 1, num_step)
     write(iwrite) (st_halo(i),  i = 1, num_step)

     ! from mrgt
     write(iwrite) (f_des(i), i = 0, num_tot)
     write(iwrite) (n_des(i), i = 0, num_tot)
     write(iwrite) (f_prg(i), i = 0, num_tot)
     write(iwrite) (numb(i),  i = 0, num_tot)
     write(iwrite) (hori(i),  i = 0, num_tot)
     write(iwrite) (mpi(i),   i = 0, num_tot)
     write(iwrite) (mhalo(i), i = 0, num_tot)

     close(iwrite)
  ENDDO


  ! --- write redshift information into "fname_z"
  iwrite = 1; fname_z = trim(base)//'sspfiles/redshifts.dat'
  open(iwrite, file=trim(fname_z), status='unknown')
  write(iwrite, '(A)') '# List of Redshifts of sspfiles for '//trim(c_nbody)//' N-body run'
  DO iz = 1, num_step
     write(cz, '(I3.3)') iz - 1
     write(iwrite, '(A, 1X, F10.5)') trim(cz), zplus1ar(iz) - 1.d0
  ENDDO
  close(iwrite)
END PROGRAM ascii2bin
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
