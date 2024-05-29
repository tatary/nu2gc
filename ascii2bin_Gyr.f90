PROGRAM main
  implicit none
  INTEGER, PARAMETER:: N=25000000,NN=200,NNN=200000
!!$    === mgrt ===
  INTEGER :: num_tot
  INTEGER :: f_des(0:N), n_des(0:N), f_prg(0:N), numb(0:N)
  INTEGER :: hori(0:N), mpi(0:N)
  DOUBLE PRECISION :: mhalo(0:N)
!!$    === mgrp ===
  INTEGER :: num_now(NN), st_halo(NN), num_step
  DOUBLE PRECISION :: times(NN), zplus1ar(NN), deltac, deltastop, zsp1
  DOUBLE PRECISION :: lumg(10,NNN)
  INTEGER   :: idum, i, loop, j, ii
  CHARACTER :: ci*3, file_mrgt*50, file_mrgp*50, file_bin*50
  CHARACTER :: nbodyrun*20
  DOUBLE PRECISION :: x

  nbodyrun  = 'w790100a'
  file_mrgt = trim(nbodyrun)//'/'//trim(nbodyrun)//'.000.mrgt'
  file_mrgp = file_mrgt(1:len_trim(file_mrgt)-1)//'p'
  ii = index(file_mrgt, '000')
  DO loop = 0, 7
!!$     do loop=0,0
     write(ci, '(I3.3)') loop

     f_des = 0; n_des = 0; f_prg = 0; numb = 0; hori = 0; mpi = 0
     mhalo = 0.d0; times = 0.d0; num_now = 0; st_halo = 0; zplus1ar = 0.d0

     file_mrgt(ii:ii+2) = ci; file_mrgp(ii:ii+2) = ci
     print '(A)', trim(file_mrgt)//', '//trim(file_mrgp)
!!$    ========================================================
     j = 1
     open(12, file = file_mrgp, status = 'old')
30   read(12, *, end = 40) times(j), x, num_now(j), st_halo(j)
     ! times(:) : [Gyr]
     if(x > 1.d0) goto 40
     zplus1ar(j) = 1.d0 / x
     j = j + 1
     goto 30
40   continue
     close(12)
     num_step = j - 1

     open(11, file = file_mrgt, status = 'old')
10   read(11, *, end=20) i, f_des(i), n_des(i), f_prg(i), numb(i),&
          hori(i), mpi(i), mhalo(i)
     goto 10
20   continue
     close(11)
     num_tot = i

     file_bin = file_mrgt(1:len_trim(file_mrgt)-1)//'b'
     write(*, *) file_bin
     open(1, file = file_bin, form = 'unformatted')
     write(1) num_step, num_tot
     print '(2I8)', num_step, num_tot

     ! from mrgp
     write(1) (times(i)*10.0**3,    i = 1, num_step)  ! Gyr -> Myr
     write(1) (zplus1ar(i), i = 1, num_step)
     write(1) (num_now(i),  i = 1, num_step)
     write(1) (st_halo(i),  i = 1, num_step)

     ! from mrgt
     write(1) (f_des(i), i = 0, num_tot)
     write(1) (n_des(i), i = 0, num_tot)
     write(1) (f_prg(i), i = 0, num_tot)
     write(1) (numb(i),  i = 0, num_tot)
     write(1) (hori(i),  i = 0, num_tot)
     write(1) (mpi(i),   i = 0, num_tot)
     write(1) (mhalo(i), i = 0, num_tot)

     close(1)
  ENDDO
end program main
