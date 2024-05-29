  print '(4(A, F10.3))', '# tots:', tots, ', totc:', totc, &
       ', toth:', toth, ', tots+totc+toth:', tots+totc+toth
  print '(2(A, F10.3))', '  totd:', totd, &
       ', (tots+totc+toth)/totd:', (tots+totc+toth)/totd

  ! Okamoto+08 UV feedback mode
  IF(param%UVfb == 2 .and. param%iO08_p+param%iO08_n0 > 0) THEN
     temp = 100.d0 / dble(param%iO08_p + param%iO08_n0)
     print '(A, G10.3)', '   % of halo w/ Macc >= 0: ', &
          dble(param%iO08_p)  * temp
     print '(A, G10.3)', '   % of halo w/ Macc  < 0: ', &
          dble(param%iO08_n0) * temp
     print '(A, G10.3)', '        --- Macc+tmp%Mhotrog '//&
          '< 0                   : ', dble(param%iO08_n1) * temp
     print '(A, G10.3)', '        --- Macc+tmp%Mhotorg+'//&
          'tmp%Mhot < 0          : ', dble(param%iO08_n2) * temp
     print '(A, G10.3)', '        --- Macc+tmp%Mhotorg+'//&
          'tmp%Mhot+Mcool_all < 0: ', dble(param%iO08_n3) * temp
     print '(A, G10.3)', '        --- Macc+Mbar_all < 0'//&
          '                      : ', dble(param%iO08_n4) * temp
  ENDIF

  IF(param%run_type /= 3) THEN
     ! --- luminosity functions written into 'ionum+2' and 'ionum+3'
     ! lf(:)%n(i,j,k) and lf_d(:)%n(i,j,k) --- i:band, j:mor(1:E,2:S0,3:S,4:all), k:mag
     DO i = 1, lf(nz)%Nbin
        ! -------------
        !    total LF
        ! -------------
        IF(param%run_type == 2) THEN !all redshift 
           iout = ionum + NFILE * (nz-1) + 2
        ELSE ! other
           iout = ionum + 2
        ENDIF

        write(iout, '(2X, F9.4, $)') lf(nz)%bin(i)
        DO j = 1, param%nwave
           lf(nz)%n(j, 4, i) = lf(nz)%n(j,1,i) + lf(nz)%n(j,2,i) + lf(nz)%n(j,3,i)
           write(iout, '(2X, E13.6, $)') lf(nz)%n(j, 4, i) * lf(nz)%invstep
        ENDDO
        DO j = 1, param%nwave
           lf_d(nz)%n(j, 4, i) = lf_d(nz)%n(j,1,i) + lf_d(nz)%n(j,2,i) &
                                 + lf_d(nz)%n(j,3,i)
           write(iout, '(2X, E13.6, $)') lf_d(nz)%n(j, 4, i) * lf(nz)%invstep
        ENDDO
        ! --- contribution from burst
        DO j = 1, param%nwave
           lf(nz)%n_brst(j, 4, i) = lf(nz)%n_brst(j,1,i) + lf(nz)%n_brst(j,2,i) &
                                    + lf(nz)%n_brst(j,3,i)
           IF(j == 1) THEN
              write(iout, '(12X, E13.6, $)') lf(nz)%n_brst(j, 4, i) * lf(nz)%invstep
           ELSE
              write(iout, '(2X, E13.6, $)') lf(nz)%n_brst(j, 4, i) * lf(nz)%invstep
           ENDIF
        ENDDO
        DO j = 1, param%nwave
           lf_d(nz)%n_brst(j, 4, i) = lf_d(nz)%n_brst(j,1,i) &
                                      + lf_d(nz)%n_brst(j,2,i) &
                                      + lf_d(nz)%n_brst(j,3,i)
           write(iout, '(2X, E13.6, $)') lf_d(nz)%n_brst(j, 4, i) * lf(nz)%invstep
        ENDDO
        write(iout, *)


        ! ----------------------------
        !    morphology-dependent LF
        ! ----------------------------
        IF(param%run_type == 2) THEN !all redshift 
           iout = ionum + NFILE * (nz-1) + 3
        ELSE ! other
           iout = ionum + 3
        ENDIF

        write(iout,'(2X, F9.4, $)') lf(nz)%bin(i)
        DO j = 1, param%nwave ! E  w/o dust
           write(iout, '(2X, E13.6, $)') lf(nz)%n(j, 1, i) * lf(nz)%invstep
        ENDDO
        DO j = 1, param%nwave ! S0 w/o dust
           write(iout, '(2X, E13.6, $)') lf(nz)%n(j, 2, i) * lf(nz)%invstep
        ENDDO
        DO j = 1, param%nwave ! S  w/o dust
           write(iout, '(2X, E13.6, $)') lf(nz)%n(j, 3, i) * lf(nz)%invstep
        ENDDO
        DO j = 1, param%nwave ! E  w/ dust
           write(iout, '(2X, E13.6, $)') lf_d(nz)%n(j, 1, i) * lf(nz)%invstep
        ENDDO
        DO j = 1, param%nwave ! S0 w/ dust
           write(iout, '(2X, E13.6, $)') lf_d(nz)%n(j, 2, i) * lf(nz)%invstep
        ENDDO
        DO j = 1, param%nwave ! S  w/ dust
           write(iout, '(2X, E13.6, $)') lf_d(nz)%n(j, 3, i) * lf(nz)%invstep
        ENDDO
        ! --- contribution from burst
        DO j = 1, param%nwave ! E  w/o dust
           IF(j == 1) THEN
              write(iout, '(12X, E13.6, $)') lf(nz)%n_brst(j, 1, i) * lf(nz)%invstep
           ELSE
              write(iout, '(2X, E13.6, $)')  lf(nz)%n_brst(j, 1, i) * lf(nz)%invstep
           ENDIF
        ENDDO
        DO j = 1, param%nwave ! S0 w/o dust
           write(iout, '(2X, E13.6, $)') lf(nz)%n_brst(j, 2, i) * lf(nz)%invstep
        ENDDO
        DO j = 1, param%nwave ! S  w/o dust
           write(iout, '(2X, E13.6, $)') lf(nz)%n_brst(j, 3, i) * lf(nz)%invstep
        ENDDO
        DO j = 1, param%nwave ! E  w/ dust
           write(iout, '(2X, E13.6, $)') lf_d(nz)%n_brst(j, 1, i) * lf(nz)%invstep
        ENDDO
        DO j = 1, param%nwave ! S0 w/ dust
           write(iout, '(2X, E13.6, $)') lf_d(nz)%n_brst(j, 2, i) * lf(nz)%invstep
        ENDDO
        DO j = 1, param%nwave ! S  w/ dust
           write(iout, '(2X, E13.6, $)') lf_d(nz)%n_brst(j, 3, i) * lf(nz)%invstep
        ENDDO
        write(iout, *)
     ENDDO


     ! --- mass functions written into 'ionum+4'
     IF(param%run_type == 2) THEN !all redshift 
        iout = ionum + NFILE * (nz-1) + 4
     ELSE ! other
        iout = ionum + 4
     ENDIF

     ! --- contribution from all galaxies
     DO i = 1, mf(nz)%Nbin
        write(iout,'(2X, F9.4, $)') mf(nz)%bin(i)
        ! mf(:)%n(i,j,k)
        ! --- i:MFtype(1:Mstar^bulge, 2:Mstar^disk, 3:Mcold, 4:MBH,
        !              5:Mhot, 6:Mstar^bulge+MBH, 7:Mstar^disk+Mcold,
        !              8:Mstar^bulge+Mstar^disk+MBH+Mcold,
        !              9:Mstar^bulge+Mstar^disk, 10:hydrogen gas, 11:Mhalo),
        !     j:mor(1:E,2:S0,3:S,4:all),
        !     k:mass
        DO j = 1, mf(nz)%N1
           mf(nz)%n(j, 4, i) = mf(nz)%n(j,1,i) + mf(nz)%n(j,2,i) + mf(nz)%n(j,3,i)
           write(iout, '(2X, E13.6, $)') mf(nz)%n(j, 4, i) * mf(nz)%invstep
        ENDDO
        DO j = 1, mf(nz)%N1 ! contribution from burst
           mf(nz)%n_brst(j, 4, i) = mf(nz)%n_brst(j,1,i) + mf(nz)%n_brst(j,2,i) &
                                    + mf(nz)%n_brst(j,3,i)
           IF(j == 1) THEN
              write(iout, '(12X, E13.6, $)') mf(nz)%n_brst(j, 4, i) * mf(nz)%invstep
           ELSE
              write(iout, '(2X, E13.6, $)')  mf(nz)%n_brst(j, 4, i) * mf(nz)%invstep
           ENDIF
        ENDDO
        write(iout, *)
     ENDDO
     write(iout, *); write(iout, *)
     ! --- contribution from the galaxies in each SF mode
     DO i = 1, mf(nz)%N1
        write(iout, '(A, I2, A)') '# index ', i, ': '//trim(c_mftype(i))
        write(iout, '(A)') '# (1)log10[Mass/(Msun/h^2)] '//&
             '(2-4)dn/dlogM[h^3/Mpc^3/dex] (2:E, 3:S0, 4:S) '//&
             '(5-7)dn/dlogM[h^3/Mpc^3/dex] for burst (5:E, 6:S0, 7:S)'
        DO j = 1, mf(nz)%Nbin
           write(iout, '(2X, F9.4, $)') mf(nz)%bin(j)
           DO k = 1, 3
              write(iout, '(2X, E13.6, $)') mf(nz)%n(i, k, j) * mf(nz)%invstep
           ENDDO
           DO k = 1, 3 ! contribution from burst
              write(iout, '(2X, E13.6, $)') mf(nz)%n_brst(i, k, j) * mf(nz)%invstep
           ENDDO
           write(iout, *)
        ENDDO
        write(iout, *); write(iout, *)
     ENDDO

     ! --- M_H / L_B distribution written into 'ionum+5'
     IF(param%ML) THEN
        IF(param%run_type == 2) THEN !all redshift 
           iout = ionum + NFILE * (nz-1) + 5
        ELSE ! other
           iout = ionum + 5
        ENDIF

        DO i = 1, ml(nz)%Nbin
           write(iout, '(2X, F13.6, $)') ml(nz)%bin(i)
           DO j = 1, ml(nz)%N2
              ! --- mean and sigma
              IF(ml(nz)%n(1,j,i) > 0.d0) THEN
                 ml(nz)%xn(1,j,i)  = ml(nz)%xn(1,j,i)  / ml(nz)%n(1,j,i)
                 ml(nz)%xxn(1,j,i) = ml(nz)%xxn(1,j,i) / ml(nz)%n(1,j,i)
                 IF(ml(nz)%xxn(1,j,i) - Square(ml(nz)%xn(1,j,i)) < 0.d0) THEN
                    ml(nz)%xxn(1, j, i) = 0.d0
                 ELSE
                    ml(nz)%xxn(1, j, i) = sqrt(ml(nz)%xxn(1,j,i) &
                                               - Square(ml(nz)%xn(1,j,i)))
                 ENDIF

                 IF(ml(nz)%xn(1,j,i) > ml(nz)%xxn(1,j,i)) THEN
                    temp = ml(nz)%xn(1,j,i) - ml(nz)%xxn(1,j,i)
                 ELSE
                    temp = 1.d-20
                 ENDIF
                 write(iout, '(3(2X, E13.6), $)') &
                      ml(nz)%xn(1,j,i), temp, ml(nz)%xn(1,j,i) + ml(nz)%xxn(1,j,i)
              ELSE
                 write(iout, '(3(2X, E13.6), $)') 1.d-20, 1.d-20, 1.d-20
              ENDIF


              ! --- median and quartiles
              bin = int(ml(nz)%n(2,j,i))
              IF(bin > 0) THEN
                 IF(bin == 1) THEN
                    IF(ml(nz)%x(1,j,i) < EPS) THEN
                       write(iout, '(3(2X, E13.6), $)') 1.d-20, 1.d-20, 1.d-20
                    ELSE
                       temp = sqrt(ml(nz)%x(1,j,i))
                       write(iout, '(3(2X, E13.6), $)') &
                            ml(nz)%x(1,j,i), ml(nz)%x(1,j,i)-temp, &
                            ml(nz)%x(1,j,i)+temp
                    ENDIF
                 ELSE
                    allocate(arr(bin))
                    DO k = 1, bin
                       arr(k) = ml(nz)%x(k,j,i)
                    ENDDO
                    call sort(bin, arr) ! sort arr(:) into ascending numerical order
                    IF(bin == 2) THEN
                       write(iout, '(3(2X, E13.6), $)') &
                            0.5d0*(arr(1)+arr(2)), arr(1), arr(2)
                    ELSEIF(bin == 3) THEN
                       write(iout, '(3(2X, E13.6), $)') arr(2), arr(1), arr(3)
                    ELSE
                       i25 = int(0.25d0*bin)
                       i50 = int(0.5d0 *bin)
                       i75 = int(0.75d0*bin)
                       IF(arr(i25) < EPS) THEN
                          temp = 1.d-20
                       ELSE
                          temp = arr(i25)
                       ENDIF
                       write(iout, '(3(2X, E13.6), $)') arr(i50), temp, arr(i75)
                    ENDIF
                    deallocate(arr)
                 ENDIF
              ELSE
                 write(iout, '(3(2X, E13.6), $)') 1.d-20, 1.d-20, 1.d-20
              ENDIF
           ENDDO
           write(iout, *)
        ENDDO
     ENDIF


     ! --- QSO LFs written into 'ionum+6'
     IF(param%run_type == 2) THEN !all redshift 
        iout = ionum + NFILE * (nz-1) + 6
     ELSE ! other
        iout = ionum + 6
     ENDIF

     DO i = 1, lf_q(nz)%Nbin
        write(iout, '(2X, F9.4, $)') lf_q(nz)%bin(i)
        DO j = 1, param%nwaveq
           write(iout, '(2X, E13.6, $)') lf_q(nz)%n(j,1,i) * lf(nz)%invstep
        ENDDO
        write(iout, *)
     ENDDO


     ! --- Write LAE related data into output files
     IF(param%LAE) THEN
        call WriteLyaLF; call WriteS03LFs
        call WriteDistributionsForLAE; call Write2DxLyaDistribution
        call CloseFilesForLAE(param%run_type)
     ENDIF


     ! --- Write SFH related data into output files
     IF(param%SFH .and. param%run_type == 1) THEN
        call WriteCSFH; call CloseFilesForSFH
     ENDIF


     ! --- Write SN related data into output files
     IF(param%SN .and. param%run_type == 1) THEN
        call WriteCosmicSNRateDensity
        call CloseFilesForSN
     ENDIF
  ENDIF
