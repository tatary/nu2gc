  DO i = 1, mcmc(iforest)%mrgp%num_now(end_step)
     me   = mcmc(iforest)%mrgp%st_halo(end_step) + i - 1
     totd = totd + mcmc(iforest)%mrgt%mhalo(me)
  ENDDO

  iout = ionum + 1
  DO_igal: DO igal = 1, endhalo
     Mssat = mcmc(iforest)%gal(igal)%Mstarb + mcmc(iforest)%gal(igal)%Mstard
     tots  = tots + Mssat; totc = totc + mcmc(iforest)%gal(igal)%Mcool
     toth  = toth + mcmc(iforest)%gal(igal)%Mhot + mcmc(iforest)%gal(igal)%Mhotorg

     totlum = 0.d0
     DO j = 1, param%tnw
        totlum = totlum + mcmc(iforest)%gal(igal)%lumg(j)
     ENDDO

     IF_empty: IF(Mssat < 1.d-13 .or. totlum <= 0.) THEN
        IF(mcmc(iforest)%mrgt%num_g(mcmc(iforest)%gal(igal)%IDhost) > 1) &
             print '(A,I2, 5(A,I7), 4(A,G10.3))', &
             'flag_c=', mcmc(iforest)%gal(igal)%flag_c, ', igal=', igal, ', ID(c_gal)=', &
             mcmc(iforest)%gal(igal)%id_cgal, ', ID(host)=', mcmc(iforest)%gal(igal)%IDhost, ', ngal(host)=', &
             mcmc(iforest)%mrgt%num_g(mcmc(iforest)%gal(igal)%IDhost), ', ID(prog)=', mcmc(iforest)%gal(igal)%IDprog, &
             ', Mprog[Msun]=', mcmc(iforest)%mrgt%mhalo(mcmc(iforest)%gal(igal)%IDprog) * param%munit, &
             ', Mhot[Msun]=', mcmc(iforest)%gal(igal)%Mhot * param%munit, &
             ', Mhotorg[Msun]=', mcmc(iforest)%gal(igal)%Mhotorg * param%munit, &
             ', Mcool[Msun]=', mcmc(iforest)%gal(igal)%Mcool * param%munit
     ELSE
        mag(param%iBband) = const%corr &
             - 2.5d0 * log10(mcmc(iforest)%gal(igal)%lumg(param%iBband) &
                             + mcmc(iforest)%gal(igal)%lumg(param%iBband+param%nwave))

!!$            IF(mag(param%iBband) <= -5.) THEN
!!$            IF(mag(param%iBband) <= -12.) THEN
        IF_bright: IF(mag(param%iBband) <= 100.) THEN
           rb0 = CalBulgeReff(mcmc(iforest)%gal(igal)%rbulge) ! [kpc/h]
           rd0 = CalDiskReff(mcmc(iforest)%gal(igal)%rdisk, rb0) ! [kpc/h]

!!$           lumtmp(:) = gal(igal)%lumg(:) / param%ups
           lumtmp(:) = mcmc(iforest)%gal(igal)%lumg(:)
           mor = DetMorType(lumtmp(param%iBband), lumtmp(param%iBband+param%nwave))
                 ! morphology w/o dust (1:E, 2:S0, 3:S)
           mag(1:param%nwave) = const%corr &
                - 2.5d0 * log10(lumtmp(1:param%nwave)&
                                + lumtmp(param%nwp1:param%tnw))

           ! --- calculate dust-extincted luminosities 'lumtmp_d(:)'
           !       and magnitudes 'mag_d(:)'
           IF(mcmc(iforest)%gal(igal)%flag_burst == 1) THEN ! starburst
              IF(mcmc(iforest)%gal(igal)%MZc_rem > 0.d0 .or. mcmc(iforest)%gal(igal)%MZc > 0.d0) THEN
                 MZctmp = max(mcmc(iforest)%gal(igal)%MZc_rem, mcmc(iforest)%gal(igal)%MZc)
                 call opt2(lumtmp, lumtmp_d, MZctmp, rb0/param%h)
                 mag_d(1:param%nwave) = const%corr &
                      - 2.5d0 * log10(lumtmp_d(1:param%nwave)&
                                      + lumtmp_d(param%nwp1:param%tnw))
              ELSE
                 lumtmp_d(:) = lumtmp(:); mag_d(:) = mag(:)
              ENDIF
           ELSE ! quiescent
              IF(mcmc(iforest)%gal(igal)%MZc > 0.d0) THEN
                 MZctmp = mcmc(iforest)%gal(igal)%MZc
                 call opt(lumtmp, lumtmp_d, MZctmp, rd0/param%h)
                 mag_d(1:param%nwave) = const%corr &
                      - 2.5d0 * log10(lumtmp_d(1:param%nwave) &
                                      + lumtmp_d(param%nwp1:param%tnw))
              ELSE
                 lumtmp_d(:) = lumtmp(:); mag_d(:) = mag(:)
              ENDIF
           ENDIF

!!$               bt_d = bt
!!$           bt_d = lumtmp_d(param%iBband) &
!!$                / (lumtmp_d(param%iBband) + lumtmp_d(param%iBband+param%nwave))
!!$           mor_d = CalMorphology(bt_d) ! morphology w/ dust (1:E,2:S0,3:S)
           mor_d = DetMorType(lumtmp_d(param%iBband), &
                              lumtmp_d(param%iBband+param%nwave))
                   ! morphology w/ dust (1:E,2:S0,3:S)
           IF(mor_d == 1 .or. mor_d == 2) THEN ! E or S0
              b_or_d = 1 ! bulge-dominated
              Vdomi  = mcmc(iforest)%gal(igal)%Vbulge; Vmino = mcmc(iforest)%gal(igal)%Vdisk
              rdomi  = rb0; rmino = rd0
           ELSE ! S
              b_or_d = 2 ! disk-dominated
              Vdomi  = mcmc(iforest)%gal(igal)%Vdisk; Vmino = mcmc(iforest)%gal(igal)%Vbulge
              rdomi  = rd0; rmino = rb0
           ENDIF


!!$     ==== lf(nz)%n(i,j,k): j=1,2,3 for E,S0,S LFs
           DO j = 1, param%nwave
              bin = int(anint((mag(j) - lf(nz)%base) * lf(nz)%invstep))
              IF(bin <= lf(nz)%Nbin .and. bin >= 1) &
                   lf(nz)%n(j, mor_d, bin) = lf(nz)%n(j, mor_d, bin) + inv_V
           ENDDO


!!$           flag_morph = mor_d
!!$           ih = gal(igal)%IDhost
!!$     ==== lf_d(nz)%n(i,j,k): j=1,2,3 for E,S0,S LFs
           DO j = 1, param%nwave
              bin = int(anint((mag_d(j) - lf_d(nz)%base) * lf_d(nz)%invstep))
              IF(bin <= lf_d(nz)%Nbin .and. bin >= 1) &
                   lf_d(nz)%n(j, mor_d, bin) = lf_d(nz)%n(j, mor_d, bin) + inv_V
           ENDDO

           temp = 1.d0 / (lumtmp(param%iVband) &
                          + lumtmp(param%iVband+param%nwave))
           zt = (mcmc(iforest)%gal(igal)%Mtb + mcmc(iforest)%gal(igal)%Mtd) * temp * param%th
!!$           zt = (gal(igal)%Mtb + gal(igal)%Mtd) / (gal(igal)%Mstarb + gal(igal)%Mstard)
           zg = (mcmc(iforest)%gal(igal)%MZb + mcmc(iforest)%gal(igal)%MZd) * temp

!!$               IF(mor_d == 3) THEN !S
!!$                  IF(lumtmp(param%iBband) <= 0.d0) THEN
!!$                     gal(igal)%Vbulge = 0.d0; rb0 = 0.d0
!!$                  ENDIF
!!$               ENDIF


!!$             IF(mor_d == 1) THEN
!!$                write(31,'(7(E13.5),2I9)') mag_d(param%iVband),mag_d(param%iUband) - mag_d(param%iVband),&
!!$                    mag_d(param%iVband) - mag_d(param%iKband),zg,zt/param%ups,gal(igal)%Mstarb + gal(igal)%Mstard,&
!!$                    gal(igal)%Vc,gal(igal)%hori,gal(igal)%mpi
!!$                write(31,'(5E13.5)') mag_d(param%iBband),mag_d(param%iIcband),zg,&
!!$                    gal(igal)%Mstarb + gal(igal)%Mstard,gal(igal)%Mcool
!!$                IF(mrgt%mhalo(gal(igal)%IDhost) >= 10.)&
!!$                 write(41,'(6E13.5)') mag_d(param%iVband),mag_d(param%iUband) - mag_d(param%iVband),&
!!$                    mag_d(param%iVband) - mag_d(param%iKband),zg,zt,mrgt%mhalo(gal(igal)%IDhost)
!!$                ENDIF
!!$             ELSEIF(mor_d == 2) THEN
!!$                write(32,'(7(E13.5))') mag_d(param%iVband),mag_d(param%iUband) - mag_d(param%iVband),&
!!$                    mag_d(param%iVband) - mag_d(param%iKband),zg,zt/param%ups,gal(igal)%Mstarb + gal(igal)%Mstard,&
!!$                    gal(igal)%Vc
!!$             ELSE
!!$                write(33,'(7(E13.5))') mag_d(param%iVband),mag_d(param%iUband) - mag_d(param%iVband),&
!!$                    mag_d(param%iVband) - mag_d(param%iKband),zg,zt/param%ups,gal(igal)%Mstarb + gal(igal)%Mstard,&
!!$                    gal(igal)%Vc
!!$                  write(29,'(7E13.5)')(mag_d(k),k = 1, param%nwave),gal(igal)%Vdisk,gal(igal)%Vc
!!$             ENDIF


!!$    ===== QSO ======
           DO j = 1, param%nwaveq
              IF(mcmc(iforest)%gal(igal)%lumq(j) > 0.d0) THEN
                 mag_q(j) = const%corr - 2.5d0 * log10(mcmc(iforest)%gal(igal)%lumq(j))
              ELSE
                 mag_q(j) = 128.d0
              ENDIF
           ENDDO

!!$     ==== lf_q(nz)%n(i,j,k): j=1 only
           DO j = 1, param%nwaveq
              bin = int(anint((mag_q(j) - lf_q(nz)%base) * lf_q(nz)%invstep))
              IF(bin <= lf_q(nz)%Nbin .and. bin >= 1 .and. mag_d(param%iBband) > mag_q(1) ) &
                   lf_q(nz)%n(j, 1, bin) = lf_q(nz)%n(j, 1, bin) + inv_V
           ENDDO

!!$----------------------------------------------------------------------------
!           write(iout, '(2I10, 4I2, I8, 3I16, 18G13.5, 12F11.4, 3G13.5)') &
!                igal, gal(igal)%id_cgal, & ! 1,2
!                mor, mor_d, gal(igal)%flag_c, gal(igal)%flag_burst, & ! 3,4,5,6
!                iforest, & ! 7
!                gal(igal)%IDhost, gal(igal)%mpi, mrgt%mpi(gal(igal)%IDhost), & ! 8,9,10
!                gal(igal)%Mstarb * param%munit, & ! 11
!                gal(igal)%Mstard * param%munit, & ! 12
!                gal(igal)%Mcool  * param%munit, & ! 13
!                mrgt%mhalo(gal(igal)%IDhost) * param%munit, & ! 14
!                gal(igal)%Mhalo * param%munit,            & ! 15
!                gal(gal(igal)%id_cgal)%Vcent, gal(igal)%Vcent, & ! 16,17
!                gal(igal)%MZc * param%munit, zg, zt, & ! 18,19,20
!                rb0, gal(igal)%Vbulge, rd0, gal(igal)%Vdisk, & ! 21,22,23,24
!                lumtmp(param%iBband), lumtmp_d(param%iBband), & ! 25,26
!                lumtmp(param%iBband+param%nwave), lumtmp_d(param%iBband+param%nwave), & ! 27,28
!                mag(param%iBband),  mag_d(param%iBband),  & ! 29,30
!                mag(param%iVband),  mag_d(param%iVband),  & ! 31,32
!                mag(param%iKpband), mag_d(param%iKpband), & ! 33,34
!                mag(param%iBband_r),  mag_d(param%iBband_r),  & ! 35,36
!                mag(param%iVband_r),  mag_d(param%iVband_r),  & ! 37,38
!                mag(param%iKpband_r), mag_d(param%iKpband_r), & ! 39,40
!                mag_d(param%iBband) + 5.d0*log10(rdomi) + 38.568d0, & ! 41
!                4.62d+5 * rdomi * Vdomi**2 &
!                  * 10.d0**(0.4d0*(mag_d(param%iBband)-5.48d0)) / param%h, & ! 42
!                4.62d+5 * rdomi * Vdomi**2 / param%h ! 43
!------------------------
!           IF((mag_d(param%iBband) < - 20.d0).and.(gal(igal)%lumq(1) > 0.d0)) THEN
!           IF(mag_q(1) <- 15.d0) THEN
!!$           IF(mag_d(param%iBband) < - 18.d0 .or. mag_q(1) < -20.0) THEN
!!$               write(iout, '(3I16, 8E13.5, 1I2)') &
!!$                gal(igal)%mpi, &
!!$                gal(igal)%IDhost, &
!!$                mrgt%mpi(gal(igal)%IDhost), &
!!$!                mrgt%mhalo(gal(igal)%IDhost) * param%munit, &
!!$                gal(igal)%Mcool * param%munit, &
!!$                gal(igal)%Mstard * param%munit, &
!!$                gal(igal)%Mstarb * param%munit, &
!!$                gal(igal)%Mbh * param%munit, &
!!$                mag_d(param%iBband), &
!!$!                   mag_d(param%iKpband), &
!!$                mag_q(1), &
!!$                gal(igal)%lumq(2), gal(igal)%lumq(3), &
!!$                gal(igal)%flag_c
!!$            ENDIF

!!$           IF (gal(igal)%Mstarb * gal(igal)%Mbh > 0.d0) THEN
!!$              write(iout, '(2E13.5, 1I4)') &
!!$                   log10(gal(igal)%Mstarb)+param%log10munit,&
!!$                   log10(gal(igal)%Mbh)+param%log10munit, mor_d
!!$           ENDIF

!!$ ==== mass functions ====
           ! mf(:)%n(i,j,k)
           ! --- i:MFtype(1:Mstar^bulge, 2:Mstar^disk, 3:Mcold, 4:MBH,
           !              5:Mhot, 6:Mstar^bulge+MBH, 7:Mstar^disk+Mcold,
           !              8:Mstar^bulge+Mstar^disk+MBH+Mcold,
           !              9:Mstar^bulge+Mstar^disk, 10:hydrogen gas, 11:Mhalo),
           !     j:mor(1:E,2:S0,3:S,4:all), k:mass
            ! 1: bulge star mass
           IF(mcmc(iforest)%gal(igal)%Mstarb > EPS) THEN
              bin = int(anint((log10(mcmc(iforest)%gal(igal)%Mstarb) + mf(nz)%base) * mf(nz)%invstep))
              IF(bin <= mf(nz)%Nbin .and. bin >= 1) THEN
                 mf(nz)%n(1, mor_d, bin) = mf(nz)%n(1, mor_d, bin) + inv_V
                 mf(nz)%n(1, 4, bin)     = mf(nz)%n(1, 4, bin)     + inv_V
              ENDIF
           ENDIF
           ! 2: disk star mass
           IF(mcmc(iforest)%gal(igal)%Mstard > EPS) THEN
              bin = int(anint((log10(mcmc(iforest)%gal(igal)%Mstard) + mf(nz)%base) * mf(nz)%invstep))
              IF(bin <= mf(nz)%Nbin .and. bin >= 1) THEN
                 mf(nz)%n(2, mor_d, bin) = mf(nz)%n(2, mor_d, bin) + inv_V
                 mf(nz)%n(2, 4, bin)     = mf(nz)%n(2, 4, bin)     + inv_V
              ENDIF
           ENDIF
           ! 3: cold gas mass
           IF(mcmc(iforest)%gal(igal)%Mcool > EPS) THEN
              bin = int(anint((log10(mcmc(iforest)%gal(igal)%Mcool) + mf(nz)%base) * mf(nz)%invstep))
              IF(bin <= mf(nz)%Nbin .and. bin >= 1) THEN
                 mf(nz)%n(3, mor_d, bin) = mf(nz)%n(3, mor_d, bin) + inv_V
                 mf(nz)%n(3, 4, bin)     = mf(nz)%n(3, 4, bin)     + inv_V
              ENDIF
           ENDIF
           ! 4: SMBH mass
           IF(mcmc(iforest)%gal(igal)%Mbh > EPS) THEN
              bin = int(anint((log10(mcmc(iforest)%gal(igal)%Mbh) + mf(nz)%base) * mf(nz)%invstep))
              IF(bin <= mf(nz)%Nbin .and. bin >= 1) THEN
                 mf(nz)%n(4, mor_d, bin) = mf(nz)%n(4, mor_d, bin) + inv_V
                 mf(nz)%n(4, 4, bin)     = mf(nz)%n(4, 4, bin)     + inv_V
              ENDIF
           ENDIF
           ! 5: hot gas mass
           IF(mcmc(iforest)%gal(igal)%Mhot > EPS) THEN
              bin = int(anint((log10(mcmc(iforest)%gal(igal)%Mhot) + mf(nz)%base) * mf(nz)%invstep))
              IF(bin <= mf(nz)%Nbin .and. bin >= 1) THEN
                 mf(nz)%n(5, mor_d, bin) = mf(nz)%n(5, mor_d, bin) + inv_V
                 mf(nz)%n(5, 4, bin)     = mf(nz)%n(5, 4, bin)     + inv_V
              ENDIF
           ENDIF
           ! 6: bulge mass
           Mbulge = mcmc(iforest)%gal(igal)%Mstarb + mcmc(iforest)%gal(igal)%Mbh
           IF(Mbulge > EPS) THEN
              bin = int(anint((log10(Mbulge) + mf(nz)%base) * mf(nz)%invstep))
              IF(bin <= mf(nz)%Nbin .and. bin >= 1) THEN
                 mf(nz)%n(6, mor_d, bin) = mf(nz)%n(6, mor_d, bin) + inv_V
                 mf(nz)%n(6, 4, bin)     = mf(nz)%n(6, 4, bin)     + inv_V
              ENDIF
           ENDIF
           ! 7: disk mass
           Mdisk = mcmc(iforest)%gal(igal)%Mstard + mcmc(iforest)%gal(igal)%Mcool
           IF(Mdisk > EPS) THEN
              bin = int(anint((log10(Mdisk) + mf(nz)%base) * mf(nz)%invstep))
              IF(bin <= mf(nz)%Nbin .and. bin >= 1) THEN
                 mf(nz)%n(7, mor_d, bin) = mf(nz)%n(7, mor_d, bin) + inv_V
                 mf(nz)%n(7, 4, bin)     = mf(nz)%n(7, 4, bin)     + inv_V
              ENDIF
           ENDIF
           ! 8: galaxy mass
           galmass = Mbulge + Mdisk
           IF(galmass > EPS) THEN
              bin = int(anint((log10(galmass) + mf(nz)%base) * mf(nz)%invstep))
              IF(bin <= mf(nz)%Nbin .and. bin >= 1) THEN
                 mf(nz)%n(8, mor_d, bin) = mf(nz)%n(8, mor_d, bin) + inv_V
                 mf(nz)%n(8, 4, bin)     = mf(nz)%n(8, 4, bin)     + inv_V
              ENDIF
           ENDIF
           ! 9: stellar mass
           IF(Mssat > EPS) THEN
              bin = int(anint((log10(Mssat) + mf(nz)%base) * mf(nz)%invstep))
              IF(bin <= mf(nz)%Nbin .and. bin >= 1) THEN
                 mf(nz)%n(9, mor_d, bin) = mf(nz)%n(9, mor_d, bin) + inv_V
                 mf(nz)%n(9, 4, bin)     = mf(nz)%n(9, 4, bin)     + inv_V
              ENDIF
           ENDIF
           ! 10: Hydrogen gas mass
           M_H = mcmc(iforest)%gal(igal)%Mcool * const%fracH ! Hydrogen gas mass in cold gas
           IF(M_H > EPS) THEN
              bin = int(anint((log10(M_H) + mf(nz)%base) * mf(nz)%invstep))
              IF(bin <= mf(nz)%Nbin .and. bin >= 1) THEN
                 mf(nz)%n(10, mor_d, bin) = mf(nz)%n(10, mor_d, bin) + inv_V
                 mf(nz)%n(10, 4, bin)     = mf(nz)%n(10, 4, bin)     + inv_V
              ENDIF
           ENDIF
           ! 11: halo mass
           IF(mcmc(iforest)%gal(igal)%flag_c == 1) THEN ! count only central galaxy
              ! "gal(igal)%Mhalo" is identical to "mrgt%mhalo(gal(i)%IDhost)" for central
!!$              temp = gal(igal)%Mhalo / mrgt%mhalo(gal(igal)%IDhost) - 1.d0
!!$              IF(temp > EPS .or. temp < -EPS) print '(I10, E12.5)', igal, temp
              bin = int(anint((log10(mcmc(iforest)%gal(igal)%Mhalo) + mf(nz)%base) * mf(nz)%invstep))
              IF(bin <= mf(nz)%Nbin .and. bin >= 1) THEN
                 mf(nz)%n(11, mor_d, bin) = mf(nz)%n(11, mor_d, bin) + inv_V
                 mf(nz)%n(11, 4, bin)     = mf(nz)%n(11, 4, bin)     + inv_V
              ENDIF
           ENDIF

!!$ ==== M_HI / L_B as a function of M_B-5logh ===
           MassLumi = const%ML * mcmc(iforest)%gal(igal)%Mcool &
                / (mcmc(iforest)%gal(igal)%lumg(param%iBband) + mcmc(iforest)%gal(igal)%lumg(param%iBband+param%nwave))
           bin = int(anint((mag_d(param%iBband) - ml(nz)%base) * ml(nz)%invstep))
           IF(bin <= ml(nz)%Nbin .and. bin >= 1) THEN
              ! for mean and sigma
              ml(nz)%n(1, mor_d, bin) = ml(nz)%n(1, mor_d, bin) + inv_V
              ml(nz)%n(1, 4,     bin) = ml(nz)%n(1, 4,     bin) + inv_V
              IF(MassLumi > 0.d0) THEN
                 temp = inv_V * MassLumi
                 ml(nz)%xn(1, mor_d, bin) = ml(nz)%xn(1, mor_d, bin) + temp
                 ml(nz)%xn(1, 4,     bin) = ml(nz)%xn(1, 4,     bin) + temp
                 temp = temp * MassLumi
                 ml(nz)%xxn(1, mor_d, bin) = ml(nz)%xxn(1, mor_d, bin) + temp
                 ml(nz)%xxn(1, 4,     bin) = ml(nz)%xxn(1, 4,     bin) + temp
              ENDIF

              ! for median and quartiles
              ml(nz)%n(2, mor_d, bin) = ml(nz)%n(2, mor_d, bin) + 1.d0
              ml(nz)%n(2, 4,     bin) = ml(nz)%n(2, 4,     bin) + 1.d0
              j = int(ml(nz)%n(2, mor_d, bin)); k = int(ml(nz)%n(2, 4, bin))
              ml(nz)%x(j, mor_d, bin) = MassLumi; ml(nz)%x(k, 4, bin) = MassLumi
           ENDIF

           ! --- for LAE related
           IF(param%LAE) &
                call MainCalcForLAE(igal, nz, iforest, mcmc(iforest)%mrgt%mhalo(mcmc(iforest)%gal(igal)%IDhost), &
                                    inv_V, rb0, rd0, mor, mor_d, Mssat)
        ENDIF IF_bright ! mag(param%iBband) <= 100.


!!$        ! --- for SFH related
!!$        IF(param%SFH) &
!!$             call SubstReffTauV(igal, mor_d, rd0, rb0, param%tauV0*gal(igal)%MZc)

!!$        ! --- for SN related
!!$        IF(param%SN) THEN
!!$           IF(Mssat*param%munit >= paramSN%thMstar) THEN
!!$              LAE%mSFR = LAE%SFR
!!$              IF(LAE%sftype == 2) THEN ! starburst
!!$                 LAE%mSFR = CalMeanSFR(LAE%tel, LAE%twind, LAE%taust_yr, &
!!$                                       LAE%Mgas_pre, LAE%Ms_b, LAE%beta, &
!!$                                       LAE%ab)
!!$                 IF(LAE%sftype2 == 3 .and. Mcool(LAE%id) > 0.d0) &
!!$                      ! multiple merger
!!$                      LAE%mSFR = LAE%mSFR + allgal(LAE%id)%sfr(1)
!!$              ENDIF
!!$              call CalAndWriteSNrate(igal, mor_d, inv_V, Mssat*param%munit, &
!!$                                     reff(igal), tauV_SFH(igal), LAE%Tmass, &
!!$                                     LAE%mSFR, LAE%magd)
!!$           ENDIF
!!$        ENDIF
     ENDIF IF_empty ! Mssat >= 1.d-13 .and. totlum > 0.d0
  ENDDO DO_igal ! i=1,endhalo
