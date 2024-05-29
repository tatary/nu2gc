  DO i = 1, mrgp%num_now(end_step)
     me   = mrgp%st_halo(end_step) + i - 1
     totd = totd + mrgt%mhalo(me)
  ENDDO

  ! galaxy catalog
  IF(param%run_type == 2) THEN ! all redshift
     iout = ionum + NFILE * nz_end + nnode * (nz-1) + inode + 1
  ELSE ! other
     iout = ionum + NFILE + inode + 1
  ENDIF

  DO_igal: DO igal = 1, endhalo ! "igal" is galaxy ID
     Mssat = gal(igal)%Mstarb + gal(igal)%Mstard
     tots  = tots + Mssat; totc = totc + gal(igal)%Mcoolb + gal(igal)%Mcoold
     toth  = toth + gal(igal)%Mhot + gal(igal)%Mhotorg

     ihost = gal(igal)%IDhost
     totlum = 0.d0
     DO iband = 1, param%tnw
        totlum = totlum + gal(igal)%lumg(iband)
     ENDDO
     IF_empty: IF((Mssat < 0.d0 .or. totlum <= 0.) &
                  .and. gal(igal)%z_col + 1.d0 > mrgp%zp1ar(mrgp%num_step-1)) THEN
        IF(mrgt%num_g(ihost) > 1 .and. gal(igal)%flag_c == 1) &
             print '(A,I2, 6(A,I7), 7(A,G10.3))', 'flag_c=', gal(igal)%flag_c, &
             ', iforest=', iforest, &
             ', igal=', igal, ', ID(c_gal)=', gal(igal)%id_cgal, &
             ', ID(host)=', ihost, ', ngal(host)=', mrgt%num_g(ihost), &
             ', ID(prog)=', gal(igal)%IDprog, &
             ', Mprog[Msun]=', mrgt%mhalo(gal(igal)%IDprog) * param%munit, &
             ', Mhot[Msun]=', gal(igal)%Mhot * param%munit, &
             ', Mhotorg[Msun]=', gal(igal)%Mhotorg * param%munit, &
             ', Mcool[Msun]=', (gal(igal)%Mcoold+gal(igal)%Mcoolb)*param%munit, &
             ', totlum =', totlum, &
             ', Mssat =', Mssat
     ELSE
        IF(gal(igal)%z_col + 1.d0 > mrgp%zp1ar(mrgp%num_step-1)) THEN
           mag(param%iBband_r) = 1000.d0
           temp = gal(igal)%lumg(param%iBband_r) &
                  + gal(igal)%lumg(param%iBband_r+param%nwave)
           IF(temp > 0.d0) &
                mag(param%iBband_r) = - 2.5d0 * log10(temp)

           IF_bright: IF(mag(param%iBband_r) <= 100.) THEN
              rb0 = CalBulgeReff(gal(igal)%rbulge) ! [kpc/h]
              rd0 = CalDiskReff(gal(igal)%rdisk, rb0) ! [kpc/h]

   !!$           lumtmp(:) = gal(igal)%lumg(:) / param%ups
              lumtmp(:) = gal(igal)%lumg(:)
              mor = DetMorType(lumtmp(param%iBband_r), &
                               lumtmp(param%iBband_r+param%nwave))
                    ! morphology w/o dust (1:E, 2:S0, 3:S)
              ! --- for debug (2014/Sep/22 by MARK)
              DO iband = 1, param%nwave
                 IF(lumtmp(iband)+lumtmp(iband+param%nwave) <= 0.d0) THEN
                    print '(A, I2, A, 2G11.3, I5, 2G11.3)', &
                         '# zero luminosity in iband = ', iband, &
                         '('//trim(ssp%bandname(param%iwave(iband)))//'): ', &
                         lumtmp(iband), lumtmp(iband+param%nwave), igal, &
                         gal(igal)%lumg(iband), gal(igal)%lumg(iband+param%nwave)
                 ENDIF
              ENDDO

              ! --- for debug (2014/Sep/22 by MARK)
              mag(1:param%nwave) = const%corr &
                   - 2.5d0 * log10(lumtmp(1:param%nwave) &
                                   + lumtmp(param%nwp1:param%tnw))
              ! --- added by MARK (2018/Mar/02)
              IF(param%wdim_band) THEN
                 DO iband = 1, param%nwave
                    IF((param%iwave(iband) >= 50 .and. param%iwave(iband) <= 58) &
                         .or. (param%iwave(iband) >= NWAVE_ALL+50 .and. &
                               param%iwave(iband) <= NWAVE_ALL+58)) THEN
                       ! NLyC, L{1216,1400,1500,1600,1700,2800,4861,6563}
                       i = param%iwave(iband) - 49
                       IF(param%iwave(iband) > NWAVE_ALL) i = i - NWAVE_ALL
                       mag(iband) = mag(iband) + const%corr_wdim(i)
                    ENDIF
                 ENDDO
              ENDIF

              ! --- calculate dust-extincted luminosities 'lumtmp_d(:)'
              !       and magnitudes 'mag_d(:)'
              IF((gal(igal)%flag_merger == 1 .or. gal(igal)%flag_di == 1) &
                 .and. gal(igal)%flag_burst == 1) THEN ! starburst
                 IF(gal(igal)%MZc_rem > 0.d0 .or. gal(igal)%MZcd > 0.d0) THEN
                    !MZctmp = max(gal(igal)%MZc_rem, gal(igal)%MZc)
                    !call opt2(lumtmp, lumtmp_d, MZctmp, rb0/param%h)
                    call opt2(lumtmp, lumtmp_d, gal(igal)%MZcd, &
                              gal(igal)%MZc_rem, rd0/param%h, rb0/param%h)
                    DO iband = 1, param%nwave
                       temp = lumtmp_d(iband) + lumtmp_d(iband+param%nwave)
                       !temp = lumtmp_d(iband)
                       mag_d(iband) = 99.d0
                       IF(temp > 0.d0) &
                            mag_d(iband) = const%corr - 2.5d0 * log10(temp)

                       ! --- added by MARK (2018/Mar/02)
                       IF(param%wdim_band) THEN
                          IF((param%iwave(iband) >= 50 .and. param%iwave(iband) <= 58) &
                               .or. (param%iwave(iband) >= NWAVE_ALL+50 .and. &
                                     param%iwave(iband) <= NWAVE_ALL+58)) THEN
                             ! NLyC, L{1216,1400,1500,1600,1700,2800,4861,6563}
                             i = param%iwave(iband) - 49
                             IF(param%iwave(iband) > NWAVE_ALL) i = i - NWAVE_ALL
                             mag_d(iband) = mag_d(iband) + const%corr_wdim(i)
                          ENDIF
                       ENDIF
                    ENDDO
                 ELSE
                    lumtmp_d(:) = lumtmp(:); mag_d(:) = mag(:)
                 ENDIF
                 ! for merging galaxies' mass ratio (shirakata; 2017/08/25)
                 IF(gal(igal)%agn(8) == gal(igal)%M2M1_av .and. &
                      gal(igal)%flag_merger == 1) THEN
                    gal(igal)%agn(8) = gal(igal)%agn(8) / gal(igal)%n_merge
                 ENDIF
              ELSE ! quiescent
                 IF(gal(igal)%MZcd > 0.d0) THEN
                    MZctmp = gal(igal)%MZcd
                    call opt(lumtmp, lumtmp_d, MZctmp, rd0/param%h)
                    DO iband = 1, param%nwave
                       temp = lumtmp_d(iband) + lumtmp_d(iband+param%nwave)
                       !temp = lumtmp_d(iband)
                       mag_d(iband) = 99.d0
                       IF(temp > 0.d0) &
                            mag_d(iband) = const%corr - 2.5d0 * log10(temp)

                       ! --- added by MARK (2018/Mar/02)
                       IF(param%wdim_band) THEN
                          IF((param%iwave(iband) >= 50 .and. param%iwave(iband) <= 58) &
                               .or. (param%iwave(iband) >= NWAVE_ALL+50 .and. &
                                     param%iwave(iband) <= NWAVE_ALL+58)) THEN
                             ! NLyC, L{1216,1400,1500,1600,1700,2800,4861,6563}
                             i = param%iwave(iband) - 49
                             IF(param%iwave(iband) > NWAVE_ALL) i = i - NWAVE_ALL
                             mag_d(iband) = mag_d(iband) + const%corr_wdim(i)
                          ENDIF
                       ENDIF
                    ENDDO
                 ELSE
                    lumtmp_d(:) = lumtmp(:); mag_d(:) = mag(:)
                 ENDIF
              ENDIF

              mor_d = DetMorType(lumtmp_d(param%iBband_r), &
                                 lumtmp_d(param%iBband_r+param%nwave))
                      ! morphology w/ dust (1:E,2:S0,3:S)
              ! --- Error message for the galaxies w/ mor_d = 0
              !      (modified by MARK, 2016/Jul/01)
              IF(param%run_type /= 3)THEN
                 IF(mor_d == 0) THEN ! morph. w/ dust ext. is not determined appropriately
                    print '(A)', '# There is a galaxy w/ mor_d = 0!!!'
                    print '(A)', '# This is because its dust-extincted '// &
                         'luminosities of disk and bulge are zero'
                    print '(A)', '# Its physical quantities are below:'

                    print '(A, 2I10, 7I2, I8, 4I16, 25G13.5, $)', '#', &
                         igal, gal(igal)%id_cgal, & ! 1,2
                         mor, mor_d, gal(igal)%flag_c, & ! 3--5
                         gal(igal)%flag_burst, gal(igal)%flag_merger, & ! 6,7
                         gal(igal)%flag_di, gal(igal)%flag_ccut, & ! 8,9
                         iforest, ihost, & ! 10,11
                         gal(igal)%mpi, mrgt%mpi(ihost), & ! 12,13
                         mrgt%f_des(ihost), gal(igal)%BT, & ! 14,15
                         gal(igal)%Mstarb  * param%munit, & ! 16
                         gal(igal)%Mstard  * param%munit, & ! 17
                         gal(igal)%Mcoolb  * param%munit, & ! 18
                         gal(igal)%Mcoold  * param%munit, & ! 19
                         mrgt%mhalo(ihost) * param%munit, & ! 20
                         gal(igal)%Mhalo   * param%munit, & ! 21
                         gal(igal)%z_col, & ! 22
                         gal(gal(igal)%id_cgal)%Vcent, gal(igal)%Vcent, & ! 23,24
                         gal(igal)%MZcb * param%munit, & ! 25
                         gal(igal)%MZcd * param%munit, & ! 26
                         zg, (const%toutGyr - zt), & ! 27,28 
                         rb0, gal(igal)%Vbulge, rd0, gal(igal)%Vdisk, & ! 29--32
                         mag_d(param%iBband_r) + 5.d0*log10(rdomi) + 38.568d0, & ! 33
                         4.62d+5 * rdomi * Square(Vdomi) &
                            * 10.d0**(0.4d0*(mag_d(param%iBband_r)-5.48d0)) &
                            / param%h, & ! 34
                         4.62d+5 * rdomi * Square(Vdomi) / param%h, & ! 35
                         gal(igal)%Mbh * param%munit, & ! 36
                         gal(igal)%MZc_rem * param%munit, & ! 37
                         0.0, 0.0    ! 38,39 (SFR, mSFR)
                    Zmassb = -99.d0; Zmassd = -99.d0
                    IF(gal(igal)%Mstarb > 0.d0) &
                         Zmassb = gal(igal)%Zmassb / gal(igal)%Mstarb
                    IF(gal(igal)%Mstard > 0.d0) &
                         Zmassd = gal(igal)%Zmassd / gal(igal)%Mstard
                    print '(3G13.5, $)', & ! 40--42
                         Zmassb, Zmassd, &
                         (gal(igal)%Zmassb + gal(igal)%Zmassd) / Mssat
                    IF(gal(igal)%flag_burst == 1 &
                         .and. gal(igal)%lumq(6) > 0.d0) THEN
                       print '(G13.5, $)', gal(igal)%lumq(6)
                    ELSE
                       print '(G13.5, $)', 0.d0
                    ENDIF ! 43
                    print '(3G13.5, $)', gal(igal)%taust * param%th_yr, & ! 44
                         gal(igal)%Lpeak, gal(igal)%dMstar_burst * param%munit ! 45,46

                    IF(param%iNLyC /= 0) THEN
                       print '(G13.5, $)', &
                            1.10626d-11 * (lumtmp(param%iNLyC) &
                                           + lumtmp(param%iNLyC+param%nwave))
                    ELSEIF(param%iNLyC_r /= 0) THEN
                       print '(G13.5, $)', &
                            1.10626d-11 * (lumtmp(param%iNLyC_r) &
                                           + lumtmp(param%iNLyC_r+param%nwave))
                    ELSE
                       print '(G13.5, $)', -99.d0
                    ENDIF
                    IF(param%wdim_band) THEN
                       DO iband = 1, param%nwave
                          IF((param%iwave(iband) >= 50 .and. param%iwave(iband) <= 58) &
                               .or. (param%iwave(iband) >= NWAVE_ALL+50 .and. &
                                     param%iwave(iband) <= NWAVE_ALL+58)) THEN
                             ! NLyC, L{1216,1400,1500,1600,1700,2800,4861,6563}
                             print '(2G13.5, $)', &
                                  lumtmp(iband)   + lumtmp(iband+param%nwave), &
                                  lumtmp_d(iband) + lumtmp_d(iband+param%nwave)
                          ELSE
                             print '(2F11.4, $)', mag(iband), mag_d(iband)
                          ENDIF
                       ENDDO
                    ELSE
                       DO iband = 1, param%nwave
                          print '(2F11.4, $)', mag(iband), mag_d(iband)
                       ENDDO
                    ENDIF
                    print *
                 ENDIF
              ENDIF

              IF(mor_d == 1 .or. mor_d == 2) THEN ! E or S0
                 b_or_d = 1 ! bulge-dominated
                 Vdomi  = gal(igal)%Vbulge; Vmino = gal(igal)%Vdisk
                 rdomi  = rb0; rmino = rd0
              ELSE ! S
                 b_or_d = 2 ! disk-dominated
                 Vdomi  = gal(igal)%Vdisk; Vmino = gal(igal)%Vbulge
                 rdomi  = rd0; rmino = rb0
              ENDIF

              ! --- lf(nz)%n(iband,mor_d,bin): mor_d=1,2,3 for E,S0,S LFs
              DO iband = 1, param%nwave
                 bin = int(anint((mag(iband) - lf(nz)%base) * lf(nz)%invstep))
                 IF(bin <= lf(nz)%Nbin .and. bin >= 1 .and. mor_d /= 0) THEN
                    lf(nz)%n(iband, mor_d, bin) &
                         = lf(nz)%n(iband, mor_d, bin) + inv_V
                    IF(gal(igal)%flag_burst == 1) &
                         lf(nz)%n_brst(iband, mor_d, bin) &
                         = lf(nz)%n_brst(iband, mor_d, bin) + inv_V
                 ENDIF
              ENDDO


              ! --- lf_d(nz)%n(i,j,k): j=1,2,3 for E,S0,S LFs
              muSB = -100.d0
              IF(param%SB) THEN ! surface brightness limit is applied
                 IF(rb0 > 0.d0) THEN
                    IF(rd0 > 0.d0) THEN ! composite of disk and bulge
                       SB%iRratio = int((log10(rb0/rd0) - SB%Rratio_base) &
                                        / SB%Rratio_step)
                       IF(SB%iRratio < 1) SB%iRratio = 1
                       IF(SB%iRratio > SB%N_Rratiop1) SB%iRratio = SB%N_Rratiop1
                       SB%BT  = lumtmp_d(SB%iband) &
                                / (lumtmp_d(SB%iband) &
                                   + lumtmp_d(SB%iband+param%nwave))
                       SB%iBT = int((SB%BT - SB%BT_base) / SB%BT_step)
                       IF(SB%iBT < 1) SB%iBT = 1
                       IF(SB%iBT > SB%N_BTp1) SB%iBT = SB%N_BTp1
                    ELSE ! pure bulge
                       SB%iRratio = SB%N_Rratiop1; SB%iBT = SB%N_BTp1
                    ENDIF
                 ELSE ! pure disk
                    SB%iRratio = 1; SB%iBT = 1
                 ENDIF
                 SB%mag_Pet = mag_d(SB%iband) &
                              + SB%Pet2TotRatio(1, SB%iRratio, SB%iBT)
                              ! M(SB%iband)-5logh --> m_Pet(SB%iband) [AB mag]
                 IF(mor_d == 1 .or. mor_d == 2) THEN ! E or S0
                    SB%r_Pet = rdomi * SB%Pet2TotRatio(2, SB%iRratio, SB%iBT)
                 ELSE ! S
                    SB%r_Pet = rdomi * SB%Pet2TotRatio(3, SB%iRratio, SB%iBT)
                 ENDIF
                 muSB = SB%mag_Pet + 5.d0 * log10(SB%r_Pet) + SB%mu0
                        ! converting Petrosian mag. and radius into
                        !  surface brightness
              ENDIF
              IF(muSB < param%SBlimit) THEN
                 ! surface brightness limit added by Makiya (2014/07/08)
                 DO iband = 1, param%nwave
                    bin = int(anint((mag_d(iband) - lf_d(nz)%base) &
                                    * lf_d(nz)%invstep))
                    IF(bin <= lf_d(nz)%Nbin .and. bin >= 1 .and. mor_d /= 0) THEN
                       lf_d(nz)%n(iband, mor_d, bin) &
                            = lf_d(nz)%n(iband, mor_d, bin) + inv_V
                       IF(gal(igal)%flag_burst == 1) &
                            lf_d(nz)%n_brst(iband, mor_d, bin) &
                            = lf_d(nz)%n_brst(iband, mor_d, bin) + inv_V
                    ENDIF
                 ENDDO
              ENDIF

              temp = 1.d0 / (lumtmp(param%iVband_r) &
                             + lumtmp(param%iVband_r+param%nwave))
              zg = (gal(igal)%MZb + gal(igal)%MZd) * temp
              zt = (gal(igal)%Mtb + gal(igal)%Mtd) * temp * param%th
!              zt = gal(igal)%Mtb * param%th / lumtmp(param%iVband_r)




   !!$    ===== QSO ======
              DO iband = 1, param%nwaveq
                 IF(gal(igal)%lumq(iband) > 0.d0) THEN
                    IF(iband == 1) THEN ! B-band
                       mag_q(iband) = - 2.5d0 * log10(gal(igal)%lumq(iband)) &
                            + 88.62d0 + const%corr
                       ! AB mag. and considering effective wavelength is 4400 \AA
                       ! (same as Marconi+04).
                    ELSE IF(iband == 2 .or. iband == 3) THEN ! UV (1450\AA)
                       mag_q(iband) = mag_q(1) + const%B2UV
                    ELSE ! X-ray, bolometric, and Eddington ratio
                       mag_q(iband) = log10(gal(igal)%lumq(iband))
                    ENDIF
                 ELSE
                    IF(iband == 1 .or. iband == 2 .or. iband == 3) THEN ! B- and UV-band
                       mag_q(iband) = 128.d0
                    ELSE IF(iband == 7 .or. iband == 8 .or. iband == 9) THEN
                       mag_q(iband) = -99.d0
                    ELSE ! X-ray and bolometric
                       mag_q(iband) = 0.d0
                    ENDIF
                 ENDIF

                 bin = int(anint((mag_q(iband) - lf_q(nz)%base) &
                                 * lf_q(nz)%invstep))
                 IF(bin <= lf_q(nz)%Nbin .and. bin >= 1) THEN
                    IF(iband == 3) THEN ! observed UV LFs
                       fobs = ObsFracAGN_UV(gal(igal)%lumq(6), param%zsp1)
                       lf_q(nz)%n(iband, 1, bin) &
                            = lf_q(nz)%n(iband, 1, bin) + fobs * inv_V
                    ELSE
                       lf_q(nz)%n(iband, 1, bin) &
                            = lf_q(nz)%n(iband, 1, bin) + inv_V
                    ENDIF
                 ENDIF
              ENDDO
!!$----------------------------------------------------------------------------
              ! --- if one changes the following expressions, the one should modify
              !      the caption of this output file written around LL749--757 in
              !      "mani_nugc.f90"
!              IF(param%run_type /= 3) THEN
              IF(param%run_type /= 3 .and. &
                 param%iSDSS_rpband_r /= 0 .and. mag_d(param%iSDSS_rpband_r) < - 15.d0) THEN
!              IF(param%run_type /= 3 .and. &
!                   (param%iGALEX_FUVband_r /= 0 .and. mag_d(param%iGALEX_FUVband_r) < - 20.d0)) THEN
!              IF(param%run_type /= 3 .and. &
!                   (param%i2MASS_Ksband_r /= 0 .and. mag_d(param%i2MASS_Ksband_r) < - 19.5d0)) THEN
!             IF(param%run_type /= 3 .and. (gal(igal)%Mstarb + gal(igal)%Mstard >= 1.d-6) &
!                .and. gal(igal)%SFR > 0.6d0) THEN
                 write(iout, '(2I10, 7I2, I8, 4I16, 25G13.5, $)') &
                      igal, gal(igal)%id_cgal, & ! 1,2
                      mor, mor_d, gal(igal)%flag_c, & ! 3--5
                      gal(igal)%flag_burst, gal(igal)%flag_merger, & ! 6,7
                      gal(igal)%flag_di, gal(igal)%flag_ccut, & ! 8,9
                      iforest, ihost, & ! 10,11
                      gal(igal)%mpi, mrgt%mpi(ihost), & ! 12,13
                      mrgt%f_des(ihost), gal(igal)%BT, & ! 14,15
                      gal(igal)%Mstarb  * param%munit, & ! 16
                      gal(igal)%Mstard  * param%munit, & ! 17
                      gal(igal)%Mcoolb  * param%munit, & ! 18
                      gal(igal)%Mcoold  * param%munit, & ! 19
                      mrgt%mhalo(ihost) * param%munit, & ! 20
                      gal(igal)%Mhalo   * param%munit, & ! 21
                      gal(igal)%z_col, & ! 22
                      gal(gal(igal)%id_cgal)%Vcent, gal(igal)%Vcent, & ! 23,24
                      gal(igal)%MZcb * param%munit, & ! 25
                      gal(igal)%MZcd * param%munit, & ! 26
                      zg, (const%toutGyr - zt), & ! 27,28
                      rb0, gal(igal)%Vbulge, rd0, gal(igal)%Vdisk, & ! 29--32
                      mag_d(param%iBband_r) + 5.d0*log10(rdomi) + 38.568d0, & ! 33
                      4.62d+5 * rdomi * Square(Vdomi) &
                         * 10.d0**(0.4d0*(mag_d(param%iBband_r)-5.48d0)) &
                         / param%h, & ! 34
                      4.62d+5 * rdomi * Square(Vdomi) / param%h, & ! 35
                      gal(igal)%Mbh * param%munit, & ! 36
                      gal(igal)%MZc_rem * param%munit, & ! 37
                      gal(igal)%SFR, gal(igal)%mSFR  ! 38,39
                 Zmassb = -99.d0; Zmassd = -99.d0
                 IF(gal(igal)%Mstarb > 0.d0) &
                      Zmassb = gal(igal)%Zmassb / gal(igal)%Mstarb
                 IF(gal(igal)%Mstard > 0.d0) &
                      Zmassd = gal(igal)%Zmassd / gal(igal)%Mstard
                 write(iout, '(3G13.5, $)') & ! 40--42
                      Zmassb, Zmassd, &
                      (gal(igal)%Zmassb + gal(igal)%Zmassd) / Mssat
                 IF(gal(igal)%flag_burst == 1 &
                      .and. gal(igal)%lumq(6) > 0.d0) THEN
                   write(iout, '(G13.5, $)') gal(igal)%lumq(6)
                 ELSE
                   write(iout, '(G13.5, $)') 0.d0
                 ENDIF ! 43
                 write(iout, '(3G13.5, $)') gal(igal)%taust * param%th_yr, & ! 44
                      gal(igal)%Lpeak, gal(igal)%dMstar_burst * param%munit ! 45,46

                 IF(param%iNLyC /= 0) THEN
                    write(iout, '(G13.5, $)') &
                         1.10626d-11 * (lumtmp(param%iNLyC) &
                                        + lumtmp(param%iNLyC+param%nwave))
                 ELSEIF(param%iNLyC_r /= 0) THEN
                    write(iout, '(G13.5, $)') &
                         1.10626d-11 * (lumtmp(param%iNLyC_r) &
                                        + lumtmp(param%iNLyC_r+param%nwave))
                 ENDIF
                 DO iband = 1, param%nwave
                    IF((param%iwave(iband) >= 50 .and. param%iwave(iband) <= 58) &
                         .or. (param%iwave(iband) >= NWAVE_ALL+50 .and. &
                               param%iwave(iband) <= NWAVE_ALL+58)) THEN
                      ! NLyC, L{1216,1400,1500,1600,1700,2800,4861,6563}
                       write(iout, '(2G13.5, $)') &
                            lumtmp(iband)   + lumtmp(iband+param%nwave), &
                            lumtmp_d(iband) + lumtmp_d(iband+param%nwave)
                    ELSE
                       write(iout, '(2F11.4, $)') mag(iband), mag_d(iband)
                    ENDIF
                 ENDDO
                write(iout, *)
              ENDIF

!             IF(param%run_type /= 3 .and. ((gal(igal)%Mstarb + gal(igal)%Mstard)>=1.d-6) &
!                .and. gal(igal)%SFR > 0.6d0) THEN
!                 write(iout, '(I2, 5G13.5, $)') &
!                      gal(igal)%flag_burst, &
!                      gal(igal)%Mstarb  * param%munit, & ! 16
!                      gal(igal)%Mstard  * param%munit, & ! 17
!                      gal(igal)%SFR, &
!                      mag_d(param%iBband_r), mag_d(param%iGALEX_FUVband_r)
!                 write(iout, *)
!              ENDIF
!
!             IF(param%run_type /= 3 .and. gal(igal)%flag_c ==1) THEN
!                 write(iout, '(2I2, 2G13.5, $)') &
!                      gal(igal)%flag_ccut, mor_d,&
!                      gal(igal)%Mbh  * param%munit, &
!                      gal(igal)%Mhalo  * param%munit
!                 write(iout, *)
!              ENDIF


!!$----------------------------------------------------------------------------
              ! --- write to *7.dat, AGN properties
              ! --- if one change the following expression,
              ! --- the one should modify the caption of this output file
              ! --- written around LL1332-1346 in "main_nugc.f90".
              IF(param%run_type == 2) THEN ! all redshift
                 iout_q = ionum + NFILE * nz_end + nnode * (nz_end + nz - 1) + inode + 1
              ELSE ! other
                 iout_q = ionum + NFILE + nnode + inode + 1
              ENDIF


              IF(param%run_type /= 3 .and. gal(igal)%lumq(4) > 1.d+41) THEN
!                 IF(gal(igal)%flag_burst == 1) THEN
                    write(iout_q, '(I10, 5I2, 30G13.5)') &
                         igal, gal(igal)%flag_burst, gal(igal)%flag_di, & ! 1,2,3
                         gal(igal)%flag_merger, mor, gal(igal)%flag_c,  & ! 4,5,6
                         max(gal(igal)%Mstard, 1.d-20) * param%munit, & ! 7
                         max(gal(igal)%Mstarb, 1.d-20) * param%munit, & ! 8
                         max(gal(igal)%Mcoold,1.d-20) * param%munit, & ! 9 reside in disk
                         max(gal(igal)%Mcoolb,1.d-20) * param%munit, & ! 10 reside in bulge
                         gal(igal)%Mbh * param%munit, & ! 11
                         max(gal(igal)%MZcd, 1.d-20) * param%munit, & ! 12
                         max(gal(igal)%MZcb, 1.d-20) * param%munit, & ! 13
                         gal(igal)%Mhalo * param%munit, & ! 14
                         max(gal(igal)%Vbulge,1.d-20), max(gal(igal)%Vdisk,1.d-20), & ! 15,16
                         max(gal(igal)%rbulge,1.d-20), max(gal(igal)%rdisk,1.d-20), & ! 17,18
                         mag(param%iBband_r), mag_d(param%iBband_r),   & ! 19,20
                         gal(igal)%SFR, gal(igal)%mSFR, & ! 21,22
                         max(gal(igal)%agn(1),1.d-20), max(gal(igal)%agn(2),1.d-20), & ! 23,24
                         max(gal(igal)%agn(3),1.d-20) * param%munit, & ! 25
                         max(gal(igal)%agn(4),1.d-20), & ! 26
                         max(gal(igal)%agn(5),1.d-20), gal(igal)%agn(6), & ! 27,28
                         max(gal(igal)%agn(7),1.d-20), max(gal(igal)%agn(10),1.d-20), & ! 29,30
                         gal(igal)%agn(8), gal(igal)%agn(9), & ! 31,32
                         gal(igal)%lumq(4), gal(igal)%lumq(6), & ! 33, 34
                         gal(igal)%lumq(7), gal(igal)%lumq(8) ! 35, 36
 !                ENDIF
              ENDIF

              !!$ ==== mass functions ====
              ! mf(:)%n(i,j,k)
              ! --- i:MFtype(1:Mstar^bulge, 2:Mstar^disk, 3:Mcold, 4:MBH,
              !              5:Mhot, 6:Mstar^bulge+MBH, 7:Mstar^disk+Mcold,
              !              8:Mstar^bulge+Mstar^disk+MBH+Mcold,
              !              9:Mstar^bulge+Mstar^disk, 10:hydrogen gas, 11:Mhalo),
              !     j:mor(1:E,2:S0,3:S,4:all), k:mass
              DO j = 1, mf(nz)%N1
                 x = 0.d0
                 IF(j ==  1) THEN ! 1: bulge star mass
                    x = gal(igal)%Mstarb
                 ELSEIF(j ==  2) THEN ! 2: disk star mass
                    x = gal(igal)%Mstard
                 ELSEIF(j ==  3) THEN ! 3: cold gas mass
                    x = gal(igal)%Mcoold + gal(igal)%Mcoolb
                 ELSEIF(j ==  4) THEN ! 4: SMBH mass
                    x = gal(igal)%Mbh
                 ELSEIF(j ==  5) THEN ! 5: hot gas mass
                    x = gal(igal)%Mhot
                 ELSEIF(j ==  6) THEN ! 6: bulge mass
                    Mbulge = gal(igal)%Mstarb + gal(igal)%Mbh + gal(igal)%Mcoolb
                    x = Mbulge
                 ELSEIF(j ==  7) THEN ! 7: disk mass
                    Mdisk = gal(igal)%Mstard + gal(igal)%Mcoold
                    x = Mdisk
                 ELSEIF(j ==  8) THEN ! 8: galaxy mass
                    galmass = Mbulge + Mdisk
                    x = galmass
                 ELSEIF(j ==  9) THEN ! 9: stellar mass
                    x = Mssat
                 ELSEIF(j == 10) THEN ! 10: Hydrogen gas mass
                    M_H = (gal(igal)%Mcoold + gal(igal)%Mcoolb) * const%fracH
                          ! Hydrogen gas mass in cold gas
                    x = M_H
                 ELSEIF(j == 11) THEN ! 11: halo mass
                    IF(gal(igal)%flag_c == 1) & ! count only central galaxy
                         ! "gal(igal)%Mhalo" is identical to
                         !   "mrgt%mhalo(gal(i)%IDhost)" for central
                         x = gal(igal)%Mhalo
                 ENDIF
                 IF(x > EPS) THEN
                    bin = int(anint((log10(x) + mf(nz)%base) * mf(nz)%invstep))
                    IF(bin <= mf(nz)%Nbin .and. bin >= 1 .and. mor_d /= 0) THEN
                       mf(nz)%n(j, mor_d, bin) = mf(nz)%n(j, mor_d, bin) + inv_V
                       IF(gal(igal)%flag_burst == 1) &
                            mf(nz)%n_brst(j, mor_d, bin) &
                            = mf(nz)%n_brst(j, mor_d, bin) + inv_V
                    ENDIF
                 ENDIF
              ENDDO


              !!$ ==== M_HI / L_B as a function of M_B-5logh ===
              IF(param%ML) THEN
                 MassLumi = const%ML * (gal(igal)%Mcoold + gal(igal)%Mcoolb) &
                      / (lumtmp(param%iBband_r) &
                         + lumtmp(param%iBband_r+param%nwave))
                 bin = int(anint((mag_d(param%iBband_r) - ml(nz)%base) &
                                 * ml(nz)%invstep))
                 IF(bin <= ml(nz)%Nbin .and. bin >= 1 .and. mor_d /= 0) THEN
                    ! for mean and sigma
                    ml(nz)%n(1, mor_d, bin) = ml(nz)%n(1, mor_d, bin) + inv_V
                    ml(nz)%n(1, 4,     bin) = ml(nz)%n(1, 4,     bin) + inv_V
                    IF(MassLumi > 0.d0) THEN
                       temp = inv_V * MassLumi
                       ml(nz)%xn(1, mor_d, bin) = ml(nz)%xn(1, mor_d, bin) + temp
                       ml(nz)%xn(1, 4,     bin) = ml(nz)%xn(1, 4,     bin) + temp
                       temp = temp * MassLumi
                       ml(nz)%xxn(1, mor_d, bin) = ml(nz)%xxn(1, mor_d, bin) &
                                                   + temp
                       ml(nz)%xxn(1, 4,     bin) = ml(nz)%xxn(1, 4,     bin) &
                                                   + temp
                    ENDIF

                    ! for median and quartiles
                    ml(nz)%n(2, mor_d, bin) = ml(nz)%n(2, mor_d, bin) + 1.d0
                    ml(nz)%n(2, 4,     bin) = ml(nz)%n(2, 4,     bin) + 1.d0
                    j = int(ml(nz)%n(2, mor_d, bin)); k = int(ml(nz)%n(2, 4, bin))
                    ml(nz)%x(j, mor_d, bin) = MassLumi
                    ml(nz)%x(k, 4, bin)     = MassLumi
                 ENDIF
              ENDIF

              !--- for M_BH - M_bulge relation
              IF(nz == 1 .and. (param%Mbh .or. param%run_type == 3)) THEN
                 IF(gal(igal)%Mstarb > 0.d0) THEN
                    bin = int(anint((log10(gal(igal)%Mstarb*param%munit) &
                                    - MbhMbulge%base) * MbhMbulge%invstep))
                    IF(bin <= MbhMbulge%Nbin .and. bin >= 1) THEN
                       ! for mean and sigma
                       MbhMbulge%n(1, mor_d, bin) &
                            = MbhMbulge%n(1, mor_d, bin) + 1.d0
                       MbhMbulge%n(1, 4,     bin) &
                            = MbhMbulge%n(1, 4,     bin) + 1.d0
                       x = gal(igal)%Mbh * param%munit
                       IF(x > 0.d0) THEN
                          temp = log10(x)
                          MbhMbulge%xn(1, mor_d, bin) = &
                               MbhMbulge%xn(1, mor_d, bin) + temp
                          MbhMbulge%xn(1, 4,     bin) = &
                               MbhMbulge%xn(1, 4,     bin) + temp
                          temp = temp**2.d0 
                          MbhMbulge%xxn(1, mor_d, bin) = &
                               MbhMbulge%xxn(1, mor_d, bin) + temp
                          MbhMbulge%xxn(1, 4,     bin) = &
                               MbhMbulge%xxn(1, 4,     bin) + temp
                       ENDIF
                    ENDIF
                 ENDIF
              ENDIF

              IF(nz == 1 .and. param%run_type == 3) THEN
                 bin = int(anint((log10(gal(igal)%Vdisk) &
                                 - DiskScale%base) * DiskScale%invstep))
                 IF(bin <= DiskScale%Nbin .and. bin >= 1) THEN
                    ! for mean and sigma
                    DiskScale%n(1, mor_d, bin) &
                         = DiskScale%n(1, mor_d, bin) + 1.d0
                    DiskScale%n(1, 4,     bin) &
                         = DiskScale%n(1, 4,     bin) + 1.d0
                    x = rd0 / param%h  
                    IF(x > 0.d0) THEN
                       temp = log10(x)
                       DiskScale%xn(1, mor_d, bin) = &
                            DiskScale%xn(1, mor_d, bin) + temp
                       DiskScale%xn(1, 4,     bin) = &
                            DiskScale%xn(1, 4,     bin) + temp
                       temp = temp**2.d0 
                       DiskScale%xxn(1, mor_d, bin) = &
                            DiskScale%xxn(1, mor_d, bin) + temp
                       DiskScale%xxn(1, 4,     bin) = &
                            DiskScale%xxn(1, 4,     bin) + temp
                    ENDIF
                 ENDIF

                 bin = int(anint((mag_d(param%i2Mass_Ksband_r) &
                                 - SphScale%base) * SphScale%invstep))
                 IF(nz == 1 .and. bin <= SphScale%Nbin .and. bin >= 1 &
                    .and. gal(igal)%Mstarb > 0.d0) THEN
                    ! for mean and sigma
                    IF(mor_d /= 3) THEN
                       SphScale%n(1, mor_d, bin) &
                            = SphScale%n(1, mor_d, bin) + 1.d0
                       SphScale%n(1, 3, bin) &
                            = SphScale%n(1, 3, bin) + 1.d0
                    ENDIF
                    SphScale%n(1, 4,     bin) &
                          = SphScale%n(1, 4,     bin) + 1.d0

                    x = 2.0 * log10(gal(igal)%Vbulge) + log10(rb0) 
                    temp = x 
                    IF(mor_d /= 3) THEN
                       SphScale%xn(1, mor_d, bin) = &
                            SphScale%xn(1, mor_d, bin) + temp
                       SphScale%xn(1, 3, bin) = &
                            SphScale%xn(1, 3, bin) + temp
                    ENDIF
                    SphScale%xn(1, 4,     bin) = &
                         SphScale%xn(1, 4,     bin) + temp
                    temp = temp**2.d0 
                    IF(mor_d /= 3) THEN
                       SphScale%xxn(1, mor_d, bin) = &
                            SphScale%xxn(1, mor_d, bin) + temp
                       SphScale%xxn(1, 3, bin) = &
                            SphScale%xxn(1, 3, bin) + temp
                    ENDIF
                    SphScale%xxn(1, 4,     bin) = &
                         SphScale%xxn(1, 4,     bin) + temp
                  ENDIF
              ENDIF

              ! --- for LAE related
              IF(param%LAE) &
                   call MainCalcForLAE(igal, nz, iforest, mrgt%mhalo(ihost), &
                                       inv_V, rb0, rd0, mor, mor_d, Mssat)

              ! --- for SFH related
              IF(param%SFH) THEN
                 call SubstReffTauV(igal, mor_d, rb0, rd0, &
                                    param%tauV0 * (gal(igal)%MZcd &
                                                   + gal(igal)%MZcb)/const%Zsun)
                 call WriteSFH(igal, nz, iforest, inv_V, rb0, rd0, mor, mor_d, &
                               Mssat)
                 IF(param%run_type == 1) &
                      call CalcCSFH(igal, gal(igal)%flag_burst, inv_V, &
                                    mrgt%mhalo(ihost))
              ENDIF
           ENDIF IF_bright ! mag(param%iBband_r) <= 100.


           ! --- for SFH related
           IF(param%SFH) &
                call SubstReffTauV(igal, mor_d, rd0, rb0, &
                                   param%tauV0 * (gal(igal)%MZcd &
                                                  + gal(igal)%MZcb)) ! ???

           ! --- for SN related
           IF(param%SN) THEN
              IF(Mssat*param%munit >= paramSN%thMstar) THEN
                 call CalAndWriteSNrate(igal, mor_d, gal(igal)%flag_burst, &
                                        inv_V, Mssat*param%munit, reff(igal), &
                                        tauV_SFH(igal), gal(igal)%Tmass, &
                                        gal(igal)%mSFR, &
                                        mrgt%mhalo(ihost) * param%munit, mag_d)
              ENDIF
           ENDIF
        ENDIF ! not newly formed haloes
     ENDIF IF_empty ! Mssat >= 1.d-13 .and. totlum > 0.d0
  ENDDO DO_igal ! igal = 1, endhalo
