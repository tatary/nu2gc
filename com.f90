!!$ ========================================================================
SUBROUTINE opt(lumtmp, lumtmp_d, MZctmp, re)
   !!$     Computing dust extinction assuming slab dust model
   !!$     lumtmp_d = fesc * lumtmp w/ fesc = (1 - exp(-tau)) / tau,
   !!$       only for disk stars
   use global_var
   implicit none

   DOUBLE PRECISION, INTENT(IN)    :: lumtmp(param%tnw)
   DOUBLE PRECISION, INTENT(INOUT) :: lumtmp_d(param%tnw)
   DOUBLE PRECISION, INTENT(IN)    :: MZctmp ! [10^14 Msun]
   DOUBLE PRECISION, INTENT(IN)    :: re     ! [kpc] not in [kpc/h]
   INTEGER :: i, j
   DOUBLE PRECISION :: tau_dust, tauV, fesc
   DOUBLE PRECISION :: Square, f_slab, f_sand ! functions

   tauV = param%tauV0 * (MZctmp / const%Zsun) / Square(re)
   tauV = tauV / param%zsp1**param%alpha_tau
   DO i = 1, param%nwave
      tau_dust = tauV * xi(i)
      fesc     = f_slab(tau_dust) ! escape fraction in slab dust geometry

      j = i + param%nwave
      lumtmp_d(i) = lumtmp(i)        ! spheroid (w/o extinction)
      lumtmp_d(j) = fesc * lumtmp(j) ! disk
   ENDDO
END SUBROUTINE opt
!!$ ========================================================================
SUBROUTINE opt2(lumtmp, lumtmp_d, MZctmpdisk, MZctmpbulge, redisk, rebulge)
   !!$     Computing dust extinction assuming screen dust model
   !!$     lumtmp_d = fesc * lumtmp w/ fesc = exp(-tau),
   !!$       for both of disk and bulge stars in case of starbursts
   use global_var
   implicit none

   DOUBLE PRECISION, INTENT(IN)    :: lumtmp(param%tnw)
   DOUBLE PRECISION, INTENT(INOUT) :: lumtmp_d(param%tnw)
   DOUBLE PRECISION, INTENT(IN)    :: MZctmpdisk, MZctmpbulge ! [10^14 Msun]
   DOUBLE PRECISION, INTENT(IN)    :: redisk, rebulge     ! [kpc] not in [kpc/h]
   INTEGER :: i, j
   DOUBLE PRECISION :: tau_dustdisk, tau_dustbulge
   DOUBLE PRECISION :: tauVdisk, tauVbulge, fescdisk, fescbulge
   DOUBLE PRECISION :: Square, f_scr, f_slab, f_sand ! functions

   ! for bulge !
   tauVbulge = param%tauV0 * (MZctmpbulge / const%Zsun) / Square(rebulge)
   tauVbulge = tauVbulge / param%zsp1**param%alpha_tau

   ! for disk !
   tauVdisk = param%tauV0 * (MZctmpdisk / const%Zsun) / Square(redisk)
   tauVdisk = tauVdisk / param%zsp1**param%alpha_tau

   DO i = 1, param%nwave
      j = i + param%nwave

      ! for bulge !
      tau_dustbulge = tauVbulge * xi(i)
      fescbulge     = f_slab(tau_dustbulge) ! escape fraction in slab geometry
      lumtmp_d(i)   = fescbulge * lumtmp(i) ! spheroid
      ! for disk !
      tau_dustdisk = tauVdisk * xi(i)
      fescdisk     = f_slab(tau_dustdisk)
      lumtmp_d(j)  = fescdisk * lumtmp(j)  ! disk
   ENDDO
END SUBROUTINE opt2
!!$ ========================================================================
SUBROUTINE star(iforest, endhalo, end_step, ionum, nz)
   use global_var; use LAErelated; use SFHrelated; use CLCOOLrelated
   implicit none
   INTEGER, INTENT(IN) :: iforest, ionum, nz
   INTEGER, INTENT(INOUT) :: endhalo, end_step
   DOUBLE PRECISION, PARAMETER :: EPS = 1.d-14 
   DOUBLE PRECISION, PARAMETER :: EPS_Z = 1.d-18 
   ! --- structure for temporally used quantities
   TYPE temporal_values
      INTEGER :: num_g
      DOUBLE PRECISION :: mhalo_sum, Morg, Mhot, MZh, Mhotorg, MZhorg, rdisk
      DOUBLE PRECISION, ALLOCATABLE :: lum(:)
      DOUBLE PRECISION, ALLOCATABLE :: lumq(:)
      DOUBLE PRECISION, ALLOCATABLE :: agn(:)
      DOUBLE PRECISION, ALLOCATABLE :: mem(:)
   END TYPE temporal_values
   TYPE(temporal_values) :: tmp
   ! --- structure for merger related quantities
   TYPE merger_related
      DOUBLE PRECISION :: Msb, Mcb, MZcb
      DOUBLE PRECISION :: MZb, Mtb ! bulge for major (MZb & Mtb)
      DOUBLE PRECISION :: Msd, Mcd, MZcd
      DOUBLE PRECISION :: MZd, Mtd ! bulge for major (MZb & Mtb)
      DOUBLE PRECISION, ALLOCATABLE :: lum(:)
      DOUBLE PRECISION :: Mbh ! BH
      DOUBLE PRECISION :: Tmassb, Tlumb, Zmassb
      DOUBLE PRECISION :: Tmassd, Tlumd, Zmassd
   END TYPE merger_related
   TYPE(merger_related) :: merger
   INTEGER :: nm, nm2, i, j, ig1, ig2, il, is
   INTEGER :: me, iprog, iprog_max, igal, igal_pre, itime, ihalo
   INTEGER :: c_gal, num_g_tot, num_g_me
   INTEGER :: flag_out, flagburst
   INTEGER :: flagdh,  flagdi, flagmerger
   INTEGER :: b_or_q ! 1:starburst, 2:quiescent
   INTEGER :: ier
   INTEGER, ALLOCATABLE :: mgt(:)
   CHARACTER(LEN=50) :: cerr_a = '# star: fail to allocate:'
   CHARACTER(LEN=50) :: cerr_d = '# star: fail to deallocate:'
   DOUBLE PRECISION :: zp1, zp1pc, t_zp1pc, t_zsp1
   DOUBLE PRECISION :: Mme, Mcentral, Msatellite, taumrg, Vsat2d, Vsat2b
   DOUBLE PRECISION :: MZd0, MZb0, Mtb0, Mtd0
   DOUBLE PRECISION :: taucol, prob, x
   DOUBLE PRECISION :: Mb1, Md1, Mc1, Mb2, Md2, Mc2, Mbh1, Mbh2
   DOUBLE PRECISION :: massb
   DOUBLE PRECISION :: y1, y1b, y2
   DOUBLE PRECISION :: e1d, e1b, e1rem, e2d, e2b, xy
   DOUBLE PRECISION :: mz, mt
   DOUBLE PRECISION :: dMch, dMZhc, dMchp, Vccool, Mcool0, dMchn, dMcool
   DOUBLE PRECISION :: xcool, fd_merge
   DOUBLE PRECISION :: Macc, zi, tidal, Mini, taudyn, tburst, Rcol, lr
   DOUBLE PRECISION :: temp
   DOUBLE PRECISION, ALLOCATABLE :: mhalo_sum(:)
   DOUBLE PRECISION :: Vgrav, Vmax, r_di
   DOUBLE PRECISION :: f_dig, f_dig_rem, f_dis, f_dis_rem
                       ! for DI (Shirakata 2016/Jun/27)
   ! --- functions
   REAL             :: ran1, gasdev
   DOUBLE PRECISION :: dens, taustar, z2t, Square, Cube, CalCST
   DOUBLE PRECISION :: CircVel, CalTidal, CalTaumrg
   DOUBLE PRECISION :: CalVbulgeMerger, CalVbulgeDI
   DOUBLE PRECISION :: CalRbulgeMerger, CalRbulgeDI, CalRdisk

   ! --- for Okamoto+08 UV feedback
   INTEGER :: flagO08_decrease, flagO08_output = 0
   INTEGER :: iO08_p, iO08_n0, iO08_n1, iO08_n2, iO08_n3, iO08_n4
   DOUBLE PRECISION :: Mc_O08, Mbar_all, Mgal_all, Mcool_all, Mratio
   DOUBLE PRECISION :: Mc_Okamoto08, fb_Okamoto08 ! function

   ! --- added by Makiya (2015/Aug/19), for fgas dependent fdiss
   DOUBLE PRECISION :: fgas, fgas_DI

   ! --- added by Shirakata (2016/Oct/12) for merger of galaxies
   DOUBLE PRECISION :: fgas_disk, fvr, fvr_rem, metallicity, fdisk
   DOUBLE PRECISION :: Mvr_star, Mvr_gas, Mvr_gasZ
   DOUBLE PRECISION :: Mvrdyn_star, Mvrdyn_gas, Mvrdyn_gasZ
   DOUBLE PRECISION :: M1, M1d, M2, M2d
   DOUBLE PRECISION :: Msat_star, Msat_gas, Msat_gasZ
   DOUBLE PRECISION :: fgas_sat, fdisk_sat, fsat_vr, fsat_vr_rem
   DOUBLE PRECISION :: RgRd, M2M1, GM2M1 ! added 2017/Jul/08



   allocate(mgt(mrgt%num_galarray),    stat=ier); call CheckIerr(ier, &
      trim(cerr_a)//' mgt')
   allocate(mhalo_sum(0:mrgt%num_tot), stat=ier); call CheckIerr(ier, &
      trim(cerr_a)//' mhalo_sum')
   allocate(tmp%lum(param%nwave),      stat=ier); call CheckIerr(ier, &
      trim(cerr_a)//' tmp%lum')
   allocate(tmp%lumq(param%nwaveq),    stat=ier); call CheckIerr(ier, &
      trim(cerr_a)//' tmp%lumq')
   allocate(tmp%agn(10),               stat=ier); call CheckIerr(ier, &
      trim(cerr_a)//' tmp%agn')
   allocate(merger%lum(param%tnw),     stat=ier); call CheckIerr(ier, &
      trim(cerr_a)//' merger%lum')
   allocate(tmp%mem(10),             stat=ier); call CheckIerr(ier, &
      trim(cerr_a)//' tmp%mem')

   num_g_tot = 0
   t_zsp1 = z2t(param%zsp1) ! age of the universe at param%zsp1
   !  added by MARK (2013/Aug/30)
   call SearchMostMassiveProgenitorHalo

   IF(param%run_type /= 3 .and. nnode == 1) &
      print '(A)', '#         z_out :zp1ar(i)-1-->  zp1pc-1 : Ngal_tot: tlife[Myr]'


   !!$  ==================================================
   !!$    ***   Main Loop for Time Snapshot Begins   ***
   !!$  ==================================================
   !!$  ID_itime: DO itime = 1, mrgp%num_step ! original
   ID_itime: DO itime = 1, mrgp%num_step
      IF(param%UVfb == 2) THEN ! Okamoto+08 UV feedback mode
         iO08_p  = 0 ! # of halo w/ Macc >= 0
         iO08_n0 = 0 ! # of halo w/ Macc < 0
         iO08_n1 = 0 ! # of halo w/ Macc+tmp%Mhotorg < 0
         iO08_n2 = 0 ! # of halo w/ Macc+tmp%Mhotorg+tmp%Mhot < 0
         iO08_n3 = 0 ! # of halo w/ Macc+tmp%Mhotorg+tmp%Mhot+Mcool_all < 0
         iO08_n4 = 0 ! # of halo w/ Macc+Mbar_all < 0
      ENDIF

      ! +------------------------------------------------------+
      ! |   calculating life time "tlife" and set "flag_out"   |
      ! +------------------------------------------------------+
      IF(itime == param%izout) THEN
         flag_out = 1; zp1pc = param%zsp1
      ELSE
         flag_out = 0; zp1pc = mrgp%zp1ar(itime+1)
      ENDIF
      t_zp1pc = z2t(zp1pc) ! added by MARK (2013/Aug/30)
      zp1 = mrgp%zp1ar(itime)
      !!$     tlife = mrgp%time(itime+1) - mrgp%time(itime)
      !!$             ! almost identical to "tlife" calculated via the following eq.
      tlife     = t_zp1pc - z2t(zp1)  ! tlife [hubble time]
      tlife_yr  = tlife * param%th_yr ! tlife [yr]
      tlife_sec = tlife_yr * const%yr2sec ! tlife [sec]
      IF(param%run_type /= 3 .and. nnode < 2) THEN
         print '(2(A, I2), 3(A,F8.3), A, I7, A, F8.3)', &
            '(', itime, '/', param%izout, ')', real(param%zsp1)-1.d0, &
            ' : ', real(mrgp%zp1ar(itime))-1.d0, ' --> ', real(zp1pc)-1.d0, &
            ' : ', num_g_tot, ' : ', tlife_yr * 1.d-6
      ENDIF
      IF(flag_out == 1 .and. param%run_type /= 3) &
         print '(A, I7)', '# num_total = ', num_total
      ! +------------------------------------------------------+

      ! --- initialize the quantities related to the # of galaxies
      !       for this timestep "itime"
      num_g_tot = 0 ! total # of galaxies at this time-step
      num_total = 0 ! cumulative # of galaxies w/ non-zero baryon in all halos
                    !   from the halo "ihalo" = 1 at this timestep "itime"
      num_cum   = 0 ! cumulative # of galaxies in all halos from the halo
                    !   "ihalo" = 1 at this timestep "itime"

      IF(param%UVfb == 2) THEN ! Okamoto+08 UV feedback mode
         Mc_O08 = Mc_Okamoto08(zp1pc - 1.d0) ! [10^14 Msun]
      ENDIF

      ! +--------------------------------------------------------------+
      ! | Create cooling function table (added by Shirakata 18/Aug/22) |
      ! +--------------------------------------------------------------+
      IF (param%CoolFN == 2) & ! Okamoto-cloudy only
         call LoadCLCoolingRedshift(zp1,clcool_table_z)

      !!$  ======================================================================
      !!$    ***   Loop for All Halos at the Time Snapshot "itime" Begings  ***
      !!$  ======================================================================
      ID_ihalo: DO ihalo = 1, mrgp%num_now(itime)
         me = mrgp%st_halo(itime) + ihalo - 1 ! (ihalo,itime) --> me
         ! corresponding consecutive number for the halo "ihalo"
         !  at the time of "itime" in all snapshots

         call InitializeLocalValiablesPerHalo
         ! initialize flagburst(=0), mhalo_sum(me), Vsat2, dMhot,
         !  tmp%xxx, merger%xxx


         ! +--------------------------------------------------------------+
         ! |   calculate total # of galaxies contained in the halo "me"   |
         ! |    and total progenitor halo mass "mhalo_sum(me)"            |
         ! +--------------------------------------------------------------+
         iprog_max = mrgt%f_prg(me) + mrgt%numb(me) - 1
                     ! the maximum ID of the "me"'s progenitor halos
         DO iprog = mrgt%f_prg(me), iprog_max ! loop for all progenitors of "me"
            tmp%num_g     = tmp%num_g     + mrgt%num_g(iprog) ! tot. # of gals
            mhalo_sum(me) = mhalo_sum(me) + mrgt%mhalo(iprog) ! tot. mass of
                                                              !  prog. halos
         ENDDO
         ! +--------------------------------------------------------------+


         ! +-------------------------------------------+
         ! |   Set Central Galaxy of Halo "me" begins  |
         ! +-------------------------------------------+
         IF(tmp%num_g == 0) THEN ! --- the halo "me" is empty:
                                 !     --- for all halos at "itime" = 1
                                 !     --- newly-collapsed halos at "itime" >= 2
            call InitializeCentral
         ELSE ! --- non-zero # of galaxies in the halo "me": for "itime" >= 2
              ! --- replace struct.: "gal_prev(igal+num_cum)" --> "gal(igal)"
            DO igal = 1, tmp%num_g
               igal_pre  = igal + num_cum
               gal(igal) = gal_prev(igal_pre)
               IF(param%LAE) allgal(igal) = allgal_prev(igal_pre)

               ! initialize SFR of this timestep "itime"
               IF(flag_out == 0) THEN
                  gal(igal)%SFR = 0.d0
                  IF(tlife > const%t0) gal(igal)%mSFR = 0.d0
                  ! added by MARK (2014/Nov/04)

                  gal(igal)%dMstar_burst = 0.d0 ! added by MARK (2017/Mar/16)

                  IF(param%LAE) allgal(igal)%SFR(:) = 0.d0
               ENDIF

               IF(param%UVfb == 2) THEN ! Okamoto+08 UV feedback mode
                  Mbar_all  = Mbar_all  + BaryonMass(igal)
                              ! tot. baryon mass (incl. hot gas) in the halo "me"
                  Mgal_all  = Mgal_all  + GalMass(igal)
                  Mcool_all = Mcool_all + ColdMass(igal)
               ENDIF
            ENDDO

            IF(tmp%num_g == 1) THEN
               igal = 1
               gal(igal)%flag_c = 1; c_gal = igal
            ELSE ! There are more than two galaxies in the same halo
               c_gal = 0
               temp = 0.d0
               DO igal = 1, tmp%num_g
                  IF(gal(igal)%flag_c == 1 .and. &
                       gal(igal)%IDhost == mrgt%c_halo(me) .and. &
                       BaryonMass(igal) > 0.d0) c_gal = igal
                  ! --- if a galaxy "igal" is the central galaxy (i.e.,
                  !      gal(igal)%flag_c = 1) resided in the most massive
                  !      prog. halo of "me" (gal(igal)%IDhost =
                  !      mrgt%c_halo(me)), the galaxy is assgined to the central
                  !      galaxy of the halo "me"
               ENDDO

               IF(c_gal == 0) THEN
                  ! in the case of the absence of central
                  ! e.g., collapse epoch of the most massive prog. is too late
                  !     to form stars (Vcirc < param%Vlow or Mhalo << Mc_O08)
                  IF(mrgt%c2_halo(me) /= -1) THEN
                     ! the case that the 2nd most massive prog. is present
                     ! in this case, the cent. gal. in the 2nd most massive
                     !   prog. is assigned to the cent. gal. in the halo "me"
                     DO igal = 1, tmp%num_g
                        IF(BaryonMass(igal) > 0.d0 .and. &
                             gal(igal)%flag_c == 1 .and. &
                             gal(igal)%IDhost == mrgt%c2_halo(me)) c_gal = igal
                     ENDDO
                  ELSE
                     ! the case that the 2nd most massive prog. is absent
                     ! in this case, the most massive gal. is assigned to
                     !   the cent. gal.
                     temp = 0.d0
                     DO igal = 1, tmp%num_g
                        IF(BaryonMass(igal) > temp) THEN
                           temp = BaryonMass(igal); c_gal = igal
                        ENDIF
                     ENDDO
                  ENDIF
               ENDIF

               IF(c_gal == 0) THEN
                  print '(2(A, I2), A, I6, 2(A, I6), A)', &
                       '# central is absence!!: (iforest, itime, me) = (', &
                       iforest, ', ', itime, ', ', me, '), (c_halo, '//&
                       'c2_halo) = (', mrgt%c_halo(me), ', ', &
                       mrgt%c2_halo(me), ')'
                  DO igal = 1, tmp%num_g
                     print'(A, I2, A, I6, 3(A, G10.3), A)', &
                          '  --- (igal, IDhost) = (', igal, ', ', &
                          gal(igal)%IDhost, '): (Mstar, Mcool, Mbar_all) '//&
                          '= (', StellarMass(igal) * param%munit, &
                          ', ', ColdMass(igal) * param%munit, &
                          ', ', BaryonMass(igal)*param%munit, ')'
                  ENDDO
               ENDIF

               DO igal = 1, tmp%num_g
                  gal(igal)%flag_c = 0
               ENDDO
               gal(c_gal)%flag_c = 1
            ENDIF

            ! --- loop for all galaxies in the halo "me" begins ---
            ID_igal1: DO igal = 1, tmp%num_g
               IF(flag_out == 0) call ResetFlags(igal) ! added by Shirakata(2016/Jul/23)
               IF(gal(igal)%flag_c == 1) THEN
                  Mcentral = GalMass(igal) ! = Mstar(disk+bulge) + Mcool + Mbh
                  ! --- add hot comp. of galaxy "igal" to tmp%{Mhot,MZh}{,org}
                  tmp%Mhot    = tmp%Mhot    + gal(igal)%Mhot
                  tmp%Mhotorg = tmp%Mhotorg + gal(igal)%Mhotorg
                  tmp%MZh     = tmp%MZh     + gal(igal)%MZh
                  tmp%MZhorg  = tmp%MZhorg  + gal(igal)%MZhorg
                  ! --- erase hot comp. of galaxy "igal"
                  gal(igal)%Mhot    = 0.d0; gal(igal)%MZh    = 0.d0
                  gal(igal)%Mhotorg = 0.d0; gal(igal)%MZhorg = 0.d0

                  gal(igal)%Mhalo = Mme; gal(igal)%clps = zp1
                  gal(igal)%mpi   = mrgt%mpi(me)
               ELSE ! satellite galaxies of the halo "me"
                  gal(igal)%Mreheat = 0

                  ! --- add hot comp. of galaxy "igal" to tmp%{Mhot,MZh}org
                  tmp%Mhotorg = tmp%Mhotorg + gal(igal)%Mhotorg + gal(igal)%Mhot
                  tmp%MZhorg  = tmp%MZhorg  + gal(igal)%MZhorg  + gal(igal)%MZh

                  ! --- erase hot components of galaxy "igal"
                  gal(igal)%Mhot    = 0.d0; gal(igal)%MZh    = 0.d0
                  gal(igal)%Mhotorg = 0.d0; gal(igal)%MZhorg = 0.d0

                  gal(igal)%Morg = 0.d0

                  ! --- if the galaxy "igal" enters the halo "me" for the first
                  !       time, gal(igal)%Telapse and gal(igal)%taumrg are
                  !       initialized as 0.d0
                  IF(gal(igal)%IDhost /= mrgt%c_halo(me)) THEN
                     gal(igal)%Telapse = 0.d0
                     gal(igal)%taumrg  = 0.d0
                  ENDIF
               ENDIF
               gal(igal)%IDprog = gal(igal)%IDhost ! prog. halo ID for "igal"
               gal(igal)%IDhost = me ! host halo ID for "igal"
               gal(igal)%hfinal = mrgt%f_des(me)
            ENDDO ID_igal1
            ! --- loop for all galaxies in the halo "me" ends ---

            num_cum = num_cum + tmp%num_g
            ! cumulative # of galaxies in all halos from the halo
            !   "ihalo" = 1 at this timestep "itime"

            ! --- calculating contrib. to hot baryon comp. from accreted DM ---
            IF(mhalo_sum(me) > 0.) THEN
               IF(param%UVfb == 2) THEN ! Okamoto+08 UV feedback mode
                  ! --------------------------------------------------------------
                  !   In this mode, baryon fraction in a given halo of mass
                  !     "Mme" depends on z, Mme, and z_reion
                  !   Newly accreting baryonic mass "Macc" followed by halo
                  !     accretion is determined by the difference of baryon
                  !     masses in the halo and its progeniters "Mbar_all"
                  !   In case of Macc < 0 (only for Mme <~ Mc_O08), the
                  !     baryonic mass in the halo should be decreased
                  !     (evaporate via UV feedback and escape from the halo)
                  !     --- in such case, the baryonic mass is decreased in the
                  !           following order:
                  !         (1)tmp%Mhotorg: hot gas in the halo for
                  !              cooling calc.
                  !         (2)tmp%Mhot: hot gas in the halo produced by SF
                  !         (3)gal(igal)%Mcool: cold gas mass of all gals in
                  !              the halo
                  ! --------------------------------------------------------------
                  flagO08_decrease = 0
                  IF(mhalo_sum(me) > Mme) THEN ! halo mass decreases from
                     !   the tot. mass of progenitors
                     temp = Mme / mhalo_sum(me) ! ~ unity
                     tmp%Mhot    = tmp%Mhot    * temp
                     tmp%Mhotorg = tmp%Mhotorg * temp
                     tmp%MZh     = tmp%MZh     * temp
                     tmp%MZhorg  = tmp%MZhorg  * temp
                     Mbar_all = (Mbar_all - Mgal_all) * temp + Mgal_all
                     flagO08_decrease = 1
                  ENDIF

                  Macc = Mme * fb_Okamoto08(Mme, zp1pc-1.d0, Mc_O08) - Mbar_all
                  ! accreting mass of baryon
                  IF(Macc > 0.d0) THEN
                     IF(flagO08_output == 1) iO08_p = iO08_p + 1
                     tmp%Mhotorg = tmp%Mhotorg + Macc
                  ELSE
                     IF(flagO08_output == 1) THEN
                        iO08_n0 = iO08_n0 + 1
                        temp = Macc + tmp%Mhotorg
                        IF(temp < 0.d0) THEN
                           iO08_n1 = iO08_n1 + 1
                           temp = temp + tmp%Mhot
                           IF(temp < 0.d0) THEN
                              iO08_n2 = iO08_n2 + 1
                              temp = temp + Mcool_all
                              IF(temp < 0.d0) THEN
                                 iO08_n3 = iO08_n3 + 1
                                 temp = Macc + Mbar_all
                                 IF(temp < 0.d0) iO08_n4 = iO08_n4 + 1
                              ENDIF
                           ENDIF
                        ENDIF
                     ENDIF
                     temp = tmp%Mhotorg + Macc
                     IF(temp > 0.d0) THEN
                        Mratio = temp / tmp%Mhotorg ! less than unity
                        tmp%Mhotorg = tmp%Mhotorg * Mratio
                        tmp%MZhorg  = tmp%MZhorg  * Mratio
                     ELSE
                        tmp%Mhotorg = temp ! negative value
                     ENDIF
                  ENDIF

                  IF(tmp%Mhotorg < 0.d0) THEN
                     ! in case that accreting mass is negative and tmp%Mhotorg
                     !   is not sufficiently large to supply the negative
                     !   accreting mass
                     ! residual negative accreting mass is assumed to be
                     !   supplied from tmp%Mhot
                     temp = tmp%Mhotorg + tmp%Mhot
                     IF(temp > 0.d0) THEN
                        Mratio = temp / tmp%Mhot ! less than unity
                        tmp%Mhot = tmp%Mhot * Mratio
                        tmp%MZh  = tmp%MZh  * Mratio
                     ELSE
                        tmp%Mhot = 0.d0; tmp%MZh = 0.d0
                     ENDIF
                     tmp%Mhotorg = 0.d0; tmp%MZhorg = 0.d0
                  ENDIF
               ELSEIF(param%UVfb == 1) THEN ! past fiducial mode
                  IF(mhalo_sum(me) < Mme) THEN
                     ! significant amount of DM is accreted
                     Macc        = Mme - mhalo_sum(me) ! accreting mass of DM
                     tmp%Mhotorg = tmp%Mhotorg + Macc * param%bar_rat
                  ELSE
                     temp = Mme / mhalo_sum(me) ! ~ unity
                     tmp%Mhot = tmp%Mhot * temp; tmp%Mhotorg = tmp%Mhotorg * temp
                     tmp%MZh  = tmp%MZh  * temp; tmp%MZhorg  = tmp%MZhorg  * temp
                  ENDIF
               ENDIF
            ENDIF
         ENDIF
         ! +-----------------------------------------+
         ! |   Set Central Galaxy of Halo "me" ends  |
         ! +-----------------------------------------+
         num_g_tot = num_g_tot + tmp%num_g
         ! total # of galaxies at this time-step


         ! --- initialize hot components of the central galaxy "c_gal" in "me"
         gal(c_gal)%Mhot    = tmp%Mhot;    gal(c_gal)%MZh    = tmp%MZh
         gal(c_gal)%Mhotorg = tmp%Mhotorg; gal(c_gal)%MZhorg = tmp%MZhorg


         ! +----------------------------------------------------------------+
         ! |   at "itime" = "param%izout", jump to the sentence labeled as  |
         ! |     2001 and skip all of the following calculations            |
         ! |   note that the sentence labeled as 2001 is still inside the   |
         ! |     DO loop for "ihalo"                                        |
         ! +----------------------------------------------------------------+
         IF(flag_out == 1) goto 2001
         ! +----------------------------------------------------------------+



!!$           *** DISTINCTION OF HALO MAJOR MERGER ***
         IF(gal(c_gal)%Mreheat <= param%freheat * Mme) THEN
            flagdh = 1

            gal(c_gal)%Mreheat = Mme ! renewal "Mreheat" as "Mme" = Mhalo(me)
            gal(c_gal)%density = dens(zp1)
            gal(c_gal)%Vc      = CircVel(Mme, gal(c_gal)%density)
            gal(c_gal)%Vcent   = gal(c_gal)%Vc
            gal(c_gal)%Telapse = tlife
            gal(c_gal)%z_form1 = zp1
            gal(c_gal)%Morg    = Mme
            gal(c_gal)%Mhot    = gal(c_gal)%Mhot + gal(c_gal)%Mhotorg
            gal(c_gal)%MZh     = gal(c_gal)%MZh  + gal(c_gal)%MZhorg
            gal(c_gal)%Mhotorg = 0.d0; gal(c_gal)%MZhorg = 0.d0
            gal(c_gal)%cst     = CalCST(zp1, Mme)
         ELSE
            flagdh = 0

            gal(c_gal)%Telapse = gal(c_gal)%Telapse + tlife
         ENDIF
         Vccool = CircVel(Mme, gal(c_gal)%density)

         ! revised by Shirakata (2016/08/05)
         call CoolingProcess
         gal(c_gal)%flag_ccut = CoolingCutoff_key


         ! ===================================================================
         !    choose the satellites which reside in halo as satellites for
         !      longer time than the timescale of dynamical friction
         !    and calculate the merger contribution
         ! ===================================================================
         nm = 1 ! # of galaxies surviving the dynamical friction
         ID_igal2: DO igal = 1, tmp%num_g
            call InitializeLocalValiablesPerGal
            IF(param%DI) THEN   ! <<DISK INSTABILITY>> (Shirakata 2017/Feb/28)
               ! --- If galactic disk become unstable,
               !       a part of the disk is destroyed.
               ! --- Destructed fraction is determined with f_di.
               ! --- Disk instability occurs only spiral galaxies.
               ! --- At first collapsed time, this process is skipped.
               IF(gal(igal)%Mstard > 0.d0 .and. gal(igal)%Mcoold > 0.d0 &
                  .and. gal(igal)%Mbh > 0.d0) THEN
                  r_di = CalRdisk(gal(igal)%Vc, gal(igal)%density)
                  IF(xcool > 0.d0 .and. gal(igal)%flag_c == 1) THEN
                     r_di = xcool * r_di
                  ENDIF
                  r_di = max(r_di, gal(igal)%rdisk)
                  r_di = 0.05d0 * r_di / log(2.d0)

                  Vmax = gal(igal)%Vmax * gal(igal)%Vmax
                  IF(gal(igal)%Mstarb > 0.d0) THEN
                     IF(r_di <= gal(igal)%rbulge) THEN
                        Vmax = Vmax &
                               + 2.d0 / 3.d0 * gal(igal)%Vbulge * gal(igal)%Vbulge
                     ELSE IF(r_di > gal(igal)%rbulge .and. gal(igal)%rbulge > 0.d0) THEN
                        Vmax = Vmax &
                               + 4.45d8 * param%h * BulgeMass(igal) / r_di 
                     ENDIF
                  ENDIF
                  Vmax = sqrt(Vmax)

                  IF(r_di > 0.d0) THEN
                     Vgrav = 3.012d+11 * sqrt(const%G * param%h * DiskMass(igal) / r_di)
                             ! km/s
                     IF(Vmax / Vgrav < param%em) THEN
                        gal(igal)%flag_di = 1; gal(igal)%flag_burst = 1
                        IF(gal(igal)%flag_c == 1) flagdi = 1; flagburst = 1
                     ENDIF
                  ENDIF
               ENDIF
            ENDIF

            IF(igal /= c_gal .and. GalMass(igal) > 0.d0) THEN
               ! prescriptions only for satellites
               IF(gal(igal)%taumrg <= 0.d0 .or. flagdh == 1) THEN
                  ! comment out by Makiya, 2015/11/01
                  ! tidal = CalTidal(gal(c_gal)%Morg/gal(igal)%Mhalo, &
                  !                  gal(c_gal)%Vc/gal(igal)%Vc)
                  ! Msatellite = gal(igal)%Mhalo * tidal

                  ! IF(Msatellite <= 0.) goto 4000
                  !   ! skip the calculations related to dyn. fric.

                  ! taumrg = param%fmerge &
                  !          * CalTaumrg(gal(c_gal)%Morg/Msatellite, gal(c_gal)%density)
                  !          ! --- timescale of dynamical friction

                  ! When we use Jiang+ model for the time scale of dynamical friction,
                  !    we should use the 'non-stripped' satellite halo mass
                  taumrg = param%fmerge &
                           * CalTaumrg(Mme/gal(igal)%Mhalo, Mme, dens(zp1), zp1)
                  gal(igal)%taumrg = taumrg
                  ! --- timescale of dynamical friction
               ELSE
                  taumrg = gal(igal)%taumrg
               ENDIF
               gal(igal)%Telapse = gal(igal)%Telapse + tlife

               IF(gal(igal)%Telapse > taumrg .and. GalMass(c_gal) > 0.d0) THEN ! <<MERGER>>
                  gal(c_gal)%n_merge = gal(c_gal)%n_merge + 1
                  Vsat2d = Vsat2d + KinEnergyDisk(igal)
                  Vsat2b = Vsat2b + KinEnergyBulge(igal)
                  ! total energy, sum[MV^2]
                  flagmerger = 1; flagburst = 1
                  fgas_disk = 0.d0; fdisk = 0.d0
                  metallicity = 0.d0
                  M1  = GalMass(c_gal); M2 = GalMass(igal)
                  M1d = DiskMass(c_gal)

                  IF(M2 > M1) THEN
                     ! --- satellite is larger than central
                     ! --- then, disk mass is renewed based on satellite's
                     !       disk property (Shirakata 2017/02/09)
                     gal(c_gal)%flag_cfalse = 1
                     M2d = DiskMass(igal); fgas_sat = 0.d0
                     Msat_star = 0.d0; Msat_gas     = 0.d0
                     fsat_vr   = 0.d0; fsat_vr_rem = 1.d0
                     IF(M2d > 0.d0) THEN
                        fgas_sat  = gal(igal)%Mcoold / M2d
                        fdisk_sat = M2d / M2
                        IF(gal(igal)%Mcoold > 0.d0) &
                           metallicity = gal(igal)%MZcd / gal(igal)%Mcoold

                        M2M1  = M1/M2
                        IF(gal(c_gal)%M2M1_max < M2M1) &
                           gal(c_gal)%M2M1_max = M2M1
                        gal(c_gal)%M2M1_av = gal(c_gal)%M2M1_av + M2M1
                        GM2M1 = 2.d0 * M2M1 / (1.d0 + M2M1)
                        RgRd  = (1.d0-fgas_sat) * fdisk_sat * 1.2d0 * GM2M1 
                        fsat_vr     = GM2M1 
                        fsat_vr_rem = 1.d0 - fsat_vr

                        Msat_gas  = gal(igal)%Mcoold * (1.d0 - (1.d0 + RgRd) &
                                    * exp(-RgRd))
                        Msat_gasZ = metallicity * Msat_gas
                        Msat_star = GM2M1 * gal(igal)%Mstard
                     ENDIF

                     merger%Msb  = merger%Msb  + gal(igal)%Mstarb + Msat_star
                     merger%Msd  = merger%Msd  + gal(igal)%Mstard - Msat_star
                     merger%Mcb  = merger%Mcb  + gal(igal)%Mcoolb + Msat_gas
                     merger%Mcd  = merger%Mcd  + gal(igal)%Mcoold - Msat_gas
                     merger%MZcb = merger%MZcb + gal(igal)%MZcb + Msat_gasZ
                     merger%MZcd = merger%MZcd + gal(igal)%MZcd - Msat_gasZ
                     merger%MZb  = merger%MZb + gal(igal)%MZb &
                                   + fsat_vr * gal(igal)%MZd
                     merger%MZd  = merger%MZd + fsat_vr_rem * gal(igal)%MZd
                     merger%Mtb  = merger%Mtb + gal(igal)%Mtb &
                                   + fsat_vr * gal(igal)%Mtd
                     merger%Mtd  = merger%Mtd + fsat_vr_rem * gal(igal)%Mtd
                     merger%lum(1:param%nwave) &
                          = merger%lum(1:param%nwave) &
                            + gal(igal)%lumg(1:param%nwave) &
                            + fsat_vr * gal(igal)%lumg(param%nwp1:param%tnw)
                     merger%lum(param%nwp1:param%tnw) &
                          = merger%lum(param%nwp1:param%tnw) &
                            + fsat_vr_rem * gal(igal)%lumg(param%nwp1:param%tnw)
                     merger%Mbh    = merger%Mbh    + gal(igal)%Mbh
                     merger%Tmassb = merger%Tmassb + gal(igal)%Tmass
                     merger%Tlumb  = merger%Tlumb  + gal(igal)%Tlum_b &
                                     + fsat_vr * gal(igal)%Tlum_d
                     merger%Tlumd  = merger%Tlumd  + fsat_vr_rem * gal(igal)%Tlum_d

                     ! --- added by MARK (2014/Nov/02)
                     merger%Zmassb = merger%Zmassb + gal(igal)%Zmassb &
                                     + fsat_vr * gal(igal)%Zmassd
                     merger%Zmassd = merger%Zmassd + fsat_vr_rem * gal(igal)%Zmassd
                  ELSE
                     Mvr_star = 0.d0; Mvr_gas = 0.d0
                     IF(M1d > 0.d0) THEN
                        fgas_disk = gal(c_gal)%Mcoold / M1d 
                        fdisk     = M1d / M1
                        IF(gal(c_gal)%Mcoold > 0.d0) &
                             metallicity = gal(c_gal)%MZcd / gal(c_gal)%Mcoold

                        M2M1  = M2/M1
                        IF(gal(c_gal)%M2M1_max < M2M1) &
                           gal(c_gal)%M2M1_max = M2M1
                        gal(c_gal)%M2M1_av = gal(c_gal)%M2M1_av + M2M1
                        GM2M1 = 2.d0 * M2M1 / (1.d0 + M2M1)
                        RgRd  = (1.d0-fgas_disk) *fdisk * 1.2d0 * GM2M1

                        Mvr_gas  = gal(c_gal)%Mcoold * (1.d0 - (1.d0 + RgRd) &
                                   * exp(-RgRd))
                        Mvr_star = GM2M1 * gal(c_gal)%Mstard
                     ENDIF
                     Mvrdyn_gas  = Mvrdyn_gas  + Mvr_gas
                     Mvrdyn_star = Mvrdyn_star + Mvr_star

                     merger%Msb  = merger%Msb  + StellarMass(igal)
                     merger%Mcb  = merger%Mcb  + ColdMass(igal)
                     merger%MZcb = merger%MZcb + ColdMZ(igal)
                     merger%MZb  = merger%MZb  + GalMZ(igal)
                     merger%Mtb  = merger%Mtb  + GalMt(igal)
                     merger%lum(1:param%nwave) = merger%lum(1:param%nwave) &
                                                 + gal(igal)%lumg(1:param%nwave) &
                                                 + gal(igal)%lumg(param%nwp1:param%tnw)
                     merger%Mbh    = merger%Mbh    + gal(igal)%Mbh
                     merger%Tmassb = merger%Tmassb + gal(igal)%Tmass
                     merger%Tlumb  = merger%Tlumb  + gal(igal)%Tlum_b &
                                     + gal(igal)%Tlum_d

                     ! --- added by MARK (2014/Nov/02)
                     merger%Zmassb = merger%Zmassb + gal(igal)%Zmassb &
                                     + gal(igal)%Zmassd
                  ENDIF

                  ! --- AGN related
                  gal(c_gal)%lumq(1:9) = 0.d0 
                  gal(c_gal)%agn(1:9)  = 0.d0
                  gal(c_gal)%agn(10)   = gal(c_gal)%agn(10) + gal(igal)%agn(10)
                  gal(igal)%agn(1:10)  = 0.d0

                  ! --- added by MARK (2014/Nov/07)
                  IF(param%SFH) THEN
                     gal(c_gal)%SFH(:,:) = gal(c_gal)%SFH(:,:) + gal(igal)%SFH(:,:)
!$                    gal(c_gal)%SFH0(:,:) = gal(c_gal)%SFH0(:,:) + gal(igal)%SFH0(:,:)
                  ENDIF

                  call EraseGalaxy(igal)

               ELSE ! <<NO MERGER>>
                  mgt(nm) = igal ! the i's galaxy ID surviving dyn. fric. are
                                 !   saved as "mgt(i)"
                  nm = nm + 1 ! # of galaxies surviving dyn. fric.
               ENDIF
            ENDIF
!!$ 4000        continue
         ENDDO ID_igal2

         gal(c_gal)%flag_burst  = flagburst
         gal(c_gal)%flag_merger = flagmerger
         gal(c_gal)%flag_di     = flagdi
         ! ======================================================================


         ! ===============================================================
         !    choose the galaxies which collide with each other via random
         !      collision and calculate the merger contribution
         ! ===============================================================
!!$           *** MERGER BY RANDOM COLLISION ***
         nm  = nm - 1 ! # of galaxies surviving dyn. fric. other than myself
         nm2 = nm ! # of galaxies surviving both dyn. fric. and rand. coll.
         !   other than myself
         ID_igal3: DO igal = 1, nm
            call InitializeLocalValiablesPerGal
            ig1 = mgt(igal) ! ID of igal's galaxy surviving dyn. fric.

            IF(nm2 == 0) THEN ! --- there is no other galaxy other than "ig1"
               goto 3000 ! skip the cals for merger via rand. coll.
            ELSE ! --- there are some galaxies other than "ig1"
               ! --- if "ig1" is empty, skip the cals for merger via rand. coll.
               IF(GalMass(ig1) <= 0.d0) goto 3000

               IF(nm2 == 1) THEN ! --- there is only one galaxy other than "ig1"
                  x = 1.d0; prob = 0.d0 ! make rand. coll. not occur forcibly
               ELSEIF(nm2 >= 2) THEN
                  IF(param%equal_mass_merger) THEN ! --- equal-mass mergers
                     Rcol = dble(nm2)
                     taucol = 1.43d+8 * Cube(Mme) * param%h &
                              * (gal(ig1)%Vc / max(gal(ig1)%Vdisk, gal(ig1)%Vbulge))**4 &
                              / (Cube(gal(c_gal)%Vc) * Square(gal(ig1)%Mhalo) * Rcol)
                     x = dble(ran1(param%idum)); prob = tlife / taucol
                  ELSE ! --- non equal-mass mergers
                     Rcol = 0.d0
                     DO j = 1, nm
                        ig2 = mgt(j)
                        IF(ig2 /= ig1) THEN
                           IF(GalMass(ig2) > 0.d0) THEN
                              lr   = (gal(ig2)%Mhalo / gal(ig1)%Mhalo)**(1.d0/3.d0)
                              temp = 1.d0 + Square(lr)
                              Rcol = Rcol + Square(1.d0+lr) * Square(temp) / 16.d0
                           ENDIF
                        ENDIF
                     ENDDO
                  ENDIF
               ENDIF
            ENDIF

            IF(param%off_collision) THEN
               x = 1.d0; prob = 0.d0 ! make rand. coll. not occur forcibly
            ENDIF

!!$     === stripping ===
!!$               IF(gal(me)%Vc > 500.d0) THEN
!!$                  gal(me)%Mhot = gal(me)%Mhot + gal(ig1)%Mcool
!!$                  gal(ig1)%Mcool = 0.d0
!!$               ENDIF

1000        continue

            IF(x <= prob) THEN ! --- random collision occurs
               nm2 = nm2 - 1 ! # of gals surviving both dyn. fric. and
               !   rand. coll. other than myself

               ! --- choose the galaxy ID of collision partner randomly
               j = int(dble(ran1(param%idum)) * dble(nm - 1)) + 1
               IF(j >= igal) j = j + 1 ! ????
               ig2 = mgt(j) ! "ig2" is the gal. ID of collision partner of "ig1"

               Mb1  = gal(ig1)%Mstarb; Mb2  = gal(ig2)%Mstarb ! bulge mass
               Md1  = gal(ig1)%Mstard; Md2  = gal(ig2)%Mstard ! disk mass
               Mc1  = ColdMass(ig1);   Mc2 = ColdMass(ig2)      ! cold gas mass
               Mbh1 = gal(ig1)%Mbh;    Mbh2 = gal(ig2)%Mbh    ! SMBH mass
               y1   = Mb1 + Md1 + Mc1 + Mbh1; y2 = Mb2 + Md2 + Mc2 + Mbh2

               ! --- gas mass fraction, added by Makiya (2015/Aug/09)
               ! --- cold gas (in disk + bulge; Shirakata 2016/Jul/24)
               ! fgas = (Mc1 + Mc2) &
               !        /(Mb1 + Mb2 + Md1 + Md2 + Mc1 + Mc2 + Mbh1 + Mbh2)

               ! --- if "ig2" is empty, make this rand. coll. not occur
               !       forcibly
               IF(y2 == 0.d0) THEN
                  x = 1.d0; prob = 0.d0; nm2 = nm2 + 1
                  goto 1000
               ENDIF

               ! --- determine which galaxy is larger in mass: ID of larger
               !      galaxy is saved as "il" and that of smaller one is as "is"
               ! --- mass ratio of "is"/"il" is saved as "xy"
               IF(y1 > y2) THEN
                  il = ig1; is = ig2; xy = y2 / y1
               ELSE
                  il = ig2; is = ig1; xy = y1 / y2
                  y1   = Mb2 + Md2 + Mc2 + Mbh2; y1 = Mb1 + Md1 + Mc1 + Mbh1
               ENDIF
               y1b  = BulgeMass(il) 
               gal(il)%flag_burst  = 1; b_or_q = 1
               gal(il)%flag_merger = 1
               gal(il)%n_merge = gal(il)%n_merge + 1

               Mini  = y1 + y2 ! initial total mass of the two merging gals
               e1b = KinEnergyBulge(il); e1d = KinEnergyDisk(il)
               e2b = KinEnergyBulge(is); e2d = KinEnergyDisk(is)
               M1  = GalMass(il);   M2 = GalMass(is)
               M1d   = DiskMass(il)

               fdisk     = M1d / M1
               fgas_disk = 0.d0
               M2M1  = M2 / M1
               IF(gal(il)%M2M1_max < M2M1) &
                  gal(il)%M2M1_max = M2M1
               gal(il)%M2M1_av = gal(il)%M2M1_av + M2M1
               GM2M1 = 2.d0 * M2M1 / (1.d0 + M2M1)

               Mvr_star = 0.d0; Mvr_gas = 0.d0
               IF(M1d > 0.d0) THEN 
                  fgas_disk   = gal(il)%Mcoold / M1d
                  RgRd = (1.d0-fgas_disk) *fdisk * 1.2d0 * GM2M1
                  Mvr_gas  = gal(il)%Mcoold * (1.d0 - (1.d0 + RgRd) &
                             * exp(-RgRd))
                  Mvr_star = GM2M1 * gal(il)%Mstard
                  IF(xy > param%fmajor) THEN
                     Mvr_gas = gal(il)%Mcoold
                     Mvr_star = gal(il)%Mstard
                  ENDIF
               ENDIF

               IF(gal(il)%flag_di == 1) THEN
                  call FracDI(gal(il)%Mcoold, gal(il)%Mstard, GalMass(il),f_dig, f_dis, f_dig_rem, f_dis_rem)
                  Mvr_star = Mvr_star + f_dis * gal(il)%Mstard 
                  Mvr_gas  = Mvr_gas + f_dig * gal(il)%Mcoold 
               ENDIF

               Mvr_star = min(Mvr_star, gal(il)%Mstard)
               Mvr_gas  = min(Mvr_gas, gal(il)%Mcoold)

               fvr      = 0.d0; fvr_rem = 1.d0
               IF(gal(il)%Mstard > 0.d0) THEN
                  fvr     = Mvr_star / gal(il)%Mstard
                  fvr_rem = 1.d0 - fvr
               ENDIF

               metallicity = 0.d0
               IF(gal(il)%Mcoold > 0.d0) &
                  metallicity = gal(il)%MZcd / gal(il)%Mcoold
               Mvr_gasZ = Mvr_gas * metallicity

               ! --- calculate bulge state --- !
               ! ----- In all merger event, starburst occur and disk does survive.
               ! ----- Disk gas mass fraction doesn't change by merger.
               Ms0 = gal(il)%Mstarb + Mvr_star + StellarMass(is)
               gal(il)%Mstard = gal(il)%Mstard - Mvr_star

               Mc0 = gal(il)%Mcoolb + Mvr_gas + ColdMass(is)
               gal(il)%Mcoold = gal(il)%Mcoold - Mvr_gas

               MZc0 = gal(il)%MZcb + Mvr_gasZ + ColdMZ(is)
               gal(il)%MZcd = gal(il)%MZcd - Mvr_gasZ

               MZb0 = fvr * gal(il)%MZd + gal(il)%MZb + GalMZ(is)
               gal(il)%MZd = fvr_rem * gal(il)%MZd

               Mtb0 = gal(il)%Mtb + fvr * gal(il)%Mtd + GalMt(is)
               gal(il)%Mtd = fvr_rem * gal(il)%Mtd

               gal(il)%lumg(1:param%nwave) &
                    = gal(il)%lumg(1:param%nwave) &                ! bulge
                      + fvr * gal(il)%lumg(param%nwp1:param%tnw) & ! disk
                      + gal(is)%lumg(1:param%nwave) &              ! bulge
                      + gal(is)%lumg(param%nwp1:param%tnw)         ! disk

               gal(il)%lumg(param%nwp1:param%tnw) &      ! disk
                    = fvr_rem * gal(il)%lumg(param%nwp1:param%tnw)


               ! --- AGN related
               Mbh0 = gal(il)%Mbh + gal(is)%Mbh
               Mbh1 = gal(il)%Mbh
               gal(il)%lumq(1:9) = 0.d0 
               gal(il)%agn(1:9)  = 0.d0
               gal(il)%agn(10)    = gal(il)%agn(10) + gal(is)%agn(10)

               ! --- added by Shirakata
               massb = Ms0 + Mc0 + Mbh0

               ! --- added by MARK (2013/Aug/29)
               Tmass0  = gal(il)%Tmass + gal(is)%Tmass
               Tlum_b0 = gal(il)%Tlum_b + fvr * gal(il)%Tlum_d &
                         + gal(is)%Tlum_b + gal(is)%Tlum_d
               gal(il)%Tlum_d = fvr_rem * gal(il)%Tlum_d

               ! --- added by MARK (2014/Nov/02)
               Zmassb0 = gal(il)%Zmassb + fvr * gal(il)%Zmassd &
                         + gal(is)%Zmassb + gal(is)%Zmassd
               gal(il)%Zmassd = fvr_rem * gal(il)%Zmassd


               ! --- SFH related
               IF(param%SFH) THEN
                  gal(il)%SFH(:,:) = gal(il)%SFH(:,:) + gal(is)%SFH(:,:)
!!$                    gal(il)%SFH0(:,:) = gal(il)%SFH0(:,:) + gal(is)%SFH0(:,:)
               ENDIF

               ! --- added by MARK (2014/Nov/07)
               SFR0  = gal(il)%SFR  + gal(is)%SFR
               mSFR0 = gal(il)%mSFR + gal(is)%mSFR

               call EraseGalaxy(is)

               ! --- added by Shirakata (2017/Oct/05)
               fd_merge = 0.d0
               IF(M1d > 0.d0) &
                  fd_merge = DiskMass(il) / M1d
               gal(il)%rdisk = fd_merge * gal(il)%rdisk
               gal(il)%diskmas = DiskMass(il)

               fgas = Mc0 / Mini
               e1rem = param%Krem * KinEnergyDisk(il)
               e1d  = e1d - e1rem
               gal(il)%Vbulge = CalVbulgeMerger(il, e1d, e1b, e1rem, e2d, e2b, y1b, y1, y2, fgas, massb)
               gal(il)%rbulge = CalRbulgeMerger(il, massb, y1b, gal(il)%Vbulge) ! kpc/h

               taust = 1.d-8
               beta  = CalBeta(gal(il)%Vbulge, b_or_q)

               x = x * tlife
               call star_formation_burst(zp1)
               call RenewalGalaxyBurst(il)

               IF(param%dyn_resp_bulge .or. param%dyn_resp_halo) THEN
                  zi    = Cal_zi(Mini / gal(il)%Mhalo, &
                                 gal(il)%Vc / gal(il)%Vbulge)
                  tidal = CalTidal(gal(c_gal)%Mhalo / gal(il)%Mhalo,&
                                   gal(c_gal)%Vc / gal(il)%Vc)
               ENDIF
               IF(param%dyn_resp_bulge) call DynamicalResponseBulge(il)
               IF(param%dyn_resp_halo)  call DynamicalResponseHalo(il)

               call EstimateMZcRem(il)

               b_or_q = 2 ! quiescent
               Ms0  = gal(il)%Mstard
               MZd0 = gal(il)%MZd
               Mtd0 = gal(il)%Mtd
               Tlum_d0 = gal(il)%Tlum_d
               Zmassd0 = gal(il)%Zmassd

               Mc0  = gal(il)%Mcoold; MZc0 = gal(il)%MZcd
               Mbh0 = gal(il)%Mbh

               ! - calculate normal star formation - !
               IF(gal(il)%rdisk <= 0.d0) THEN
                    gal(il)%rdisk = CalRdisk(gal(il)%Vc, gal(il)%density)
                    gal(il)%Vdisk = gal(il)%Vc
                 ENDIF
                 
               taust = taustar(gal(il)%Vdisk, gal(il)%rdisk, gal(il)%clps, Mc0, MZc0, b_or_q)
               beta  = CalBeta(gal(il)%Vdisk, b_or_q)
               
               IF(Mc0 > 0.d0 .and. MZc0 / Mc0 > 1.d0) print '(A)', '1'

               call star_formation
               call RenewalGalaxyQuiescent(il)


               IF(param%dyn_resp_disk .or. param%dyn_resp_halo) THEN
                  zi = Cal_zi(Mini / gal(il)%Mhalo, &
                              gal(il)%Vc / gal(il)%Vdisk)

               ! comment out by Makiya, 2015/11/01
               !tidal = CalTidal(gal(c_gal)%Mhalo / gal(il)%Mhalo,&
               !                 gal(c_gal)%Vc / gal(il)%Vc)

                  tidal = 1.d0
               ENDIF
               IF(param%dyn_resp_disk) call DynamicalResponseDisk(il)
               IF(param%dyn_resp_halo) call DynamicalResponseHalo(il)

            ELSE ! --- random collision does not occur
               ! --- calculate ordinary star formation (disk formation) for
               !       the satellite galaxy "ig1"
               IF(gal(ig1)%flag_di == 1 .and. gal(ig1)%flag_merger == 0) THEN
                  ! --- disk instability without mergers
                  b_or_q = 1 ! starburst via disc instability
                  call FracDI(gal(ig1)%Mcoold, gal(ig1)%Mstard, GalMass(ig1), &
                              f_dig, f_dis, f_dig_rem, f_dis_rem)

                  y1 = DiskMass(ig1); y2 = BulgeMass(ig1) 
                  Mini = y1 + y2; M1d = DiskMass(ig1)

                  e1d = KinEnergyDisk(ig1) 
                  e1b = KinEnergyBulge(ig1) 

                  Mb1 = gal(ig1)%Mstarb; Md1 = f_dis * gal(ig1)%Mstard 
                  Ms0 = Mb1 + Md1
                  gal(ig1)%Mstard = f_dis_rem * gal(ig1)%Mstard

                  Mc1  = f_dig * gal(ig1)%Mcoold 
                  Mc0  = Mc1 + gal(ig1)%Mcoolb 
                  gal(ig1)%Mcoold = f_dig_rem * gal(ig1)%Mcoold
                  MZc0 = f_dig * gal(ig1)%MZcd + gal(ig1)%MZcb
                  gal(ig1)%MZcd = f_dig_rem * gal(ig1)%MZcd

                  MZb0 = gal(ig1)%MZb + f_dis * gal(ig1)%MZd
                  gal(ig1)%MZd = f_dis_rem * gal(ig1)%MZd
                  Mtb0 = gal(ig1)%Mtb + f_dis * gal(ig1)%Mtd
                  gal(ig1)%Mtd = f_dis_rem * gal(ig1)%Mtd

                  gal(ig1)%lumg(1:param%nwave) &
                       = gal(ig1)%lumg(1:param%nwave) &
                         + f_dis * gal(ig1)%lumg(param%nwp1:param%tnw)
                  gal(ig1)%lumg(param%nwp1:param%tnw) &
                       = f_dis_rem * gal(ig1)%lumg(param%nwp1:param%tnw)

                  ! --- AGN related
                  Mbh0 = gal(ig1)%Mbh
                  Mbh1 = Mbh0
                  gal(ig1)%lumq(1:9) = 0.d0
                  gal(ig1)%agn(1:9)  = 0.d0
                  gal(ig1)%agn(10)   = gal(ig1)%agn(10)

                  Tmass0 = f_dis * gal(ig1)%Tmass
                  Tlum_b0 = gal(ig1)%Tlum_b + f_dis * gal(ig1)%Tlum_d
                  gal(ig1)%Tlum_d = f_dis_rem * gal(ig1)%Tlum_d

                  Zmassb0 = gal(ig1)%Zmassb + f_dis * gal(ig1)%Zmassd
                  gal(ig1)%Zmassd = f_dis_rem * gal(ig1)%Zmassd

                  SFR0 = gal(ig1)%SFR; mSFR0 = gal(ig1)%mSFR

                  ! --- added by Shirakata (2017/Oct/05)
                  fd_merge = 0.d0
                  IF(M1d > 0.d0) &
                     fd_merge = DiskMass(ig1) / M1d
                  gal(ig1)%rdisk = fd_merge * gal(ig1)%rdisk
                  gal(ig1)%diskmas = DiskMass(ig1)


                  fgas_DI = Mc0 / Mini 

                  e1rem = param%Krem * KinEnergyDisk(ig1)
                  e1d = e1d - e1rem 
                  massb = Ms0 + Mc0 + Mbh0

                  gal(ig1)%Vbulge = CalVbulgeDI(e1d, e1b, y1, y2, fgas_DI, massb)
                  gal(ig1)%rbulge = CalRbulgeDI(massb, gal(ig1)%Vbulge) ! kpc/h

                  taust = 1.d-8
                  beta  = CalBeta(gal(ig1)%Vbulge, b_or_q)

                  call star_formation_burst(zp1)
                  x = dble(ran1(param%idum)) * tlife
                  call RenewalGalaxyBurst(ig1)

                  IF(param%dyn_resp_bulge .or. param%dyn_resp_halo) THEN
                     zi = Cal_zi(Mini / gal(ig1)%Mhalo, &
                                 gal(ig1)%Vc / gal(ig1)%Vbulge)
                     tidal = 1.d0
                  ENDIF
                  IF(param%dyn_resp_bulge) call DynamicalResponseBulge(ig1)
                  IF(param%dyn_resp_halo)  call DynamicalResponseHalo(ig1)

                  call EstimateMZcRem(ig1)
               ENDIF ! flag_di = 0 and random collision does not occur

               b_or_q = 2 ! quiescent

               Ms0  = gal(ig1)%Mstard; MZd0 = gal(ig1)%MZd; Mtd0 = gal(ig1)%Mtd
               Mc0  = gal(ig1)%Mcoold; MZc0 = gal(ig1)%MZcd
               Mini = GalMass(ig1)

               Mbh0 = gal(ig1)%Mbh

               ! --- added by MARK (2013/Sep/02)
               Tmass0 = gal(ig1)%Tmass; Tlum_d0 = gal(ig1)%Tlum_d

               ! --- added by MARK (2014/Nov/02)
               Zmassd0 = gal(ig1)%Zmassd

               ! --- added by MARK (2014/Nov/07)
               SFR0 = gal(ig1)%SFR; mSFR0 = gal(ig1)%mSFR

               IF(gal(ig1)%rdisk <= 0.d0) THEN
                  gal(ig1)%rdisk = CalRdisk(gal(ig1)%Vc, gal(ig1)%density)
                  gal(ig1)%Vdisk = gal(ig1)%Vc
               ENDIF

               taust = taustar(gal(ig1)%Vdisk, gal(ig1)%rdisk, gal(ig1)%clps, Mc0, MZc0, b_or_q)
               beta  = CalBeta(gal(ig1)%Vdisk, b_or_q)

               IF(Mc0 > 0.d0 .and. MZc0 / Mc0 > 1.d0) print '(A)', '2'

               call star_formation
               call RenewalGalaxyQuiescent(ig1)

               IF(param%dyn_resp_disk .or. param%dyn_resp_halo) THEN
                  IF(gal(ig1)%Vdisk <= 0.d0) gal(ig1)%Vdisk = gal(ig1)%Vc
                  zi = Cal_zi(Mini / gal(ig1)%Mhalo, &
                              gal(ig1)%Vc / gal(ig1)%Vdisk)

                  ! commnet out by Makiya, 2015/11/01
                  ! tidal = CalTidal(gal(c_gal)%Mhalo / gal(ig1)%Mhalo,&
                  !                  gal(c_gal)%Vc / gal(ig1)%Vc)

                  tidal = 1.d0

               ENDIF
               IF(param%dyn_resp_disk) call DynamicalResponseDisk(ig1)
               IF(param%dyn_resp_halo) call DynamicalResponseHalo(ig1)
            ENDIF
3000        continue
         ENDDO ID_igal3
! ======================================================================



         ! ===============================================================
         ! ||   calculations of star formation for the central galaxy   ||
         ! ||     "c_gal" in "me" in the timestep "itime":              ||
         ! ||   --- (1) starburst via dyn, fric. of satellites          ||
         ! ||   --- (2) disk instability without mergers                ||
         ! ||   --- (3) disk size estimation                            ||
         ! ||   --- (4) ordinary star formation (disk formation)        ||
         ! ===============================================================


         call InitializeLocalValiablesPerGal 

         ! +-------------------------------------------------+
         ! |  (1) starburst via dyn. fric. of satellites     |
         ! +-------------------------------------------------+
         IF(gal(c_gal)%flag_merger == 1) THEN
            b_or_q = 1 ! starburst
            y1 = GalMass(c_gal); y1b = BulgeMass(c_gal)
            y2 = merger%Msb + merger%Mcb + merger%Mbh &
                 + merger%Msd + merger%Mcd
            Mini = y1 + y2; M1d = DiskMass(c_gal)

            e1d = KinEnergyDisk(c_gal); e1b = KinEnergyBulge(c_gal) 
            e2d = Vsat2d;               e2b = Vsat2b 

            IF(gal(c_gal)%flag_di == 1) THEN
                  call FracDI(gal(c_gal)%Mcoold, gal(c_gal)%Mstard, GalMass(c_gal), &
                       f_dig, f_dis, f_dig_rem, f_dis_rem)
               Mvrdyn_star = Mvrdyn_star + f_dis * gal(c_gal)%Mstard
               Mvrdyn_gas = Mvrdyn_gas + f_dig * gal(c_gal)%Mcoold
            ENDIF

            IF(gal(c_gal)%M2M1_max > param%fmajor) THEN
               Mvrdyn_star = gal(c_gal)%Mstard
               Mvrdyn_gas  = gal(c_gal)%Mcoold
            ENDIF

            Mvrdyn_star = min(Mvrdyn_star, gal(c_gal)%Mstard)
            Mvrdyn_gas = min(Mvrdyn_gas, gal(c_gal)%Mcoold)

            fvr = 0.d0; fvr_rem = 1.d0
            IF(gal(c_gal)%Mstard > 0.d0) THEN
               fvr     = Mvrdyn_star / gal(c_gal)%Mstard
               fvr_rem = 1.d0 - fvr
            ENDIF
            metallicity = 0.d0
            IF(gal(c_gal)%Mcoold > 0.d0) &
                 metallicity = gal(c_gal)%MZcd / gal(c_gal)%Mcoold
            Mvrdyn_gasZ = Mvrdyn_gas * metallicity

            Mb1 = gal(c_gal)%Mstarb
            Md1 = Mvrdyn_star 
            Ms0 = Mb1 + Md1 + merger%Msb

            gal(c_gal)%Mstard = fvr_rem * gal(c_gal)%Mstard &
                                + merger%Msd

            Mc1 = gal(c_gal)%Mcoolb + Mvrdyn_gas
            Mc0 = Mc1 + merger%Mcb
            gal(c_gal)%Mcoold = gal(c_gal)%Mcoold - Mvrdyn_gas &
                                + merger%Mcd

            ! --- AGN related
            Mbh1 = gal(c_gal)%Mbh
            Mbh0 = Mbh1 + merger%Mbh
            gal(c_gal)%Mbh = 0.d0
            merger%Mbh = 0.d0

            MZc0 = gal(c_gal)%MZcb + Mvrdyn_gasZ + merger%MZcb
            gal(c_gal)%MZcd = gal(c_gal)%MZcd - Mvrdyn_gasZ &
                              + merger%MZcd

            IF(Mc0 > 0.d0 .and. MZc0 / Mc0 > 1.d0) print '(A)', '3'

            MZb0 = gal(c_gal)%MZb + fvr * gal(c_gal)%MZd + merger%MZb
            gal(c_gal)%MZd = max(fvr_rem * gal(c_gal)%MZb, 0.d0) + merger%MZd
            Mtb0 = gal(c_gal)%Mtb + fvr * gal(c_gal)%Mtd + merger%Mtb
            gal(c_gal)%Mtd = max(fvr_rem * gal(c_gal)%Mtd, 0.d0) + merger%Mtd

            gal(c_gal)%lumg(1:param%nwave) &
                 = gal(c_gal)%lumg(1:param%nwave) &                ! bulge
                   + fvr * gal(c_gal)%lumg(param%nwp1:param%tnw) & ! disk
                   + merger%lum(1:param%nwave)
            gal(c_gal)%lumg(param%nwp1:param%tnw) &
                 = fvr_rem * gal(c_gal)%lumg(param%nwp1:param%tnw) &
                   + merger%lum(param%nwp1:param%tnw)


            ! --- added by Shirakata
            massb = Ms0 + Mc0 + Mbh0

            ! --- added by MARK (2013/Aug/29)
            Tmass0  = gal(c_gal)%Tmass + merger%Tmassb + merger%Tmassd
            Tlum_b0 = gal(c_gal)%Tlum_b + fvr * gal(c_gal)%Tlum_d &
                      + merger%Tlumb
            gal(c_gal)%Tlum_d = fvr_rem * gal(c_gal)%Tlum_d + merger%Tlumd

            ! --- added by MARK (2014/Nov/02)
            Zmassb0 = gal(c_gal)%Zmassb + fvr * gal(c_gal)%Zmassd &
                      + merger%Zmassb
            gal(c_gal)%Zmassd = fvr_rem * gal(c_gal)%Zmassd + merger%Zmassd

            ! --- added by MARK (2014/Nov/08)
            SFR0 = gal(c_gal)%SFR; mSFR0 = gal(c_gal)%mSFR

            ! --- added by Shirakata (2017/Oct/05)
            fd_merge = 0.d0
            IF(M1d > 0.d0) &
               fd_merge = DiskMass(c_gal) / M1d
            gal(c_gal)%rdisk = fd_merge * gal(c_gal)%rdisk
            gal(c_gal)%diskmas = DiskMass(c_gal)

            ! --- gas mass fraction, added by Makiya (2015/Aug/09)
            fgas = Mc0 / Mini

            e1rem = param%Krem * KinEnergyDisk(c_gal)
            e1d  = MAX(0.d0, e1d - e1rem)
            gal(c_gal)%Vbulge = CalVbulgeMerger(c_gal, e1d, e1b, e1rem, e2d, e2b, y1b, y1, y2, fgas, massb)
            gal(c_gal)%rbulge = CalRbulgeMerger(c_gal, massb, y1b, gal(c_gal)%Vbulge) ! kpc/h
            taust = 1.d-8
            beta  = CalBeta(gal(c_gal)%Vbulge, b_or_q)

            call star_formation_burst(zp1)
            x = dble(ran1(param%idum)) * tlife
            call RenewalGalaxyBurst(c_gal)

            IF(param%dyn_resp_bulge .or. param%dyn_resp_halo) THEN
               zi    = Cal_zi(Mini / gal(c_gal)%Mhalo, &
                              gal(c_gal)%Vc / gal(c_gal)%Vbulge)
               tidal = 1.d0
            ENDIF
            IF(param%dyn_resp_bulge) call DynamicalResponseBulge(c_gal)
            IF(param%dyn_resp_halo)  call DynamicalResponseHalo(c_gal)

            call EstimateMZcRem(c_gal)

            ! +----------------------------------------+
            ! |  (2) disk instability without mergers  |
            ! +----------------------------------------+
         ELSE IF(gal(c_gal)%flag_merger == 0 .and. gal(c_gal)%flag_di == 1) THEN ! DI only
            b_or_q = 1 ! starburst via disc instability

            call FracDI(gal(c_gal)%Mcoold, gal(c_gal)%Mstard, GalMass(c_gal), &
                 f_dig, f_dis, f_dig_rem, f_dis_rem)

            y1 = DiskMass(c_gal); y2 = BulgeMass(c_gal)
            Mini = y1 + y2; M1d = DiskMass(c_gal)

            e1d = KinEnergyDisk(c_gal) 
            e1b = KinEnergyBulge(c_gal) 

            Mb1 = gal(c_gal)%Mstarb; Md1 = f_dis * gal(c_gal)%Mstard
            Ms0 = Mb1 + Md1
            gal(c_gal)%Mstard = f_dis_rem * gal(c_gal)%Mstard

            Mc1  = f_dig * gal(c_gal)%Mcoold 
            Mc0  = Mc1 + gal(c_gal)%Mcoolb 
            gal(c_gal)%Mcoold = f_dig_rem * gal(c_gal)%Mcoold
            MZc0 = f_dig * gal(c_gal)%MZcd + gal(c_gal)%MZcb
            gal(c_gal)%MZcd = f_dig_rem * gal(c_gal)%MZcd

            MZb0 = gal(c_gal)%MZb + f_dis * gal(c_gal)%MZd
            gal(c_gal)%MZd = f_dis_rem * gal(c_gal)%MZd
            Mtb0 = gal(c_gal)%Mtb + f_dis * gal(c_gal)%Mtd
            gal(c_gal)%Mtd = f_dis_rem * gal(c_gal)%Mtd

            IF(Mc0 > 0.d0 .and. MZc0 / Mc0 > 1.d0) print '(A)', '4'
            gal(c_gal)%lumg(1:param%nwave) &
                 = gal(c_gal)%lumg(1:param%nwave) &
                   + f_dis * gal(c_gal)%lumg(param%nwp1:param%tnw)
            gal(c_gal)%lumg(param%nwp1:param%tnw) &
                 = f_dis_rem * gal(c_gal)%lumg(param%nwp1:param%tnw)

            ! --- AGN related
            Mbh0 = gal(c_gal)%Mbh
            Mbh1 = Mbh0
            gal(c_gal)%lumq(1:9) = 0.d0 
            gal(c_gal)%agn(1:9)  = 0.d0
            gal(c_gal)%agn(10)   = gal(c_gal)%agn(10)

            Tmass0  = f_dis * gal(c_gal)%Tmass
            Tlum_b0 = gal(c_gal)%Tlum_b + f_dis * gal(c_gal)%Tlum_d
            gal(c_gal)%Tlum_d = f_dis_rem * gal(c_gal)%Tlum_d

            Zmassb0 = gal(c_gal)%Zmassb + f_dis * gal(c_gal)%Zmassd
            gal(c_gal)%Zmassd = f_dis_rem * gal(c_gal)%Zmassd

            SFR0 = gal(c_gal)%SFR; mSFR0 = gal(c_gal)%mSFR

            ! --- added by Shirakata (2017/Oct/05)
            fd_merge = 0.d0
            IF(M1d > 0.d0) &
               fd_merge = DiskMass(c_gal) / M1d
            gal(c_gal)%rdisk = fd_merge * gal(c_gal)%rdisk
            gal(c_gal)%diskmas = DiskMass(c_gal)

            massb = Ms0 + Mc0 + Mbh0
            fgas_DI = Mc0 / Mini 

            e1rem = param%Krem * KinEnergyDisk(c_gal)
            e1d = e1d - e1rem 
            gal(c_gal)%Vbulge = CalVbulgeDI(e1d, e1b, y1, y2, fgas_DI, massb)
            gal(c_gal)%rbulge = CalRbulgeDI(Mc0+Ms0+Mbh0, gal(c_gal)%Vbulge) ! kpc/h

            taust = 1.d-8
            beta  = CalBeta(gal(c_gal)%Vbulge, b_or_q)

            call star_formation_burst(zp1)
            x = dble(ran1(param%idum)) * tlife
            call RenewalGalaxyBurst(c_gal)

            IF(param%dyn_resp_bulge .or. param%dyn_resp_halo) THEN
               zi    = Cal_zi(Mini / gal(c_gal)%Mhalo, &
                              gal(c_gal)%Vc / gal(c_gal)%Vbulge)
               tidal = 1.d0
            ENDIF
            IF(param%dyn_resp_bulge) call DynamicalResponseBulge(c_gal)
            IF(param%dyn_resp_halo)  call DynamicalResponseHalo(c_gal)

            call EstimateMZcRem(c_gal)
         ENDIF


         ! +-----------------------------+
         ! |  (3) disk size estimation   |
         ! +-----------------------------+
         ! --- Criteria for estimation:
         !     (a) When r_disk = 0: assign r_disk and Vdisk
         !     (b) When r_disk > 0: if Vc < param%Vcut & disk mass increases &
         !                          r_disk increases, assign new r_disk and Vdisk
         ! +-----------------------------+
         x = DiskMass(c_gal)
         IF(x > 0.d0) THEN
            IF(gal(c_gal)%rdisk <= 1.d-5) THEN
               gal(c_gal)%Vdisk   = Vccool 
               gal(c_gal)%diskmas = x
               gal(c_gal)%rdisk   = CalRdisk(Vccool, gal(c_gal)%density) ! kpc/h
               IF(xcool > 0.) THEN
                  gal(c_gal)%rdisk = gal(c_gal)%rdisk * xcool
               ENDIF
            ELSE IF(gal(c_gal)%flag_cfalse == 1) THEN
               tmp%rdisk = CalRdisk(Vccool, gal(c_gal)%density)
               IF(xcool > 0.d0) THEN
                  tmp%rdisk = xcool * tmp%rdisk
               ENDIF
               gal(c_gal)%Vdisk   = Vccool
               gal(c_gal)%diskmas = x
               gal(c_gal)%rdisk   = tmp%rdisk
            ELSE IF (x > gal(c_gal)%diskmas) THEN
               tmp%rdisk = CalRdisk(gal(c_gal)%Vc, gal(c_gal)%density) ! kpc/h
               IF(xcool > 0.d0) THEN
                  tmp%rdisk = xcool * tmp%rdisk
               ENDIF
               IF(gal(c_gal)%rdisk < tmp%rdisk) THEN
                  gal(c_gal)%Vdisk   = gal(c_gal)%Vc
                  gal(c_gal)%diskmas = x
                  gal(c_gal)%rdisk   = tmp%rdisk
               ENDIF
            ENDIF
         ENDIF
         ! +-----------------------------+

         ! +-------------------------------------------------+
         ! |  (4) ordinary star formation (disk formation)   |
         ! +-------------------------------------------------+
         b_or_q = 2 ! quiescent

         Mini = GalMass(c_gal) ! Mstard + Mstarb + Mcool + Mbh
         IF(gal(c_gal)%Mcoold <= 0.d0) goto 1001

         Ms0  = gal(c_gal)%Mstard; Mc0  = gal(c_gal)%Mcoold
         MZc0 = gal(c_gal)%MZcd;   MZd0 = gal(c_gal)%MZd; Mtd0 = gal(c_gal)%Mtd
         ! --- AGN related
         Mbh0 = gal(c_gal)%Mbh
         ! --- added by MARK (2013/Aug/30)
         Tmass0 = gal(c_gal)%Tmass; Tlum_d0 = gal(c_gal)%Tlum_d
         ! --- added by MARK (2014/Nov/02)
         Zmassd0 = gal(c_gal)%Zmassd
         ! --- added by MARK (2014/Nov/07)
         SFR0 = gal(c_gal)%SFR; mSFR0 = gal(c_gal)%mSFR


         taust = taustar(gal(c_gal)%Vdisk, gal(c_gal)%rdisk, gal(c_gal)%clps, Mc0, MZc0, b_or_q)
         beta  = CalBeta(gal(c_gal)%Vdisk, b_or_q)
         
         IF(Mc0 > 0.d0 .and. MZc0 / Mc0 > 1.d0) print '(A)', '5'

         call star_formation
         call RenewalGalaxyQuiescent(c_gal)

         IF(param%dyn_resp_disk .or. param%dyn_resp_halo) THEN
            zi    = Cal_zi(Mini / gal(c_gal)%Mhalo, &
                           gal(c_gal)%Vc / gal(c_gal)%Vdisk)
            tidal = 1.d0
         ENDIF
         gal(c_gal)%diskmas = DiskMass(c_gal)
         IF(gal(c_gal)%diskmas > 0.d0) THEN
            IF(param%dyn_resp_disk) call DynamicalResponseDisk(c_gal)
            IF(param%dyn_resp_halo) call DynamicalResponseHalo(c_gal)
         ENDIF

1001     continue
         ! +-------------------------------------------------+
         ! ===============================================================


2001     continue
!!$     *** WRITING DATA OF GALAXIES WITHIN THIS HALO ***

         num_g_me = 0; i = num_total
         DO igal = 1, tmp%num_g ! loop for all galaxies in the halo w/ ID of "me"
            temp = GalMass(igal) + gal(igal)%Mhot + gal(igal)%Mhotorg
            IF(temp > 0.d0) THEN ! if the galaxy "igal" is not empty, its
               !  information "gal(igal)" is transfered into
               !  next timestep via "gal_next(igal)"
               i = i + 1
               gal_next(i) = gal(igal)
               IF(param%LAE) allgal_next(i) = allgal(igal)

               ! --- save the new ID of central in the halo "me" as "c_gal" temporally
               !      (added by MARK on 2013/Sep/25)
               IF(gal(igal)%flag_c == 1) c_gal = i

               num_g_me = num_g_me + 1
            ENDIF
            ! --- for ID trace (added by Shirakata 2018/Feb/08)
            IF(param%run_type /= 3 .and. param%traceIDs == 1 &
               .and. gal(igal)%Mbh > 0.d0 .and. gal(igal)%flag_c == 1) THEN
               DO j = 0, num_tar
                  IF(iforest == targ(j)%iforest) THEN
                     IF(gal(igal)%hfinal == targ(j)%fdes &
                        .and. gal(igal)%hstart == targ(j)%hstart) THEN
                        write(ionum+NFILE+nnode+nnode+inode+1, '(F6.3, I8, 3I16, 3I2, 15G13.5)')&
                           zp1, iforest, gal(igal)%mpi, &
                           gal(igal)%hstart, gal(igal)%hfinal, &
                           gal(igal)%flag_merger, gal(igal)%flag_di, gal(igal)%flag_ccut, &
                           gal(igal)%BT, gal(igal)%Mcoold, gal(igal)%Mstard, &
                           gal(igal)%Mstarb, gal(igal)%Mbh, &
                           gal(igal)%mem(1), gal(igal)%mem(2), &
                           gal(igal)%mem(3), gal(igal)%mem(4), &
                           gal(igal)%mem(5), gal(igal)%mem(6), &
                           gal(igal)%mem(7), gal(igal)%mem(8), &
                           gal(igal)%mem(9), gal(igal)%mem(10)

                           gal(igal)%mem(:) = 0.d0
                        exit
                     ENDIF
                  ENDIF
               ENDDO
               gal(igal)%mem(:) = 0.d0
            ENDIF
         ENDDO
         ! --- save the new ID of central in the halo "me" as "gal_next(igal)%id_cgal"
         !      (added by MARK on 2013/Sep/25)
         DO igal = num_total + 1, num_total + num_g_me
            gal_next(igal)%id_cgal = c_gal
         ENDDO
         mrgt%num_g(me) = num_g_me ! # of galaxies w/ non-zero baryon component
         !   in the halo "me"
         num_total = i ! cumulative # of galaxies w/ non-zero baryon in all halos
         !   from the halo "ihalo" = 1 at this timestep "itime"

         ! CSFR and SMF, added by Makiya
         IF(param%run_type == 1 &
            .or. (param%run_type == 2 .and. nz == 1) &
            .or. (param%run_type == 3 .and. nz == 1)) THEN
            DO igal = 1, num_g_me
               CSFR(2,itime) = CSFR(2,itime) + gal(igal)%SFR

               x = gal(igal)%Mstarb + gal(igal)%Mstard
               IF(x > 0.d0) THEN
                  SMF_base =  param%log10munit + 2.d0 * log10(param%h)
                  SMF_bin  = int(anint((log10(x) + SMF_base) / stepMF))
                  IF(SMF_bin > 0) &
                  SMF_z(itime, SMF_bin) = SMF_z(itime, SMF_bin) + 1.d0
               ENDIF
            ENDDO
         ENDIF
      ENDDO ID_ihalo
!!$  ===================================================================
!!$    ***   Loop for All Halos at the Time Snapshot "itime" Ends  ***
!!$  ===================================================================


      mrgp%num_tot_gal(itime) = num_total

      IF(flag_out == 1) THEN
         endhalo = num_total; end_step = itime
         IF(param%run_type /= 3) &
            print '(2(A, I8))', '# endhalo = ', endhalo, ', end_step = ', end_step

         ! --- added by MARK (2014/Nov/08)
         temp = param%munit / param%th_yr ! conv. factor of SFR to [Msun/yr]
         gal_next(1:num_total)%SFR  = gal_next(1:num_total)%SFR  * temp
         gal_next(1:num_total)%mSFR = gal_next(1:num_total)%mSFR * temp

         num_g_tot = 0; tmp%mhalo_sum = 0.d0 ! initialize
         iprog_max = mrgp%st_halo(itime) + mrgp%num_now(itime) - 1
         ! the max. ID of halo among all halos at the time "itime"
         DO iprog = mrgp%st_halo(itime), iprog_max
            num_g_tot     = num_g_tot     + mrgt%num_g(iprog) ! total # of gals
            tmp%mhalo_sum = tmp%mhalo_sum + mrgt%mhalo(iprog) ! total halo mass
         ENDDO
         IF(param%run_type /= 3) THEN
            print '(A, I8, 2(A, G10.3), A, G10.3, A)', &
               '# num_g_tot = ', num_g_tot, ', Mhalo_sum =', &
               tmp%mhalo_sum * param%munit, '[Msun]: ', &
               2.7755 * 7.**3 * 0.3 / 0.7, ', num_g_tot/mrgt%num_tot = ', &
               100.d0*dble(num_g_tot)/dble(mrgt%num_tot), '[%]'
         ENDIF

         deallocate(mgt,       stat=ier); call CheckIerr(ier, trim(cerr_d)//' mgt')
         deallocate(mhalo_sum, stat=ier); call CheckIerr(ier, trim(cerr_d)//' mhalo_sum')
         deallocate(tmp%lum,   stat=ier); call CheckIerr(ier, trim(cerr_d)//' tmp%lum')
         deallocate(tmp%lumq,  stat=ier); call CheckIerr(ier, trim(cerr_d)//' tmp%lumq')
         deallocate(tmp%agn,   stat=ier); call CheckIerr(ier, trim(cerr_d)//' tmp%agn')
         deallocate(merger%lum,stat=ier); call CheckIerr(ier, trim(cerr_d)//' merger%lum')
         deallocate(tmp%mem,   stat=ier); call CheckIerr(ier, trim(cerr_d)//' tmp%mem')

         RETURN
      ENDIF

      IF(num_total > 0 .and. flag_out == 0) THEN
         gal_prev(1:num_total) = gal_next(1:num_total)
         IF(param%LAE) allgal_prev(1:num_total) = allgal_next(1:num_total)
      ENDIF

      ! Okamoto+08 UV feedback mode
      IF(param%UVfb == 2 .and. flagO08_output == 1 .and. iO08_p+iO08_n0 > 0) THEN
         param%iO08_p  = param%iO08_p  + iO08_p
         param%iO08_n0 = param%iO08_n0 + iO08_n0
         param%iO08_n1 = param%iO08_n1 + iO08_n1
         param%iO08_n2 = param%iO08_n2 + iO08_n2
         param%iO08_n3 = param%iO08_n3 + iO08_n3
         param%iO08_n4 = param%iO08_n4 + iO08_n4

         temp = 100.d0 / dble(iO08_p + iO08_n0)
         print '(A, G10.3)', '   % of halo w/ Macc >= 0: ',   dble(iO08_p)  * temp
         print '(A, G10.3)', '   % of halo w/ Macc  < 0: ',   dble(iO08_n0) * temp
         print '(A, G10.3)', '     --- Macc+tmp%Mhotrog '//&
            '< 0                   : ', dble(iO08_n1) * temp
         print '(A, G10.3)', '     --- Macc+tmp%Mhotorg+'//&
            'tmp%Mhot < 0          : ', dble(iO08_n2) * temp
         print '(A, G10.3)', '     --- Macc+tmp%Mhotorg+'//&
            'tmp%Mhot+Mcool_all < 0: ', dble(iO08_n3) * temp
         print '(A, G10.3)', '     --- Macc+Mbar_all < 0'//&
            '                      : ', dble(iO08_n4) * temp
      ENDIF
   ENDDO ID_itime
!!$  ================================================
!!$    ***   Main Loop for Time Snapshot Ends   ***
!!$  ================================================
!!$ ========================================================================
 CONTAINS
!!$ ========================================================================
   SUBROUTINE SearchMostMassiveProgenitorHalo
     ! --- search for the most and 2nd most massive progenitor halo
     !     + substitute their consecutive IDs into mrgt%c_halo(me) and mrgt%c2_halo(me)
     ! --- mrgt%c_halo(me) corresponds to "childnum" in Mitaka model
     INTEGER :: itime, ihalo, iprog, iprog_max
     DOUBLE PRECISION :: m_big, m_big2

     DO itime = 1, mrgp%num_step ! loop for time snapshot
        DO ihalo = 1, mrgp%num_now(itime) ! loop for all halos at "itime"
           me = mrgp%st_halo(itime) + ihalo - 1 ! (ihalo,itime) --> me
           ! corresponding consecutive number for the halo "ihalo"
           !  at the time of "itime" in all snapshots
           mrgt%c_halo(me) = -1; mrgt%c2_halo(me) = -1; m_big = 0.d0; m_big2 = 0.d0

           iprog_max = mrgt%f_prg(me) + mrgt%numb(me) - 1
           ! the maximum ID of the "me"'s progenitor halos
           DO iprog = mrgt%f_prg(me), iprog_max ! loop for all prog. of "me"
              IF(mrgt%mhalo(iprog) > m_big2) THEN
                 IF(mrgt%mhalo(iprog) > m_big) THEN
                    m_big2 = m_big;             mrgt%c2_halo(me) = mrgt%c_halo(me)
                    m_big  = mrgt%mhalo(iprog); mrgt%c_halo(me)  = iprog
                 ELSE
                    m_big2 = mrgt%mhalo(iprog); mrgt%c2_halo(me) = iprog
                 ENDIF
              ENDIF
           ENDDO
        ENDDO
     ENDDO
   END SUBROUTINE SearchMostMassiveProgenitorHalo
!!$ ========================================================================
   SUBROUTINE InitializeLocalValiablesPerHalo
      ! --- These valiables are initialized 1 time per 1 halo
     flagburst = 0; flagdi = 0; flagmerger = 0
     c_gal = 0

     mhalo_sum(me) = 0.d0; Vsat2d = 0.d0; Vsat2b = 0.d0; dMhot  = 0.d0
     Mc0    = 0.d0; dMch  = 0.d0; dMchn = 0.d0; dMchp = 0.d0 
     dMcool = 0.d0; dMZhc = 0.d0
     xcool  = 0.d0 
     
     Mvrdyn_star = 0.d0; Mvrdyn_gas = 0.d0; Mvrdyn_gasZ = 0.d0
     Mme = mrgt%mhalo(me) ! halo mass of the halo "me"

     IF(param%UVfb == 2) THEN ! Okamoto+08 UV feedback mode
        Mbar_all  = 0.d0 ! total baryonic mass (incl. hot gas) in the halo
        Mgal_all  = 0.d0 ! total stellar + cold gas + BH mass in the halo
        Mcool_all = 0.d0 ! total cold mass in the halo
     ENDIF

     ! --- structure "tmp"
     tmp%num_g = 0 ! total # of galaxies in this halo
     tmp%Mhot  = 0.d0; tmp%MZh = 0.d0
     tmp%Morg  = 0.d0; tmp%Mhotorg = 0.d0; tmp%MZhorg = 0.d0

     ! --- structure "merger"
     merger%Msb   = 0.d0; merger%Msd  = 0.d0
     merger%Mcb   = 0.d0; merger%Mcd  = 0.d0
     merger%MZcb  = 0.d0; merger%MZcd = 0.d0
     merger%MZb   = 0.d0; merger%MZd  = 0.d0
     merger%Mtb   = 0.d0; merger%Mtd  = 0.d0
     merger%lum(:) = 0.d0
     merger%Mbh    = 0.d0
     merger%Tmassb = 0.d0; merger%Tmassd = 0.d0
     merger%Tlumb  = 0.d0; merger%Tlumd  = 0.d0
     merger%Zmassb = 0.d0; merger%Zmassd = 0.d0
   END SUBROUTINE InitializeLocalValiablesPerHalo
!!$ ========================================================================
   SUBROUTINE InitializeLocalValiablesPerGal
      ! --- These valiables are initialized 1 time per 1 galaxy.
     Vgrav   = 0.d0; Vmax  = 0.d0; r_di      = 0.d0; 
     f_dis   = 0.d0; f_dig = 0.d0; f_dis_rem = 1.d0; f_dig_rem = 1.d0
     fgas_DI = 0.d0

     fgas_sat  = 0.d0; fdisk_sat   = 0.d0; metallicity = 0.d0
     fvr       = 0.d0; fvr_rem     = 0.d0 
     fsat_vr   = 0.d0; fsat_vr_rem = 0.d0
     fgas_disk = 0.d0; fgas        = 0.d0 

     Mvr_star  = 0.d0; Mvr_gas     = 0.d0; Mvr_gasZ  = 0.d0
     Msat_star = 0.d0; Msat_gas   = 0.d0; Msat_gasZ = 0.d0
     massb = 0.d0; xy   = 0.d0
     M1    = 0.d0; M2   = 0.d0
     M1d   = 0.d0; M2d  = 0.d0
     M2M1  = 0.d0; RgRd = 0.d0
   END SUBROUTINE InitializeLocalValiablesPerGal
!!$ ========================================================================

   ! --- initialize the quantities related to the newly formed galaxy
   !       in the halo "me" at the time "itime"
   SUBROUTINE InitializeCentral
!!$    DOUBLE PRECISION :: dens, CircVel ! functions

     ! --- set and initialize temporary quantities
     tmp%num_g = 1 ! 0 --> 1 because of 1 newly formed gal in the halo "me"
     tmp%Mhot  = Mme * param%bar_rat
     IF(param%UVfb == 2) & ! Okamoto+08 UV feedback mode
          tmp%Mhot = Mme * fb_Okamoto08(Mme, zp1pc-1.d0, Mc_O08)
     tmp%MZh = 0.d0; tmp%Mhotorg = 0.d0; tmp%MZhorg = 0.d0

     ! --- set and initialize the quantities for gal(c_gal)
     c_gal  = 1 ! c_gal = igal = 1 because this galaxy is only 1 in the halo "me"
     !   at the time "itime"
     gal(c_gal)%flag_c = 1 ! flag for central galaxy is on because this galaxy is
     !  only one galaxy in the halo "me"
     gal(c_gal)%IDhost = me; gal(c_gal)%IDprog = -1
     gal(c_gal)%hstart = me ! added by Shirakata (2018/Feb/14)
     gal(c_gal)%hori   = mrgt%hori(me)
     gal(c_gal)%mpi    = mrgt%mpi(me) ! mpi of "me" is substituted into "gal(c_gal)%mpi"
     gal(c_gal)%Mhalo  = Mme; gal(c_gal)%clps = zp1; gal(c_gal)%Mreheat = 0.d0

     gal(c_gal)%flag_burst  = 0
     gal(c_gal)%flag_di     = 0
     gal(c_gal)%flag_merger = 0
     gal(c_gal)%flag_seed   = 0
     gal(c_gal)%flag_dead   = 0
     gal(c_gal)%flag_cfalse = 0

     gal(c_gal)%z_col   = 0.d0; gal(c_gal)%z_form1 = 0.d0
     gal(c_gal)%Mstarb  = 0.d0; gal(c_gal)%MZb = 0.d0; gal(c_gal)%Vbulge = 0.d0
     gal(c_gal)%Mstard  = 0.d0; gal(c_gal)%MZd = 0.d0; gal(c_gal)%Vdisk  = 0.d0
     gal(c_gal)%Vmax    = 0.d0
     gal(c_gal)%Mtb     = 0.d0; gal(c_gal)%rbulge = 0.d0
     gal(c_gal)%Mtd     = 0.d0; gal(c_gal)%rdisk  = 0.d0
     gal(c_gal)%Mcoold  = 0.d0; gal(c_gal)%MZcd   = 0.d0; gal(c_gal)%diskmas = 0.d0
     gal(c_gal)%Mcoolb  = 0.d0; gal(c_gal)%MZcb   = 0.d0
     gal(c_gal)%MZc_rem = 0.d0 ! added by MARK (2014/Sep/17)

     gal(c_gal)%lumg(:) = 0.d0
     ! -- added by MARK (2013/Aug/30)
     gal(c_gal)%Tmass = 0.d0; gal(c_gal)%Tlum_b = 0.d0; gal(c_gal)%Tlum_d = 0.d0

     ! --- added by MARK (2014/Nov/02)
     gal(c_gal)%Zmassb = 0.d0; gal(c_gal)%Zmassd = 0.d0

     ! --- added by MARK (2014/Nov/07)
     gal(c_gal)%SFR = 0.d0; gal(c_gal)%mSFR = 0.d0

     ! --- added by Shirakata (2016/03/24)
     gal(c_gal)%BT = 0.d0

     ! --- added by Shirakata (2017/02/27)
     gal(c_gal)%taumrg = 0.d0

     ! --- added by Shirakata (2017/08/24)
     gal(c_gal)%M2M1_av = 0.d0; gal(c_gal)%M2M1_max = 0.d0

     ! --- added by MARK (2017/Mar/16)
     gal(c_gal)%beta = 0.d0; gal(c_gal)%taust = 0.d0
     gal(c_gal)%dMstar_burst = 0.d0

     ! --- AGN related
     gal(c_gal)%Mbh     = 0.d0
     gal(c_gal)%lumq(:) = 0.d0
     gal(c_gal)%agn(:)  = 0.d0
     gal(c_gal)%Lpeak   = 0.d0

     IF(param%traceIDs == 1) THEN
        gal(c_gal)%mem(:) = 0.d0
     ENDIF
     ! --- added by Shirakata (2016/10/14)
     IF(param%LAE) call EraseGalaxyForLAE(c_gal)
     IF(param%SFH) call EraseGalaxyForSFH(c_gal)
   END SUBROUTINE InitializeCentral
!!$ ========================================================================
   ! --- return the total galaxy mass of the galaxy "id"
   DOUBLE PRECISION FUNCTION GalMass(id)
     INTEGER, INTENT(IN) :: id

     GalMass = gal(id)%Mstard + gal(id)%Mstarb + gal(id)%Mbh &
               + gal(id)%Mcoold + gal(id)%Mcoolb
   END FUNCTION GalMass
!!$ ========================================================================
   ! --- return the total baryon mass of the galaxy "id"
   DOUBLE PRECISION FUNCTION BaryonMass(id)
     INTEGER, INTENT(IN) :: id

     BaryonMass = gal(id)%Mstard + gal(id)%Mstarb + gal(id)%Mbh &
                  + gal(id)%Mcoold + gal(id)%Mcoolb &
                  + gal(id)%Mhot + gal(id)%Mhotorg
   END FUNCTION BaryonMass
!!$ ========================================================================
   ! --- return the bulge + SMBH mass of the galaxy "id"
   DOUBLE PRECISION FUNCTION BulgeMass(id)
     INTEGER, INTENT(IN) :: id
     BulgeMass = gal(id)%Mstarb + gal(id)%Mbh + gal(id)%Mcoolb
   END FUNCTION BulgeMass
!!$ ========================================================================
   ! --- return the disk mass (star + gas) of the galaxy "id"
   DOUBLE PRECISION FUNCTION DiskMass(id)
     INTEGER, INTENT(IN) :: id

     DiskMass = gal(id)%Mstard + gal(id)%Mcoold
   END FUNCTION DiskMass
!!$ ========================================================================
   ! --- return the stellar mass (bulge + disk) of the galaxy "id"
   DOUBLE PRECISION FUNCTION StellarMass(id)
     INTEGER, INTENT(IN) :: id

     StellarMass = gal(id)%Mstarb + gal(id)%Mstard
   END FUNCTION StellarMass
!!$ ========================================================================
   ! --- return the cold gas mass (bulge + disk) of the galaxy "id"
   DOUBLE PRECISION FUNCTION ColdMass(id)
     INTEGER, INTENT(IN) :: id

     ColdMass = gal(id)%Mcoolb + gal(id)%Mcoold
   END FUNCTION ColdMass
!!$ ========================================================================
   ! --- return the cold gas mass (bulge + disk) of the galaxy "id"
   DOUBLE PRECISION FUNCTION HotMass(id)
     INTEGER, INTENT(IN) :: id

     HotMass = gal(id)%Mhot + gal(id)%Mhotorg
   END FUNCTION HotMass
!!$ ========================================================================
   ! --- return the mass of heavy elements (bulge + disk) of the galaxy "id"
   DOUBLE PRECISION FUNCTION GalMZ(id)
     INTEGER, INTENT(IN) :: id

     GalMZ = gal(id)%MZb + gal(id)%MZd
   END FUNCTION GalMZ
!!$ ========================================================================
   ! --- return the mass of heavy elements (bulge + disk) in the gas phase
   ! --- of the galaxy "id"
   DOUBLE PRECISION FUNCTION ColdMZ(id)
     INTEGER, INTENT(IN) :: id

     ColdMZ = gal(id)%MZcb + gal(id)%MZcd
   END FUNCTION ColdMZ
!!$ ========================================================================
   DOUBLE PRECISION FUNCTION GalMt(id)
     INTEGER, INTENT(IN) :: id

     GalMt = gal(id)%Mtb + gal(id)%Mtd
   END FUNCTION GalMt
!!$ ========================================================================
   DOUBLE PRECISION FUNCTION KinEnergyBulge(id)
     INTEGER, INTENT(IN) :: id
     DOUBLE PRECISION :: Square ! function

     KinEnergyBulge = (gal(id)%Mstarb + gal(id)%Mbh + gal(id)%Mcoolb) &
                 * Square(gal(id)%Vbulge) 
   END FUNCTION KinEnergyBulge
!!$ ========================================================================
   DOUBLE PRECISION FUNCTION KinEnergyDisk(id)
     INTEGER, INTENT(IN) :: id
     DOUBLE PRECISION :: Square ! function

     KinEnergyDisk = (gal(id)%Mstard + gal(id)%Mcoold) * Square(gal(id)%Vdisk)
   END FUNCTION KinEnergyDisk
!!$ ========================================================================
   ! --- clean up all of the quantities of the galaxy "id"
   SUBROUTINE EraseGalaxy(id)
     INTEGER, INTENT(IN) :: id

     gal(id)%Mstarb = 0.d0; gal(id)%Mstard = 0.d0
     gal(id)%Mcoold = 0.d0; gal(id)%Mcoolb = 0.d0
     gal(id)%MZcd   = 0.d0; gal(id)%MZcb   = 0.d0
     gal(id)%MZb    = 0.d0; gal(id)%MZd    = 0.d0
     gal(id)%Mtb    = 0.d0; gal(id)%Mtd    = 0.d0
     gal(id)%MZc_rem = 0.d0 ! added by MARK (2014/Sep/17)

     gal(id)%lumg(:) = 0.d0

     ! --- added by MARK (2013/Aug/30)
     gal(id)%Tmass = 0.d0; gal(id)%Tlum_b = 0.d0; gal(id)%Tlum_d = 0.d0

     ! --- added by MARK (2014/Nov/02)
     gal(id)%Zmassb = 0.d0; gal(id)%Zmassd = 0.d0

     ! --- added by MARK (2014/Nov/07)
     gal(id)%SFR = 0.d0; gal(id)%mSFR = 0.d0

     ! --- added by Shirakata (2016/03/24)
     gal(id)%BT = 0.d0

     ! --- added by Shirakata (2017/02/27)
     gal(id)%taumrg = 0.d0
     gal(id)%z_col  = 0.d0; gal(id)%z_form1 = 0.d0

     ! --- added by MARK (2017/Mar/16)
     gal(id)%beta = 0.d0; gal(id)%taust = 0.d0
     gal(id)%dMstar_burst = 0.d0

     ! --- AGN related
     gal(id)%lumq(:) = 0.d0
     gal(id)%Mbh = 0.d0; gal(id)%agn(:) = 0.d0
     gal(id)%Lpeak = 0.d0

     gal(id)%flag_merger = 0; gal(id)%flag_di   = 0; gal(id)%flag_burst = 0
     gal(id)%flag_seed   = 0; gal(id)%flag_ccut = 0; gal(id)%flag_dead  = 0
     gal(id)%flag_cfalse = 0; gal(id)%n_merge = 0
     gal(id)%M2M1_av = 0.d0; gal(id)%M2M1_max = 0.d0

     IF(param%traceIDs == 1) gal(id)%mem(:) = 0.d0
     IF(param%LAE) call EraseGalaxyForLAE(id)
     IF(param%SFH) call EraseGalaxyForSFH(id)
   END SUBROUTINE EraseGalaxy
!!$ ========================================================================
   DOUBLE PRECISION FUNCTION Cal_rfri(mfmi, D)
     DOUBLE PRECISION, INTENT(IN) :: mfmi, D
     DOUBLE PRECISION :: tmp

     tmp = 0.5d0 * D
     Cal_rfri = (1.d0 + tmp) / (mfmi + tmp)
   END FUNCTION Cal_rfri
!!$ ========================================================================
   DOUBLE PRECISION FUNCTION Cal_rfri_disk(mfmi, D)
     DOUBLE PRECISION, INTENT(IN) :: mfmi, D
     DOUBLE PRECISION :: tmp, Square

     tmp = 1.d0 + mfmi / D
     Cal_rfri_disk = (tmp + sqrt(Square(tmp) - 4.d0 / D)) / 2.d0
   END FUNCTION Cal_rfri_disk
!!$ ========================================================================
   DOUBLE PRECISION FUNCTION Cal_zi(Mratio, Vratio)
     DOUBLE PRECISION, INTENT(IN) :: Mratio, Vratio
     DOUBLE PRECISION :: Square

     Cal_zi = 0.5d0 * Mratio * Square(Vratio)
   END FUNCTION Cal_zi
!!$ ========================================================================
   DOUBLE PRECISION FUNCTION Cal_vfvi(x, D, zi, zf)
     DOUBLE PRECISION, INTENT(IN) :: x, D, zi, zf
     DOUBLE PRECISION :: fxi, fxf, tmp

     fxi  = log(1.d0 + zi) / zi + log(1.d0 + 1.d0 / zi)
     fxf  = log(1.d0 + zf) / zf + log(1.d0 + 1.d0 / zf)
     tmp = 0.5d0 * D
     Cal_vfvi = sqrt((x + tmp * fxf) / (1.d0 + tmp * fxi))
   END FUNCTION Cal_vfvi
!!$ ========================================================================
   DOUBLE PRECISION FUNCTION Cal_vfvi_disk(x, D, zi, zf, cst)
     DOUBLE PRECISION, INTENT(IN) :: x, D, zi, zf, cst
     DOUBLE PRECISION :: fxi, fxf, tmp

     fxi = f_resp_disk(zi, cst); fxf = f_resp_disk(zf, cst)

     tmp = 1.d0/(0.25d0 * D)
     Cal_vfvi_disk = sqrt((x + tmp * fxf) / (1.d0 + tmp * fxi))
   END FUNCTION Cal_vfvi_disk
!!$ ========================================================================
   DOUBLE PRECISION FUNCTION f_resp_disk(z, c)
     DOUBLE PRECISION, INTENT(IN) :: z, c
     DOUBLE PRECISION :: cz  ! = c * z
     DOUBLE PRECISION :: cz1 ! = sqrt(cz**2 + 1.0)
     DOUBLE PRECISION :: cp1 ! = c + 1.0
     DOUBLE PRECISION :: Square ! function

     cp1 = c + 1.d0; cz  = c * z
     cz1 = sqrt(Square(cz) + 1.d0)
     f_resp_disk = c * (log(cz/2.d0) / cz &
                  + log((cz1 + 1.d0) * (cz1 / cz + 1.d0) * sqrt(cz1)) / (cz * cz1)) &
                  / (log(cp1) - c / cp1)
   END FUNCTION f_resp_disk
!!$ ========================================================================
   ! --- return the SN feedback intensity beta
   DOUBLE PRECISION FUNCTION CalBeta(velocity, b_or_q)
     DOUBLE PRECISION, INTENT(IN) :: velocity
     INTEGER, INTENT(IN) :: b_or_q ! 1:starburst, 2:quiescent

     CalBeta = (param%Vhot(b_or_q) / velocity)**param%alphot(b_or_q)
     beta    = CalBeta
   END FUNCTION CalBeta
!!$ ========================================================================
   ! --- calculate dynamical responce for bulge
   !      and renewal bulge velocity dispersion of the galaxy "id"
   SUBROUTINE DynamicalResponseBulge(id)
     INTEGER, INTENT(IN) :: id
     DOUBLE PRECISION :: D, mfmi, rfri, zf, vfvi, temp
     DOUBLE PRECISION :: Square ! function

     mfmi = BulgeMass(id) / Mini
     D    = 2.d0 * Square(gal(id)%Vcent / gal(id)%Vbulge)
     rfri = Cal_rfri(mfmi, D)
     gal(id)%rbulge = gal(id)%rbulge * rfri

     temp = zi / tidal
     zf   = rfri * temp
     vfvi = Cal_vfvi(mfmi/rfri, D, temp, zf)

     gal(id)%Vbulge = gal(id)%Vbulge * vfvi
   END SUBROUTINE DynamicalResponseBulge
!!$ ========================================================================
   ! --- renewal circular velocity of the host halo of the galaxy "id"
   !       via dynamical response
   SUBROUTINE DynamicalResponseHalo(id)
     INTEGER, INTENT(IN) :: id

     gal(id)%Vcent = gal(id)%Vcent &
                     * RespDM(GalMass(id), Mini, gal(id)%Mhalo * zi * tidal)
   END SUBROUTINE DynamicalResponseHalo
!!$ ========================================================================
   DOUBLE PRECISION FUNCTION RespDM(Mfinal, Minitial, x)
     DOUBLE PRECISION, INTENT(IN) :: Mfinal, Minitial, x
     DOUBLE PRECISION :: tmp

     tmp    = param%fdm * x
     RespDM = (Mfinal + tmp) / (Minitial + tmp)
   END FUNCTION RespDM
!!$ ========================================================================
   ! --- calculate dynamical responce for disk
   !      and renewal disk rotation velocity of the galaxy "id"
   SUBROUTINE DynamicalResponseDisk(id)
     INTEGER, INTENT(IN) :: id
     DOUBLE PRECISION :: D, mfmi, rfri, zf, vfvi, temp
     DOUBLE PRECISION :: cst ! concentration of dark halo
     DOUBLE PRECISION :: czi ! cst * zi
     DOUBLE PRECISION :: cstp1 ! = cst + 1.0
     DOUBLE PRECISION :: log_czih ! = log(czi/2)
     DOUBLE PRECISION :: D1, D2, D3, D4, D5
     DOUBLE PRECISION :: Square, Cube ! functions
     DOUBLE PRECISION :: CalCST ! concentration param. of DM halo

     zi = 0.1d0
     IF(gal(id)%rdisk > 0.d0) &
          zi = (0.06d0 * gal(id)%rdisk / sqrt(2.d0)) &
               / (0.5d0 * 6.67d0 * gal(id)%Mhalo * 2.d+8 &
                  / (Square(gal(id)%Vc) * 3.08d0))

     ! cst = 10.d0
     cst   = CalCST(gal(id)%clps, gal(id)%Mhalo)
     cstp1 = cst + 1.d0

     mfmi = GalMass(id) / Mini; czi = cst * zi; log_czih = log(czi * 0.5d0)

     D1 = cst * gal(id)%Mhalo / Mini
     D2 = log(cstp1) - cst / cstp1
     D3 = czi * zi * (3.d0 + 2.d0 * log_czih)
     D4 = Square(czi) * zi * (-16.d0 / 3.d0)
     D5 = Cube(czi) * zi * (-33.d0 / 8.d0 - 4.5d0 * log_czih)
     D  = D1 * (D3 + D4 + D5) / D2

     rfri = Cal_rfri_disk(mfmi, D)
     gal(id)%rdisk = gal(id)%rdisk * rfri
     temp = zi / tidal
     zf   = rfri * temp
     vfvi = Cal_vfvi_disk(mfmi/rfri, Mini/(gal(id)%Mhalo*zi), temp, zf, cst)
     gal(id)%Vdisk = gal(id)%Vdisk * vfvi
   END SUBROUTINE DynamicalResponseDisk
!!$ ========================================================================
   ! --- calculate cooling process and renewal the hot and cold gas components
   !       of the galaxy "id"
   SUBROUTINE CoolingProcess
     ! --- AGN feedback related
     DOUBLE PRECISION :: dMbh_radio, dMcool_orig ! for Croton like AGN feedback
     DOUBLE PRECISION :: t_ff, Lcool, Tvir, Ledd, t_dyn, t_cool, t_coolp
     DOUBLE PRECISION :: rrs, rrsn, rrsp, xcoolp
     DOUBLE PRECISION :: cst0
     ! for Bower like AGN feedback
     DOUBLE PRECISION :: dMbh_tmp
     ! for hot gas accretion (added by Shirakata 2016/Nov/10)
     DOUBLE PRECISION :: Mr2h, MZr2h, frac_ret
     DOUBLE PRECISION :: Mhalo, Vc, density
     DOUBLE PRECISION :: Square, Cube ! functions
     DOUBLE PRECISION :: M_DM

     ! --- Current halo properties
     !   Mhalo   = Mme
     !   density = dens(zp1)
     !   Vc      = CircVel(Mhalo, density)
     !   cst0    = CalCST(Mhalo,zp1)

     ! --- Halo properties at formation time
     Mhalo   = MIN(gal(c_gal)%Morg, Mme)
     density = gal(c_gal)%density
     Vc      = gal(c_gal)%Vc
     cst0    = gal(c_gal)%cst

     ! --- Calculate Vmax for disk instability
     gal(c_gal)%Vmax = 0.465d0 * sqrt(cst0 / (log(1+cst0) - cst0/(1.d0+cst0))) * Vc

     Zh0 = 0.d0
     IF(flagdh == 1 .or. mrgt%numb(me) == 0) THEN
        IF(gal(c_gal)%Mhot > 0.) THEN
           Zh0 = min(1.d0, gal(c_gal)%MZh / gal(c_gal)%Mhot)
           gal(c_gal)%Mratio = gal(c_gal)%Mhot / Mhalo
           call cool1(gal(c_gal)%Telapse, gal(c_gal)%Mratio, Vc,  density, &
                      dMch, zp1, Mhalo, Tvir, rrs, xcool, cst0, t_cool)
           gal(c_gal)%Mhotorg = gal(c_gal)%Mhot
           gal(c_gal)%MZhorg  = gal(c_gal)%MZh
           gal(c_gal)%Mhot = 0.d0; gal(c_gal)%MZh = 0.d0; dMchn = dMch
        ELSE
           dMch = 0.d0; rrs = 0.d0; xcool = 0.d0
        ENDIF
     ELSE
        ! --- even if a halo does not major merge
        ! --- the hot gas accretion occur (Shirakata 2016/Nov/10)

        t_dyn = const%t_dyn0 / sqrt(density) ! [yr]
        !    const%t_dyn0 = sqrt(3.d0 * const%PI / (32.d0 * 6.674d-11 * 1.88d-26
        !                                                 * Square(param%h)))
        !                   / const%yr2sec defined in main_nugc.f90

        frac_ret = param%alp_ret

        Mr2h  = frac_ret * gal(c_gal)%Mhot * min(tlife_yr / t_dyn, 1.d0)
        MZr2h = frac_ret * gal(c_gal)%MZh  * min(tlife_yr / t_dyn, 1.d0)

        gal(c_gal)%Mhot    = gal(c_gal)%Mhot - Mr2h
        gal(c_gal)%MZh     = gal(c_gal)%MZh  - MZr2h
        gal(c_gal)%Mhotorg = gal(c_gal)%Mhotorg + Mr2h
        gal(c_gal)%MZhorg  = gal(c_gal)%MZhorg  + MZr2h
        gal(c_gal)%Mratio  = gal(c_gal)%Mhotorg / Mhalo

        IF(gal(c_gal)%Mhotorg > 0.) THEN 
             Zh0 = min(1.d0, gal(c_gal)%MZhorg / gal(c_gal)%Mhotorg)

           dMch = 0.d0
           call cool1(gal(c_gal)%Telapse, gal(c_gal)%Mratio, Vc, density, &
                      dMchn, zp1, Mhalo, Tvir, rrsn, xcool, cst0, t_cool)
           rrs = rrsn

           IF(dMchn * Mhalo <= Mme) THEN
              call cool1(gal(c_gal)%Telapse - tlife, gal(c_gal)%Mratio, Vc, density, &
                         dMchp, zp1, Mhalo, Tvir, rrsp, xcoolp, cst0, t_coolp)
              dMch = dMchn - dMchp
           ELSE
              dMch = 0.d0; xcool = 0.d0; rrs = 0.d0
           ENDIF
        ELSE
           dMch = 0.d0; xcool = 0.d0; rrs = 0.d0
        ENDIF
     ENDIF

     dMcool = gal(c_gal)%Mratio * Mhalo * dMch

     IF(gal(c_gal)%Mhotorg < dMcool) THEN
        dMcool = gal(c_gal)%Mhotorg; dMZhc = gal(c_gal)%MZhorg
     ENDIF

     !========= for AGN feedback ==========================================
     CoolingCutoff_key = 0; dMbh_tmp = 0.d0
     IF(dMcool > 0.d0) THEN
        IF(param%AGNFB_key == 1) THEN ! Vcut model
           IF(Vccool >= param%Vcut) THEN
              CoolingCutoff_key = 1
              dMcool = 0.d0; dMch = 0.d0
           ENDIF
        ELSE IF(param%AGNFB_key == 2) THEN ! Croton+06
           dMcool_orig = dMcool

           ! eq.(10) of Croton+06. kappa is 6e-6 Msun/yr in their model
           dMbh_radio = param%kappa_croton * (gal(c_gal)%Mbh*param%munit / 1.d+8) &
                        * (gal(c_gal)%Mhotorg / (gal(c_gal)%Mhalo * 0.1d0)) &
                        * Cube(Vccool / 200.d0) / param%munit
           dMcool = max(0.d0, dMcool - (2.d0 * param%eta_croton * dMbh_radio &
                                        * Square(const%c / (Vccool * 1.d+3))) * tlife_yr)
           dMch  = dMcool / dMcool_orig * dMch
           dMchn = dMchp + dMch

           dMbh_tmp = dMbh_radio * tlife_yr
           gal(c_gal)%Mbh = gal(c_gal)%Mbh + dMbh_tmp
        ELSE IF(param%AGNFB_key == 3) THEN ! Bower+06
!!$          Lcool = 1.5d0 * dMcool * param%munit / (param%mu * const%mp) &
!!$                        * const%kB * Tvir / tlife_sec
!!$          Lcool = Lcool * 1.d+41 ! [erg/s]
!!$          Lcool = param%Lcool0 * dMcool * Tvir / tlife_sec ! [erg/s]
           Lcool = param%Lcool0 * dMcool * Tvir &
                   / (gal(c_gal)%Telapse * param%th_yr * 3.15d+7) ! [erg/s]
         ! param%Lcool0 = 1.5d+41 * param%munit * const%kB
         !                / (param%mu * const%mp) defined in main_nugc.f90
!!$          Ledd = 1.26d0 * 1.d+38 * (gal(c_gal)%Mbh * param%munit)
           Ledd = param%Ledd0  * gal(c_gal)%Mbh
         ! Eddington lum. of central BH [erg/s]
         ! param%Ledd0 = 1.26d+38 * param%munit defined in main_nugc.f90

!!$         t_ff = 5.02d+4 * sqrt(HotMass(c_gal) / Mhalo) &
!!$                * sqrt((rrs*CalRdisk(Vc,density)/(cst0*param%h))**3.d0) &
!!$                / sqrt(dMcool)  ! [yr]

           M_DM = Mhalo * (1.d0 - param%bar_rat) &
                  * (log(1.d0+rrs) + 1.d0 / (1.d0+rrs) - 1.d0) &
                  / (log(1.d0+cst0) + 1.d0 / (1.d0+cst0) - 1.d0)
           t_ff = 5.02d+4 * sqrt((rrs * CalRdisk(Vc,density) &
                                 / (cst0 * param%h))**3.d0 / M_DM) ! [yr]

           IF(param%alp_bower * min(t_cool, gal(c_gal)%Telapse) * param%th_yr >= t_ff &
                .and. param%eps_bower * Ledd > Lcool) THEN

              CoolingCutoff_key = 1
!!$            dMbh_tmp = Lcool * tlife_sec * 1.d-7 / (Square(const%c) * 2.d-1 &
!!$                                                    * const%Msolar * param%munit)

              dMbh_tmp = const%dMbh0 * Lcool * tlife_sec
            ! const%dMbh0 = 1.d-7 / (Square(const%c) * 2.d-1 * const%Msolar
            !                        * param%munit) defined in main_nugc.f90
              gal(c_gal)%Mbh = gal(c_gal)%Mbh + dMbh_tmp
              dMcool = 0.d0; dMch = 0.d0; xcool = 0.d0
           ENDIF
        ENDIF
     ELSE  ! dMcool <= 0.d0
        dMcool = 0.d0; dMch = 0.d0; xcool = 0.d0
     ENDIF
!==========================
     dMZhc = dMch * gal(c_gal)%MZhorg
     gal(c_gal)%Mhotorg = max(0.d0, gal(c_gal)%Mhotorg - dMcool - dMbh_tmp)
     gal(c_gal)%Mcoold  = gal(c_gal)%Mcoold + dMcool
     gal(c_gal)%MZhorg = max(0.d0, gal(c_gal)%MZhorg - dMZhc)
     gal(c_gal)%MZcd   = gal(c_gal)%MZcd + dMZhc
   END SUBROUTINE CoolingProcess
!!$ ========================================================================
   ! --- renewal the quantities of the galaxy "id" via quiescent star-formation
   SUBROUTINE RenewalGalaxyQuiescent(id)
     INTEGER, INTENT(IN) :: id
     DOUBLE PRECISION :: ab, taueff, Zmass, u, eu, ts, us, eus
     INTEGER          :: nseed,clock
     INTEGER, ALLOCATABLE :: seed_random(:)
     DOUBLE PRECISION :: fseed, Mseed
     ! DOUBLE PRECISION :: temp, temp1, temp2 ! for check


     gal(id)%Mstard  = Ms0 + dMstar
     gal(id)%Mcoold  = Mc0 - dMcold
     gal(id)%MZcd    = gal(id)%Mcoold * Zc
     gal(c_gal)%Mhot = gal(c_gal)%Mhot + dMhot
     gal(c_gal)%MZh  = gal(c_gal)%MZh  + dMZh

     IF(gal(id)%Mbh < param%Mbhseed/param%munit) THEN
!!$      call random_seed(size=nseed)
!!$      allocate(seed_random(nseed))
!!$      call system_clock(count=clock)
!!$      seed_random = clock
!!$      call random_seed(put=seed_random)
!!$      call random_number(fseed)
!!$      fseed = 2.d0 * fseed
!!$      Mseed = fseed + log10(param%Mbhseed)
!!$      Mseed = 10 ** Mseed   ! Msun

!!$      dMbh = min(Mseed/param%munit, gal(c_gal)%Mhot)
        dMbh = min(param%Mbhseed/param%munit, gal(c_gal)%Mhot)
        gal(c_gal)%Mhot = gal(c_gal)%Mhot - dMbh
        gal(id)%Mbh = gal(id)%Mbh + dMbh
        dMbh = 0.d0
        gal(id)%flag_seed = 1; gal(id)%z_col = zp1 - 1.d0 
     ENDIF
     fseed = 0.d0; Mseed = 0.d0
     Mcool0 = max((Mc0-dMbh_norm), 0.d0) * param%munit
     IF(param%LAE) call CalAllGalForLAE(b_or_q, id, Mcool0, zp1pc)

     x = 0.d0
     call lum(zp1pc, zp1, tmp%lum, mz, mt, x)

     gal(id)%lumg(param%nwp1:param%tnw) &
          = gal(id)%lumg(param%nwp1:param%tnw) &
            + tmp%lum(1:param%nwave) * Mcool0
     gal(id)%MZd = MZd0 + mz * Mcool0
     gal(id)%Mtd = Mtd0 + mt * Mcool0

     gal(id)%SFR = SFR0 + dMstar / (tlife * ssp%alp(b_or_q))

     ! --- added by Shirakata (2016/03/24)
     gal(id)%BT = 0.d0
     IF(gal(id)%lumg(param%iBband_r) + gal(id)%lumg(param%iBband_r+param%nwave) &
          > 0.d0) THEN
        gal(id)%BT = gal(id)%lumg(param%iBband_r) &
                     / (gal(id)%lumg(param%iBband_r) &
                        + gal(id)%lumg(param%iBband_r+param%nwave))
     ENDIF

     ! --- added by MARK (2014/Nov/05)
     IF(tlife < const%t0) THEN
        ! tlife is shorter than the timescale during which mSFR is evaluated
        gal(id)%mSFR = gal(id)%SFR
     ELSE
        !  in the case that tlife is longer than the timescale during which mSFR
        !  is evaluated const%t0, only the contribution of the stellar mass newly
        !  formed during the past const%t0 is added into gal(id)%mSFR
        ab = ssp%alp(b_or_q) + beta; taueff = taust / ab
!!$      IF(taueff < 1.d-3) taueff = 1.d-3 ! to avoid divergence of mSFR
!!$                                        !  (added by Makiya, 2015/Nov/16)
!!$      gal(id)%mSFR = mSFR0 &
!!$         + gal(id)%Mcoold * (exp(const%t0 / taueff) - 1.d0) / (ab * const%t0)
        ! --- modified by MARK (2017/Mar/15) to avoid divergence of mSFR without
        !      modification of taueff
        gal(id)%mSFR = mSFR0 &
             + Mc0 * (exp(-(tlife-const%t0) / taueff) - exp(-tlife / taueff)) &
               / (ab * const%t0)
!!$      ! --- for check
!!$      IF(taust / ab < 1.d-3) THEN
!!$         taueff = taust / ab
!!$         temp = param%munit / param%th_yr
!!$         temp1 = mSFR0 + gal(id)%Mcoold * (exp(const%t0 / taueff) - 1.d0) / (ab * const%t0)
!!$         temp2 = mSFR0 + Mc0 * (exp(-(tlife-const%t0) / taueff) - exp(-tlife / taueff)) &
!!$                             / (ab * const%t0)
!!$         print '(A, I10, 6(G13.5,1X))', '# taueff < 1d-3: ', id, taueff, &
!!$              temp1 * temp, gal(id)%mSFR * temp, gal(id)%mSFR / temp1 - 1.d0, &
!!$              temp2 * temp, temp1 / temp2 - 1.d0
!!$      ENDIF
     ENDIF

     ! --- added by MARK (2017/Mar/16)
     gal(id)%taust = taust

     ! --- added by MARK (2013/Aug/30)
     gal(id)%beta   = beta; ab = ssp%alp(b_or_q) + beta; taueff = taust / ab
     gal(id)%Tlum_d = Tlum_d0 + mt * Mcool0
     gal(id)%Tmass  = Tmass0 + ssp%alp(b_or_q)/ab * Mc0 * tlife &
                             + dMstar * (t_zsp1 - taueff - t_zp1pc)

     ! --- added by MARK (2014/Nov/02); see eq.(21) in Nagashima+05
     !     mass-weighted mean metallicity of the newly formed stars in this SF
     u = tlife / taueff; eu = exp(-u)
     IF(Zc0 + ssp%p(b_or_q) * tlife / taust > 1.d0) THEN
        ts = taust * (1.d0 - Zc0) / ssp%p(b_or_q) ! the time at which Zc(ts) = 1
        us = ts / taueff; eus = exp(-us)

        Zmass = ((Zc0 + ssp%p(b_or_q) / ab) * (1.e0 - eus) &
                - ssp%p(b_or_q) * ts * eus / taust + (eus - eu)) / (1.d0 - eu)
     ELSE
        Zmass = Zc0 + ssp%p(b_or_q) / ab * (1.d0 - eu - u * eu) / (1.d0 - eu)
     ENDIF
     gal(id)%Zmassd = Zmassd0 + Zmass * dMstar

     ! --- added by Shirakata (2017/Jul/13)
     IF(gal(id)%flag_burst ==0) gal(id)%MZc_rem = 0.d0
     
     IF(param%SFH) call CalSFHrelatedInCom(b_or_q, id)
   END SUBROUTINE RenewalGalaxyQuiescent
!!$ ========================================================================
   ! --- renewal the quantities of the galaxy "id" via starburst
   SUBROUTINE RenewalGalaxyBurst(id)
     INTEGER, INTENT(IN) :: id
     DOUBLE PRECISION :: ab
     DOUBLE PRECISION, PARAMETER :: EPS = 1.d-10
     DOUBLE PRECISION :: seedmass


     gal(id)%Mstarb  = Ms0 + dMstar
     gal(id)%Mcoolb  = Mc0 - dMcold ! dMcold = Mc0 in star_formation_burst
     gal(id)%MZcb    = (Mc0 - dMcold) * Zc
     gal(c_gal)%Mhot = gal(c_gal)%Mhot + dMhot
     gal(c_gal)%MZh  = gal(c_gal)%MZh  + dMZh

     Mcool0 = (Mc0 - dMbh) * param%munit
     IF(param%LAE) call CalAllGalForLAE(b_or_q, id, Mcool0, zp1pc)

     call lum_burst(zp1pc, zp1, tmp%lum, mz, mt, x)
     ! "x" is determined before this subroutine is called

     gal(id)%lumg(1:param%nwave) &
          = gal(id)%lumg(1:param%nwave) + tmp%lum(1:param%nwave) * Mcool0
     gal(id)%MZb = MZb0 + mz * Mcool0
     gal(id)%Mtb = Mtb0 + mt * Mcool0
     gal(id)%SFR = SFR0 + dMstar / (tlife * ssp%alp(b_or_q))
                   ! mean SFR during the timestep

     ! --- added by Shirakata (2016/03/24)
     gal(id)%BT = 0.d0
     IF(gal(id)%lumg(param%iBband_r) + gal(id)%lumg(param%iBband_r+param%nwave) &
          > 0.d0) THEN
        gal(id)%BT = gal(id)%lumg(param%iBband_r) &
                     / (gal(id)%lumg(param%iBband_r) &
                        + gal(id)%lumg(param%iBband_r+param%nwave))
     ENDIF

     ! --- added by MARK (2014/Nov/05)
     gal(id)%mSFR = mSFR0
     IF(tlife - x < const%t0) &
          ! gal(id)%mSFR = gal(id)%mSFR + dMstar / ((tlife - x) * ssp%alp(b_or_q))
          gal(id)%mSFR = gal(id)%mSFR + dMstar / (const%t0 * ssp%alp(b_or_q))
     ! starburst occurs at the look-back time of tlife - x smaller than the
     !  timescale during which mSFR is evaluated

     ! --- added by MARK (2017/Mar/16)
     gal(id)%taust = taust
     gal(id)%dMstar_burst = dMstar / ssp%alp(b_or_q)
                            ! stellar mass newly formed in this starburst
                            !  (returned mass fraction corrected)

     ! --- added by MARK (2013/Aug/30)
     gal(id)%beta   = beta; ab = ssp%alp(b_or_q) + beta
     gal(id)%Tlum_b = Tlum_b0 + mt * Mcool0
     IF(abs(zp1pc - param%zsp1) <= EPS .and. tlife > 0.d0) THEN
        gal(id)%Tmass = Tmass0
     ELSE
        gal(id)%Tmass = Tmass0 + dMstar * (t_zsp1 - t_zp1pc + tlife - x)
     ENDIF

!!$    ! --- added by MARK (2014/Nov/02)
!!$    ! --- all of the newly formed stars are assumed to have the identical
!!$    !      metallicity of Zc0
!!$    gal(id)%Zmassb = Zmassb0 + Zc0 * dMstar
     ! --- modified by MARK (2015/May/27) complained by Nagashima & Okamoto
     gal(id)%Zmassb = Zmassb0 + (Zc0  + ssp%p(b_or_q) / ab) * dMstar

     ! --- AGN related
     gal(id)%Mbh = Mbh0
     IF(gal(id)%Mbh <= param%Mbhseed/param%munit) THEN
        seedmass = min(param%Mbhseed/param%munit, gal(c_gal)%Mhot)
        gal(id)%Mbh     = gal(id)%Mbh + seedmass
        gal(c_gal)%Mhot = gal(c_gal)%Mhot - seedmass
        gal(id)%flag_seed = 1; gal(id)%z_col = zp1 - 1.d0
     ENDIF
     call lumagn(itime, id, zp1, zp1pc, x, tmp%lumq, tmp%agn, tmp%mem)
     gal(id)%lumq(1:9) = tmp%lumq(1:9)
     gal(id)%agn(1:7)  = tmp%agn(1:7)
     gal(id)%agn(8:9)  = 0.d0
     gal(id)%agn(10)   = tmp%agn(10)
     IF(gal(id)%flag_merger == 1) THEN
        gal(id)%agn(8)   = gal(id)%M2M1_av
        gal(id)%agn(9)  = gal(id)%M2M1_max
     ENDIF

     IF(param%traceIDs == 1) THEN
        gal(id)%mem(:) = tmp%mem(:)
     ENDIF
     IF(param%SFH) call CalSFHrelatedInCom(b_or_q, id)
   END SUBROUTINE RenewalGalaxyBurst
!!$ ========================================================================
   SUBROUTINE EstimateMZcRem(id)
     INTEGER, INTENT(IN) :: id

     taudyn = gal(id)%rbulge / gal(id)%Vbulge / 10.d0
     ! taudyn = H0*t_dyn = 0.1*(r[kpc/h])/(V[km/s])
     tburst = param%alpha_burst * taudyn
     temp   = (tlife - x) / tburst
     gal(id)%MZc_rem = Mc0 * exp(-(ssp%alp(1) + beta) * temp) &
                           * min(1.d0, Zc0 + ssp%p(1) * temp) + gal(id)%MZcb
     IF(gal(id)%MZc_rem < 1.d-50) gal(id)%MZc_rem = 0.d0
   END SUBROUTINE EstimateMZcRem
!!$ ========================================================================
   SUBROUTINE ResetFlags(id)
     INTEGER, INTENT(IN) :: id

     gal(id)%flag_burst  = 0
     gal(id)%flag_di     = 0
     gal(id)%flag_merger = 0; gal(id)%flag_seed = 0
     gal(id)%flag_cfalse = 0; gal(id)%n_merge   = 0
     gal(id)%M2M1_av     = 0.d0; gal(id)%M2M1_max = 0.d0
     IF(gal(id)%flag_c == 0) gal(id)%flag_ccut = 0
     IF(param%traceIDs == 1) gal(id)%mem(:) = 0.d0
   END SUBROUTINE ResetFlags
!!$ ========================================================================
 END SUBROUTINE star
!!$=======================================================================
!!$     ------------------ SUBROUTINE cool --------------------------
SUBROUTINE cool1(tl, Mratio0, Vctmp, density0, dMch, zp1, Mhalo, Tvir, rrs1, &
                 xcoolvir, cst0, t_cool)
   use global_var; use CLCOOLrelated 
   implicit none

   DOUBLE PRECISION, INTENT(IN)    :: tl, Mratio0, Vctmp, density0
   DOUBLE PRECISION, INTENT(IN)    :: zp1, Mhalo, cst0
   DOUBLE PRECISION, INTENT(INOUT) :: dMch, rrs1, xcoolvir, t_cool
   DOUBLE PRECISION, INTENT(INOUT) :: Tvir
   DOUBLE PRECISION, PARAMETER :: EPS = 1.d-1
   INTEGER :: i_yy, i_yyp1, jj
   DOUBLE PRECISION :: dMchcool, yy, lami, lamT, lamis, lamTs
   DOUBLE PRECISION :: Ac, rrc2, Bc, rrc1, rffrc2
   DOUBLE PRECISION :: f1, f2, fmid, delx, cst
   DOUBLE PRECISION :: x1, x2, xmid
   DOUBLE PRECISION :: nh,ne
   DOUBLE PRECISION :: tlconst
   DOUBLE PRECISION :: tmp, Square, Cube
   DOUBLE PRECISION :: FFR

!!$         Ac = 1.d0 / 1.40253d-3 ! c = 10,  Ac = 1/A

   !--- concentration parameter of DM halo
   cst  = cst0

   !--- initialization of local variables
   rrc1   = 0.d0; rrs1 = 0.d0; xcoolvir = 0.d0
   t_cool = tl

   tmp = cst / 0.22d0
   tmp = 1.d0  / (tmp - atan(tmp))
   Ac  = (31.3d0 * Cube(cst)) * tmp
!!$         Bc   = 2.27761d-2
   Bc   = tmp
   Tvir = 60.5d0 * param%mu * Square(Vctmp)

   dMchcool = 0.d0
   rffrc2   = 1.d+10

   IF(param%CoolFN == 1) THEN ! SD93 cooling function
      IF(Tvir  > param%Tvir_max) Tvir = param%Tvir_max
      IF(Tvir <= param%Tvir_min) THEN
         dMch = 0.d0
         goto 1112
      ENDIF

      yy    = 20.d0 * log10(Tvir) - 79.d0
      i_yy  = int(yy); i_yyp1 = i_yy + 1
      lami  = lambda(i_yy); lamis = lambdas(i_yy)
      tmp   = yy - aint(yy)
      lamT  = (lambda(i_yyp1)  - lami)  * tmp + lami
      lamTs = (lambdas(i_yyp1) - lamis) * tmp + lamis
      lamT  = -(lamT - lamTs) * Zh0 / const%Zsun + lamT

      tlconst = 7.11d4
      ! Assuming fully ionized gas with solar metallicity
      ! n_H = (12/27) * ngas, n_e = 1.209 * n_H
      ! (See Sutherland & Dopita 1993 Table 6. 283p)
      ! Their definition: Cooling func. = cool. rate / (n_H*n_E)
   ELSE IF(param%CoolFN == 2) THEN ! Okamoto-cloudy
      nh    = 5.29d-6 * param%h * param%h * Ac * Mratio0 * density0 / param%mu
      ! --- Assuming rho_gas = f_gas * delta_c * rho_crit(z)
      ! --- and nh = (12/27) * n_gas (ngas = rho_hot / mu / m_p)
      ! --- Shirakata (2018/Aug/27)
      call CalcHeatCoolRate(Tvir, Zh0, nh, zp1, clcool_table_z, lamT)
      tlconst = 5.95d4
      ! Definition: Cooling func. = cool. rate / (n_H*n_H)
   ENDIF

!!$     dMchcool = 302. * (tlife * param%h * Mratio * density / Tvir)**(0.5) &
!!$         * 10.**(0.5 * (lamT + 23.))
!!$     rrc2 = 2.81d+5 * param%h * (10.**(lamT + 23.)) * density0 * Mratio0 &
!!$         * tl * Ac / Tvir - 1.
   rrc2 = MAX(tlconst * param%h * (10.d0**(lamT + 23.d0)) * density0 * Mratio0 &
          * tl * Ac / Tvir - 1.d0, 1e-5)
   
   IF(rrc2 > 0.d0) THEN
      ! --- Cooling time at r = R_vir [hubble time]
      ! --- If Telapse > t_cool (Rvir), the cooling time 
      !       which is used for Bower AGN feedback condition
      !       should set t_cool (R_vir) (Shirakata 2017/Feb/22)
      t_cool = Tvir * (1.d0 + cst / 0.22d0) * (10.d0**(-lamT - 23.d0)) &
               / (tlconst * param%h * density0 * Mratio0 * Ac)

      ! --- compare r_cool with r_ff (Shirakata 2016/Nov/10)
      ! --- Search the rffrc2 within cooling radius
      delx = EPS
      IF(FFR(0.22d0*sqrt(rrc2), tl, density0, cst, Vctmp, Mhalo) < 0.d0) &
           goto 1111
      ! rffrc2 > rrc2

      rffrc2 = rrc2 ! Initialize

      x1 = delx; x2 = x1  ! x1, x2, xmid = r/rs
      f1 = FFR(x1, tl, density0, cst, Vctmp, Mhalo)
      DO WHILE(x1 <= 0.22d0*sqrt(rrc2))
         f2 = FFR(x1 + delx, tl, density0, cst, Vctmp, Mhalo)
         IF(f1 * f2 < 0.0) THEN
            x2 = x1 + delx
            exit
         ENDIF
         f1 = f2; x1 = x1 + delx
      ENDDO

      x1 = x2 - delx; xmid = (x1 + x2) * 0.5d0
      f1 = FFR(x1, tl, density0, cst, Vctmp, Mhalo)
      f2 = FFR(x2, tl, density0, cst, Vctmp, Mhalo)
      jj = 0
      DO WHILE(dabs(x1 - x2) > EPS)
         jj = jj + 1
         xmid = (x1 + x2) * 0.5d0
         fmid = FFR(xmid, tl, density0, cst, Vctmp, Mhalo)
         IF(fmid * f1 >= 0.d0) THEN
            x1 = xmid; f1 = fmid
         ELSE IF(fmid * f2 >= 0.d0) THEN
            x2 = xmid; f2 = fmid
         ELSE
            print *, "2175: ERROR!!"
            stop
         ENDIF

         IF(jj > 1000) THEN
            print *, "2177: itteration exceeds 1000!"
            stop
         ENDIF
      ENDDO
      rffrc2 = 20.67d0 * xmid * xmid ! free fall radius
      ! xmid = r / rs; r/rc = xmid / 0.22

      1111        continue

      ! --- Cooling radius is set is the minimum value
      ! --- of cooling radius, free-fall radius and virial radius
      rrc1 = sqrt(min(rrc2, rffrc2))
      rrc1 = min(rrc1, cst/0.22d0)
      dMchcool = Bc * (rrc1 - atan(rrc1))
      rrs1 = 0.22d0 * sqrt(rrc2)
             ! r_cool / rs; for free fall time measurament
             ! even if rrc2 > rffrc2, rrs1 is set to sqrt(rrc2)
             ! (Shirakata; 2016/Nov/15)
      rrs1 = min(rrs1, cst)
      xcoolvir = rrc1 * 0.22d0 / cst
      xcoolvir = min(xcoolvir, 1.d0)         ! Rvir / rc = cst / 0.22
                                             ! need for disk size measurament
                                             ! (Shirakata; 2016/Nov/15)
   ENDIF
!!$     ==========================
   IF(dMch > 1.d0) dMch = 1.d0
   IF(dMchcool > 1.d0) dMchcool = 1.d0
!!$        === z-depENDently evolution type ===
!!$            IF((zp1now < zonplus1) .and. (zp1now >= zoffplus1)) THEN
!!$               fuv_pre = 1.32d+4 * Square(param%mu) * mrgt%mhalo**(5.d0/3.d0) &
!!$                                 * param%h**(-10.d0/3.d0) / Square(Mhotini)
!!$               fuv = fuv_pre * (1.d0 / dens(zp1now))**(5.d0/3.d0)
!!$               IF(zp1now < 3.d0) THEN
!!$                  dMch = dMchcool / (1.d0 + Cube(dMchcool) * i21 &
!!$                         (zp1now / 3.d0)**4 * fuv)
!!$               ELSE
!!$                  dMch = dMchcool / (1.d0 + Cube(dMchcool) * i21 &
!!$                         (1.d0 / (zp1now / 3.d0)) * fuv)
!!$               ENDIF
!!$            ELSE
   dMch = dMchcool

!!$            ENDIF
!!$     ==========================
1112 continue
END SUBROUTINE cool1
!!$ ========================================================================
!!$     ---------------------- SUBROUTINE lum --------------
SUBROUTINE lum(zp1pc, zp1, lumtmp, mz, mt, tel)
   use global_var; use SFHrelated
   implicit none

   DOUBLE PRECISION, INTENT(IN) :: zp1pc, zp1, tel
   DOUBLE PRECISION, INTENT(INOUT) :: lumtmp(param%nwave), mz, mt
   DOUBLE PRECISION, PARAMETER :: EPS  = 1.d-8
   DOUBLE PRECISION, PARAMETER :: EPS2 = 1.d-14
   INTEGER :: i, tt, zli, k, ttp1, zlip1
   DOUBLE PRECISION :: tstart, tfinal, tstart_lk, tfinal_lk, ttm1, ttm2, ttm3
   DOUBLE PRECISION :: tauz, tstep, tbranch
   DOUBLE PRECISION :: zt1, zt2, zt, iz
   DOUBLE PRECISION :: sfrt, sfrtt
   DOUBLE PRECISION :: fxyi0, fxyi1, fxy(param%nwave)
   DOUBLE PRECISION :: tout, tauinv, cc, cc_taust
   DOUBLE PRECISION :: tmp
   INTEGER, PARAMETER :: b_or_q = 2
   ! --- functions
   INTEGER :: IntTimeSSP, IntMetalSSP
   DOUBLE PRECISION :: z2t

   lumtmp(:) = 0.d0; mz = 0.d0; mt = 0.d0
   IF(param%SFH) mSFH(:,:) = 0.d0
   IF(zp1pc == zp1) RETURN

   sfrtt = 0.d0
   tauinv = ssp%p(b_or_q) / taust
   tauz   = 1.d0 / tauinv
   tout   = const%tout; tstart = z2t(zp1) + tel; tfinal = z2t(zp1pc)
   tstart_lk = tout - tstart; tfinal_lk  = tout - tfinal
   tbranch = tstart_lk - tfinal_lk

   cc = 1.d0 + beta - ssp%R(b_or_q); cc_taust = cc / taust

   iz = 1; zli = 1; ttm1 = 0.d0; ttm2 = 0.d0; zt2 = Zc0
   DO WHILE(abs(tbranch - ttm2) > EPS)
      tt    = IntTimeSSP(tstart_lk - ttm1)
      tstep = ssp%time(tt + 1) - ssp%time(tt)

      IF(ttm2 + tstep > tbranch) tstep = tbranch - ttm2
      ttm1 = ttm2; ttm2 = ttm2 + tstep

      zt1 = zt2; zt2 = Zc0 + ttm2 * tauinv; tmp = zt1 * zt2; zt = 1.d-8
      IF(tmp > EPS) zt = sqrt(tmp)

      IF(zt > ssp%y(b_or_q)) zt = ssp%y(b_or_q)
      sfrt = exp(-cc_taust * ttm1) * (1.d0 - exp(-cc_taust * tstep)) / cc
      IF(sfrt <= EPS2) goto 20
      sfrtt = sfrtt + sfrt

      ttm3 = tstart_lk - ttm1 - 0.5d0 * tstep
      zli  = IntMetalSSP(zt); zlip1 = zli + 1
      tt   = IntTimeSSP(ttm3); ttp1 = tt + 1
      tmp  = (ttm3 - ssp%time(tt)) / (ssp%time(ttp1) - ssp%time(tt))
      DO k = 1, param%nwave
         i = param%iwave(k)
         fxyi0 = ssp%lumi(b_or_q, i, tt, zli) &
            + (ssp%lumi(b_or_q, i, ttp1, zli)   - ssp%lumi(b_or_q, i, tt, zli))   * tmp
         fxyi1 = ssp%lumi(b_or_q, i, tt, zlip1) &
            + (ssp%lumi(b_or_q, i, ttp1, zlip1) - ssp%lumi(b_or_q, i, tt, zlip1)) * tmp
         !!$            fxy(k) = (fxyi1 / fxyi0)**(zl - dble(zli)) * fxyi0
         IF(zt > ssp%chem(NCHEM)) THEN
            fxy(k) = fxyi1
         ELSEIF(zt < ssp%chem(1)) THEN
            fxy(k) = fxyi0
         ELSE
            !!$               fxy(k) = (fxyi1 / fxyi0)**(zl - dble(zli)) * fxyi0
            fxy(k) = (fxyi1 - fxyi0) * (zt - ssp%chem(zli)) &
               / (ssp%chem(zlip1) - ssp%chem(zli)) + fxyi0
         ENDIF
         lumtmp(k) = lumtmp(k) + sfrt * fxy(k)
      ENDDO
      tmp = sfrt * fxy(param%iVband_r)
      mz = mz + tmp * zt; mt = mt + tmp * (tstart + ttm1)
   ENDDO
   20 continue

   IF(param%SFH) call CalcSFH(b_or_q, tstart, tfinal, taust)
END SUBROUTINE lum
!!$     ---------------------- SUBROUTINE lum_burst --------------
SUBROUTINE lum_burst(zp1pc, zp1, lumtmp, mz, mt, tel)
   use global_var; use SFHrelated
   implicit none

   DOUBLE PRECISION, INTENT(IN) :: zp1pc, zp1, tel
   DOUBLE PRECISION, INTENT(INOUT) :: lumtmp(param%nwave), mz, mt
   DOUBLE PRECISION, PARAMETER :: EPS  = 1.d-8
   DOUBLE PRECISION, PARAMETER :: EPS2 = 1.d-14
   INTEGER :: i, tt, zli, k, ttp1, zlip1
   DOUBLE PRECISION :: tstart, tfinal, tstart_lk, tfinal_lk, ttm1, ttm2, ttm3
   DOUBLE PRECISION :: tauz, tstep, tbranch
   DOUBLE PRECISION :: zt1, zt2, zt, iz
   DOUBLE PRECISION :: sfrt, sfrtt
   DOUBLE PRECISION :: fxyi0, fxyi1, fxy(param%nwave)
   DOUBLE PRECISION :: tout, tauinv, cc, cc_taust
   DOUBLE PRECISION :: tmp
   INTEGER, PARAMETER :: b_or_q = 1
   ! --- functions
   INTEGER :: IntTimeSSP, IntMetalSSP
   DOUBLE PRECISION :: z2t

   lumtmp(:) = 0.d0; mz = 0.d0; mt = 0.d0
   IF(param%SFH) mSFH(:,:) = 0.d0
   IF(zp1pc == zp1) RETURN

   sfrtt = 0.d0
   tauinv = ssp%p(b_or_q) / taust
   tauz   = 1.d0 / tauinv
   tout   = const%tout; tstart = z2t(zp1) + tel; tfinal = z2t(zp1pc)
   tstart_lk = tout - tstart; tfinal_lk  = tout - tfinal
   tbranch = tstart_lk - tfinal_lk

   cc = 1.d0 + beta - ssp%R(b_or_q); cc_taust = cc / taust

   iz = 1; zli = 1; ttm1 = 0.d0; ttm2 = 0.d0; zt2 = Zc0
   DO WHILE(abs(tbranch - ttm2) > EPS)
      tt    = IntTimeSSP(tstart_lk - ttm1)
      tstep = ssp%time(tt + 1) - ssp%time(tt)

      IF(ttm2 + tstep > tbranch) tstep = tbranch - ttm2
      ttm1 = ttm2; ttm2 = ttm2 + tstep

      zt1 = zt2; zt2 = Zc0 + ttm2 * tauinv; tmp = zt1 * zt2; zt = 1.d-8
      IF(tmp > EPS) zt = sqrt(tmp)

      IF(zt > ssp%y(b_or_q)) zt = ssp%y(b_or_q)


      sfrt = exp(-cc_taust * ttm1) * (1.d0 - exp(-cc_taust * tstep)) / cc
      IF(sfrt <= EPS2) goto 20
      sfrtt = sfrtt + sfrt

      ttm3 = tstart_lk - ttm1 - 0.5d0 * tstep
      zli  = IntMetalSSP(zt); zlip1 = zli + 1
      tt   = IntTimeSSP(ttm3); ttp1 = tt + 1
      tmp  = (ttm3 - ssp%time(tt)) / (ssp%time(ttp1) - ssp%time(tt))
      DO k = 1, param%nwave
         i = param%iwave(k)
         fxyi0 = ssp%lumi(b_or_q, i, tt, zli) &
            + (ssp%lumi(b_or_q, i, ttp1, zli)   - ssp%lumi(b_or_q, i, tt, zli))   * tmp
         fxyi1 = ssp%lumi(b_or_q, i, tt, zlip1) &
            + (ssp%lumi(b_or_q, i, ttp1, zlip1) - ssp%lumi(b_or_q, i, tt, zlip1)) * tmp

         IF(zt > ssp%chem(NCHEM)) THEN
            fxy(k) = fxyi1
         ELSEIF(zt < ssp%chem(1)) THEN
            fxy(k) = fxyi0
         ELSE
            fxy(k) = (fxyi1 - fxyi0) * (zt - ssp%chem(zli)) &
               / (ssp%chem(zlip1) - ssp%chem(zli)) + fxyi0
         ENDIF
         lumtmp(k) = lumtmp(k) + sfrt * fxy(k)
      ENDDO
      tmp = sfrt * fxy(param%iVband_r)
      mz = mz + tmp * zt; mt = mt + tmp * (tstart + ttm1)
   ENDDO
   20 continue

   IF(param%SFH) call CalcSFH(b_or_q, tstart, tfinal, taust)
END SUBROUTINE lum_burst
!!$     ---------------------- SUBROUTINE lumagn --------------
SUBROUTINE lumagn(itime, id, zp1, zp1pc, tel, lumqtmp, agntmp, memtmp)
   ! --- calculating the QSO bolometric luminosity
   ! --- Revised by Shirakata (2016/Aug/11)
   ! ----- Mdot = const. = dMbh / tacc

   ! --- DATA STRUCTURE (agntmp)--- !
   ! ----- must update for all output time.
   !
   ! agntmp(1) : Lbol [10^12 Lsun], at output time
   ! agntmp(2) : Lbol [10^12 Lsun], peak luminosity
   ! agntmp(3) : Macc_eff [10^14 Msun]
   ! agntmp(4) : Mcold [10^14 Msun], at tstart
   ! agntmp(5) : MZcold [10^14 Msun], at tstart
   ! agntmp(6) : t_ad (accretion disk timescale) [yr]
   ! agntmp(7) : t_gal (bulge dynamical time) [Gyr]
   ! agntmp(10) : Mrem [10^14 Msun]
   !
   use global_var
   implicit none
   INTEGER, INTENT(in)  :: id, itime
   DOUBLE PRECISION, INTENT(in) :: zp1, zp1pc, tel
   DOUBLE PRECISION, INTENT(out) :: lumqtmp(param%nwaveq)
   DOUBLE PRECISION, INTENT(out) :: agntmp(10)
   DOUBLE PRECISION, INTENT(out) :: memtmp(10)
   DOUBLE PRECISION, PARAMETER :: EPS = 1.d-30
   DOUBLE PRECISION :: tstart, tout, tage, tacc, tout2, tage2, tdelay
   DOUBLE PRECISION :: t_ad, t_gal ! dyn. time of accretion disk and bulge (Shirakata 16/09/21)
   DOUBLE PRECISION :: Mdotpeak, Mdot
   DOUBLE PRECISION :: Lbol, Lbol2, log10Lbol, Lbolpeak
   DOUBLE PRECISION :: Ledd0, Ledd1, Ledd2, Medd0, Medd1, Medd2 !@ tstart, tout, tout2
   DOUBLE PRECISION :: Mdot2, dMbh2, dMbh1, lumqtmp2B, lumqtmp2HX, lumqtmp2SX
   DOUBLE PRECISION :: dMbhrem, Macc_eff, Mbhinit
   DOUBLE PRECISION :: tmp
   DOUBLE PRECISION :: z2t ! functions
   DOUBLE PRECISION :: Edd_mratio, Edd_lratio     ! mass, luminosity
   DOUBLE PRECISION :: Edd_mratiop, Edd_lratiop   ! mass, luminosity
   DOUBLE PRECISION :: Edd_mratio2, Edd_lratio2   ! mass, luminosity
   DOUBLE PRECISION :: Mdot2Lbol, MdotLbol2L_HX, MdotLbol2L_SX, MdotLbol2L_B !functions

   Mbhinit  = gal(id)%Mbh - gal(id)%agn(10)
   Macc_eff = dMbh + gal(id)%agn(10)
   lumqtmp(:) = 0.d0; agntmp(:) = 0.d0
   memtmp(:)  = 0.d0

   t_gal = 0.9523 * gal(id)%rbulge / gal(id)%Vbulge / param%h ! [Gyr]

   t_ad  = 0.d0
   IF(Macc_eff * param%munit > param%Mbhseed) THEN
      t_ad  = param%tad_0 &
              *((Macc_eff *param%munit)** param%alpha_ad) &
              *(Mbhinit*param%munit) ** param%gamma_ad
   ENDIF

   tacc  = param%tacc_0 * t_gal + param%tacc_1 * t_ad ! [Gyr]
   ! tacc   = 3e-2 * zp1 **(-1.5d0)   ! Kauffmann & Haenelt 00
   tdelay = 0.d0
   IF(param%delay) tdelay = t_gal / param%th    ![hubble time]

   tout  = const%tout; tstart = z2t(zp1) + tel + tdelay ! [hubble time]
   tout2 = z2t(zp1pc) ! [hubble time]
   tage  = max((tout  - tstart), 0.d0) * param%th ! the elapsed time from the onset of burst to
   ! the output time [Gyr]
   tage2 = max((tout2 - tstart), 0.d0) * param%th ! the elapsed time from the onset of burst to
   ! the end of this time slice [Gyr]

   IF(tacc <= 13.8) THEN
      Mdotpeak =  Macc_eff * 1.d+5 / tacc
      Mdot = 0.d0; Mdot2 = 0.d0

      IF(tage >= param%acc * tacc) THEN
         tmp   = exp(-(tage - param%acc * tacc) / ((1.d0 - param%acc) * tacc))
         Mdot  = Mdotpeak * tmp
         dMbh1 = Macc_eff * (1.d0 - param%acc) * tmp
         tmp = 0.d0
         IF(tage2 >= param%acc * tacc) THEN
            Mdot2 = 0.d0; dMbh2 = 0.d0
           IF(param%acc < 1.d0) THEN
              tmp = exp(-(tage2 - param%acc * tacc) / ((1.d0 - param%acc) * tacc))
              Mdot2 = Mdotpeak * tmp
              dMbh2 = Macc_eff * (1.d0 - param%acc) * tmp
           ENDIF
         ELSE
            Mdot2 = Mdotpeak
            dMbh2 = Macc_eff * (1.d0 - tage2 / tacc)
         ENDIF
      ELSE
         Mdot  = Mdotpeak; Mdot2 = Mdotpeak
         dMbh2 = Macc_eff * (1.d0 - tage2 / tacc)
         dMbh1 = Macc_eff * (1.d0 - tage / tacc)
      ENDIF
   ELSE
      Mdotpeak = 0.d0; Mdot = 0.d0; Mdot2 = 0.d0
      dMbh1 = Macc_eff; dMbh2 = Macc_eff
   ENDIF

   dMbhrem = dMbh2
   gal(id)%Mbh = gal(id)%Mbh + dMbh

   ! --- Eddington limit @ t = tout --- !
   Ledd1 = 2.6d-46 * param%Ledd0 * (gal(id)%Mbh - dMbh1) ! at tout [10^12 Lsun]
   Medd1 = 6.83d-2 * Ledd1

   Ledd2 = 2.6d-46 * param%Ledd0 * (gal(id)%Mbh - dMbh2) ! at tout2 [10^12 Lsun]
   Medd2 = 6.83d-2 * Ledd2

   Ledd0 = 2.6d-46 * param%Ledd0 * Mbhinit               ! at tstart [10^12 Lsun]
   Medd0 = 6.83d-2 * Ledd0
   ! Ledd / c^2
   ! critical accretion rate [Msun/yr]
   ! at tout, Msun/yr = 14.9d0 * 10^12 Lsun (same as Mdot2Lbol)

   Edd_mratio = 0.d0
   IF(Medd1 > 0.d0) THEN ! at t = tout
      Edd_mratio = max(Mdot / Medd1, EPS)
   ENDIF
   Edd_mratio2 = 0.d0
   IF(Medd2 > 0.d0) THEN ! at t = tout2
      Edd_mratio2 = max(Mdot2 / Medd2, EPS)
   ENDIF
   Edd_mratiop = 0.d0
   IF(Medd0 > 0.d0) THEN ! at t = tstart
      Edd_mratiop = max(Mdotpeak / Medd0, EPS)
   ENDIF

   Lbol  = Mdot2Lbol(Edd_mratio,  Mdot,  Ledd1) ! bolometric luminosity [10^12 Lsun]
   Lbol2 = Mdot2Lbol(Edd_mratio2, Mdot2, Ledd2) ! bolometric luminosity [10^12 Lsun]
   Lbolpeak = Mdot2Lbol(Edd_mratiop, Mdotpeak,Ledd0) ! peak bolometric
                                                        ! luminosity
   Edd_lratio  = Lbol / Ledd1
   Edd_lratiop = Lbolpeak / Ledd0
   Edd_lratio2 = Lbol2 / Ledd2

   ! --- calculate luminosity (HX, SX, B)
   IF(Lbol > 1.d-12) THEN
      log10Lbol  = log10(Lbol) ! L used in eq.(21) of Marconi+04
      lumqtmp(1) = MdotLbol2L_B(Lbol,  log10Lbol) ! Mdot(   B  ) [erg/s]
      lumqtmp(2) = 1.d0                           ! Mdot(  UV  ) ! not necessary
      lumqtmp(3) = 1.d0                           ! Mdot(  UV  ) ! not necessary
      lumqtmp(4) = MdotLbol2L_HX(Lbol, log10Lbol,Edd_mratio,Edd_lratio) ! Mdot(hard X) [erg/s]
      lumqtmp(5) = MdotLbol2L_SX(Lbol, log10Lbol) ! Mdot(soft X) [erg/s]
      lumqtmp(6) = Lbol * 3.826d+45               ! bolometric [erg/s]

      lumqtmp(7) = Edd_mratio                     ! Eddington ratio (mass)
      lumqtmp(8) = Edd_lratio                     ! Eddington ratio (lum)
      lumqtmp(9) = gal(id)%Mbh                    ! BH mass [Msun]
   ELSE
      log10Lbol  = 0.d0 ! L used in eq.(21) of Marconi+04
      lumqtmp(1:9) = 0.d0
   ENDIF



   IF(Lbol2 > 1.d-12) THEN
      log10Lbol  = log10(Lbol2) ! L used in eq.(21) of Marconi+04
      lumqtmp2B  = MdotLbol2L_B(Lbol2,  log10Lbol) ! Mdot(B-band) [erg/s]
      lumqtmp2HX = MdotLbol2L_HX(Lbol2, log10Lbol,Edd_mratio2, Edd_lratio2) ! Mdot(hard X) [erg/s]
      lumqtmp2SX = MdotLbol2L_SX(Lbol2, log10Lbol) ! Mdot(soft X) [erg/s]
                   ! absolute magnitude in B-band at tout2 [ABmag]
                   ! not M_B-5log10h
   ELSE
      log10Lbol  = 0.d0 ! L used in eq.(21) of Marconi+04
      lumqtmp2B  = 0.d0 ! L_B  [Msun/yr]
      lumqtmp2HX = 0.d0 ! L_HX [erg/s]
      lumqtmp2SX = 0.d0 ! L_SX [erg/s]
   ENDIF

   ! assign agntmp
   agntmp(1) = Lbol               ! Lbol (at tout) [10^12 Lsun]
   agntmp(2) = Lbolpeak           ! Lbol (at highest) [10^12 Lsun]
   agntmp(3) = Macc_eff
   agntmp(4) = Mc0  * param%munit ! cold gas mass exhausted by an SB
   agntmp(5) = MZc0 * param%munit ! heavy element mass in cold gas
   agntmp(6) = t_ad
   agntmp(7) = t_gal              ! [Gyr]
   agntmp(10) = dMbhrem            ! [10^14 Msun]

   IF(gal(id)%Lpeak < Lbolpeak) gal(id)%Lpeak = Lbolpeak
   !Added by Shirakata (2017/10/29)

   IF(param%traceIDs == 1) THEN
      memtmp(1) = Lbolpeak
      memtmp(2) = Lbol2
      memtmp(3) = Edd_mratiop
      memtmp(4) = Edd_mratio2
      memtmp(5) = Edd_lratiop
      memtmp(6) = Edd_lratio2
      memtmp(7) = Mbhinit
      memtmp(8) = Macc_eff
      memtmp(9) = tstart           ! hubble time
      memtmp(10) = tacc / param%th ! hubble time
   ENDIF

END SUBROUTINE lumagn
!!$     -------------------FUNCTION Mdot2Lbol--------------------
DOUBLE PRECISION FUNCTION Mdot2Lbol(EddRatio, Mdot, Ledd) RESULT(Lbol) 
   ! --- converting accretion rate of Mdot [Msun/yr] into bolometric luminosity
   !       of Lbol [10^12 Lsun]
   use global_var
   implicit none
   DOUBLE PRECISION, INTENT(IN) :: EddRatio, Mdot, Ledd
   DOUBLE PRECISION :: A, B, C
   DOUBLE PRECISION :: Lbol0 = 14.9d0 ! conversion factor from [Msun/yr]
   DOUBLE PRECISION :: EPS = 1.d-30

   Lbol = 0.d0 ! modified by MARK (2017/Mar/21)
   IF(EddRatio > EPS) THEN
      IF(param%ML_AGN == 1) THEN
         Lbol = Lbol0 * param%eps_agn * Mdot ! [10^12 Lsun]
      ELSEIF(param%ML_AGN == 2) THEN
         IF(EddRatio >= 20.d0) THEN
            Lbol = Ledd * 2.d0 * (1.d0 + log(EddRatio / 20.d0))  ! [10^12 Lsun]
         ELSE IF(0.03d0 < EddRatio .and. EddRatio < 20.d0) THEN
            Lbol = Ledd  * EddRatio / 10.d0                      ! [10^12 Lsun]
         ELSE IF(EddRatio < 0.03d0) THEN
            Lbol = Ledd * 10.d0 * EddRatio * EddRatio/ 3.d0 ! [10^12 Lsun]
         ENDIF
      ELSEIF(param%ML_AGN == 3) THEN ! T. Kawaguchi in private communication
         C = dtanh(dlog10(EddRatio/param%Edd_max))
         A = 1.d0 / (1.d0 + 3.5d0 * (1.d0 + C))
         B = 1.d0 / (EddRatio / param%Edd_max)
         Lbol = Ledd * 1.d0 / (A + B)
!         ! Added by Shirakata (2018/Aug/12)
!         IF(EddRatio < 0.01 * param%Edd_max) THEN
!            Lbol = Ledd * 0.01d0 * (100.d0 * EddRatio / param%Edd_max)**2.d0
!         ENDIF
      ENDIF
   ENDIF
END FUNCTION Mdot2Lbol
!!!$     -------------------FUNCTION MdotLbol2L_HX--------------------
DOUBLE PRECISION FUNCTION MdotLbol2L_HX(Lbol, log10Lbol, mdot, lamb) RESULT(L_HX)
   ! --- converting accretion rate Mdot [Msun/yr] and log10Lbol = log10(Lbol)-12
   !       into Mdot(2--10 keV) [Msun/yr] using a third-degree polynomial relation
   !       between Lbol and L(2--10 keV) provided in Marconi+04 [eq.(21)]
   ! --- log10[Lbol/L(2--10 keV)] = 1.54 + 0.24L + 0.012L^2 - 0.0015L^3
   !       where L = log10(Lbol[Lsun]) - 12
   use global_var
   DOUBLE PRECISION, INTENT(IN) :: Lbol ! [10^12 Lsun]
   DOUBLE PRECISION, INTENT(IN) :: log10Lbol ! L used in eq.(21) of Marconi+04
   DOUBLE PRECISION, INTENT(IN) :: mdot, lamb 
   DOUBLE PRECISION :: frac_HX ! = log10[Lbol / L(2--10 keV)]
   DOUBLE PRECISION :: coeff(4) = (/1.54d0, 0.24d0, 0.012d0, -0.0015d0/)
   DOUBLE PRECISION :: Lsun2L0 = 3.826d+45 ! conversion factor from [10^12 Lsun] to [erg/s]
   IF(param%BolCor == 1) THEN ! Marconi+ 04
      frac_HX = coeff(1) + (coeff(2) + (coeff(3) + coeff(4) * log10Lbol) * log10Lbol) &
                           * log10Lbol
      L_HX = Lbol * Lsun2L0 / 10.d0**frac_HX ! [erg/s]
!      IF(mdot < 0.01d0 * param%Edd_max) THEN
!         L_HX = 200.d0 * lamb * 10.d0**frac_HX * L_HX
!         print *, 200.d0 * lamb * 10.d0**frac_HX
!      ENDIF
   ENDIF
END FUNCTION MdotLbol2L_HX
!!!$     -------------------FUNCTION MdotLbol2L_SX--------------------
DOUBLE PRECISION FUNCTION MdotLbol2L_SX(Lbol, log10Lbol) RESULT(L_SX)
   ! --- converting accretion rate Mdot [Msun/yr] and log10Lbol = log10(Lbol)-12
   !       into Mdot(0.5--2 keV) [Msun/yr] using a third-degree polynomial relation
   !       between Lbol and L(0.5--2 keV) provided in Marconi+04 [eq.(21)]
   ! --- log10[Lbol/L(0.5--2 keV)] = 1.65 + 0.22L + 0.012L^2 - 0.0015L^3
   !       where L = log10(Lbol[Lsun]) - 12
   use global_var
   DOUBLE PRECISION, INTENT(IN) :: Lbol ! [10^12 Lsun]
   DOUBLE PRECISION, INTENT(IN) :: log10Lbol ! L used in eq.(21) of Marconi+04
   DOUBLE PRECISION :: frac_SX ! = log10[Lbol / L(0.5--2 keV)]
   DOUBLE PRECISION :: coeff(4) = (/1.65d0, 0.22d0, 0.012d0, -0.0015d0/)
   DOUBLE PRECISION :: Lsun2L0 = 3.826d+45 ! conversion factor from [Msun/yr] to [erg/s]
   IF(param%BolCor == 1) THEN ! Marconi+ 04
      frac_SX = coeff(1) + (coeff(2) + (coeff(3) + coeff(4) * log10Lbol) * log10Lbol) &
                           * log10Lbol
      L_SX = Lbol * Lsun2L0 / 10.d0**frac_SX ! [erg/s]
   ENDIF
END FUNCTION MdotLbol2L_SX
!!!$     -------------------FUNCTION MdotLbol2L_B--------------------
DOUBLE PRECISION FUNCTION MdotLbol2L_B(Lbol, log10Lbol) RESULT(L_B)
   ! --- converting accretion rate Mdot [Msun/yr] and log10Lbol = log10(Lbol)-12
   !       into normalized B-band luminosity normLB using a third-degree polynomial
   !       relation between Lbol and nu_B L_B provided in Marconi+04 [eq.(21)]
   ! --- log10[Lbol/nu_B L_B] = 0.80 - 0.067L + 0.017L^2 - 0.0023L^3
   !       where L = log10(Lbol[Lsun]) - 12
   use global_var
   DOUBLE PRECISION, INTENT(IN) :: Lbol ! [Msun/yr]
   DOUBLE PRECISION, INTENT(IN) :: log10Lbol ! L used in eq.(21) of Marconi+04
   DOUBLE PRECISION :: frac_B ! = log10[Lbol / nu_B L_B]
   DOUBLE PRECISION :: coeff(4) = (/0.80d0, -0.067d0, 0.017d0, -0.0023d0/)
   DOUBLE PRECISION :: Lsun2L0 = 3.826d+45 ! conversion factor from [Msun/yr] to [erg/s]
   IF(param%BolCor == 1) THEN ! Marconi+ 04
      frac_B = coeff(1) + (coeff(2) + (coeff(3) + coeff(4) * log10Lbol) &
                                      * log10Lbol) * log10Lbol
      L_B = Lbol * Lsun2L0 / 10.d0**frac_B ! [erg/s]
   ENDIF
END FUNCTION MdotLbol2L_B
!!$     -------------------------------------------------------
INTEGER FUNCTION IntTimeSSP(time) RESULT(tt)
   use global_var
   implicit none
   DOUBLE PRECISION, INTENT(IN) :: time

   IF(param%type_ssp == 2) THEN ! PEGASE
      IF(time > ssp%time(500)) THEN
         tt = 500
      ELSEIF(time > ssp%time(400)) THEN
         tt = 400
      ELSEIF(time > ssp%time(300)) THEN
         tt = 300
      ELSEIF(time > ssp%time(200)) THEN
         tt = 200
      ELSEIF(time > ssp%time(100)) THEN
         tt = 100
      ELSE
         tt = 1
      ENDIF
   ELSEIF(param%type_ssp == 3) THEN ! BC03
      IF(time > ssp%time(110)) THEN
         tt = 110
      ELSE
         tt = 1
      ENDIF
   ELSEIF(param%type_ssp == 1) THEN ! KA97 (original)
      tt = 1
   ENDIF

   DO WHILE(time > ssp%time(tt+1) .and. tt < NSFR)
   tt = tt + 1
   ENDDO
   tt = min(tt, NSFR - 1)
END FUNCTION IntTimeSSP
!!$     -------------------------------------------------------
INTEGER FUNCTION IntMetalSSP(metal) RESULT(zli)
   use global_var
   implicit none
   DOUBLE PRECISION, INTENT(IN) :: metal

   IF(metal < ssp%chem(1)) THEN
      zli = 1
   ELSEIF(metal > ssp%chem(NCHEM)) THEN
      zli = NCHEM - 1
   ELSE
      zli = 1
      DO WHILE(zli < NCHEM .and. metal > ssp%chem(zli+1))
      zli = zli + 1
      ENDDO
      zli = min(zli, NCHEM - 1)
   ENDIF
END FUNCTION IntMetalSSP
!!$     ------------------ FUNCTION gasdev ---------------------------
REAL FUNCTION gasdev(idum)
   use global_var

   INTEGER :: idum
   REAL :: fac, gset, rsq, v1, v2, ran1
   SAVE :: gset

   IF(param%iset == 0) THEN
      1    v1 = 2. * (ran1(idum)) - 1.
      v2 = 2. * (ran1(idum)) - 1.
      rsq = v1**2 + v2**2
      IF(rsq >= 1. .or. rsq == 0.) goto 1
      fac  = sqrt(-2. * log(rsq) / rsq)
      gset = v1 * fac
      gasdev = v2 * fac
      param%iset=1
   ELSE
      gasdev=gset
      param%iset=0
   ENDIF
   return
END FUNCTION gasdev
!!$     ------------------- FUNCTION ran1 ------------------------
REAL FUNCTION ran1(idum)
   use global_var

   INTEGER :: idum
   INTEGER, PARAMETER :: IA = 16807
   INTEGER, PARAMETER :: IM = 2147483647
   INTEGER, PARAMETER :: IQ = 127773
   INTEGER, PARAMETER :: IR = 2836
   INTEGER, PARAMETER :: NTAB = 32
   INTEGER, PARAMETER :: NDIV = 1 + (IM - 1) / NTAB
   INTEGER :: j, k
   REAL, PARAMETER :: AM = 1. / IM
   REAL, PARAMETER :: EPS = 1.2e-7
   REAL, PARAMETER :: RNMX = 1. - EPS


   IF(idum <= 0 .or. param%iy == 0) THEN
      idum = max(-idum, 1)
      DO j = NTAB + 8, 1, -1
      k = idum / IQ
      idum = IA * (idum - k * IQ) - IR * k
      IF(idum < 0) idum = idum + IM
      IF(j <= NTAB) param%iv(j) = idum
      ENDDO
      param%iy = param%iv(1)
   ENDIF
   k = idum / IQ
   idum = IA * (idum - k * IQ) - IR * k
   IF(idum < 0) idum = idum + IM
   j = 1 + param%iy / NDIV
   param%iy = param%iv(j)
   param%iv(j) = idum
   ran1 = min(AM * param%iy, RNMX)

   return
END FUNCTION ran1
!!$============================================================================
DOUBLE PRECISION FUNCTION Square(x)
   implicit none
   DOUBLE PRECISION, INTENT(IN) :: x

   Square = x * x
END FUNCTION Square
!!$============================================================================
DOUBLE PRECISION FUNCTION Cube(x)
   implicit none
   DOUBLE PRECISION, INTENT(IN) :: x

   Cube = x * x * x
END FUNCTION Cube
!!$ ========================================================================
INTEGER FUNCTION Factorial(n)
   implicit none
   INTEGER, INTENT(IN) :: n
   INTEGER :: i, A
   A = 1
   DO i = 1, n
      A = A * i
   END DO
   Factorial = A
END FUNCTION Factorial
!!$ ========================================================================
SUBROUTINE FracDI(Mgas, Mstar, Mgal, fg, fs, fgrem, fsrem)
   use global_var
   implicit none
   DOUBLE PRECISION, INTENT(IN) :: Mgas, Mstar, Mgal
   DOUBLE PRECISION, INTENT(INOUT) :: fg, fs, fgrem, fsrem
   DOUBLE PRECISION :: Mdisk, fgas, fdisk, RgRd
   Mdisk = Mgas + Mstar
   IF(Mdisk > 0.d0) THEN
      fdisk = Mdisk / Mgal
      fgas = MIN(Mgas / Mdisk, 1.d0)
      RgRd = (1.d0-fgas) * fdisk * param%fdi
      fg = (1.d0 - (1.d0 + RgRd) * exp(-RgRd)) ! same as the merger case
      fs = param%fdi
   ELSE
      fg = 0.d0; fs = 0.d0
   ENDIF
   fgrem = 1.d0 - fg
   fsrem = 1.d0 - fs 
END SUBROUTINE FracDI
!!$ ========================================================================
!!$  Fitting function for the characteristic mass of Mc(z)
!!$  --- Mc(z) = 6.5e9 * exp(-0.604*z) * exp[-(z / 8.37)**17.6]
!!$      --- functional form is given by Nagashima on 2014/03/07
!!$      --- numerical constants are determined via gnuplot fitting tool
!!$            at zreion=9 by MARK on 2014/03/12
!!$  CAUTION!!!: this fitting formula is valid only for param%zreion = 9
DOUBLE PRECISION FUNCTION Mc_Okamoto08(redshift) RESULT(Mc)
   use global_var
   implicit none
   DOUBLE PRECISION, INTENT(IN) :: redshift
   DOUBLE PRECISION :: tmp1, tmp2

   Mc = 0.d0
   IF(redshift <= param%zreion) THEN
      tmp1 = 0.604d0 * redshift; tmp2 = (redshift / 8.37d0)**17.6d0
      Mc = 6.5d+9 * exp(-tmp1) * exp(-tmp2) / param%munit ! [Msun] --> [10^14 Msun]
   ENDIF
END FUNCTION Mc_Okamoto08
!!$ ========================================================================
!!$ baryon fraction after reionization (Gnedin00 and Okamoto+08)
DOUBLE PRECISION FUNCTION fb_Okamoto08(Mhalo, redshift, Mc) RESULT(fb)
   use global_var
   implicit none
   DOUBLE PRECISION, INTENT(IN) :: Mhalo, redshift, Mc ! Mhalo,Mc [10^14 Msun]
   DOUBLE PRECISION :: tmp

   fb = param%bar_rat
   IF(redshift <= param%zreion) THEN
      tmp = (Mc / Mhalo)**const%alpha_UV
      fb  = param%bar_rat * (1.d0 + const%UV1 * tmp)**const%UV2
   ENDIF
END FUNCTION fb_Okamoto08
!!$ ========================================================================
!!$     === Return the Circular Velocity [km/s] by Inputting
!!$         Mass [param%munit] and density [rho_crit@z=0] ===
DOUBLE PRECISION FUNCTION CircVel(mass, density)
   use global_var
   DOUBLE PRECISION, INTENT(IN) :: mass, density
   DOUBLE PRECISION :: constV = 3.13d+2 ! [km/s]

   CircVel = constV * (param%h * mass)**const%V1 * density**const%V2
   ! const%V1 = 1/3, const%V2 = 1/6: substituted at "SetConstants" in main_nugc.f90
END FUNCTION CircVel
!!$ ========================================================================
DOUBLE PRECISION FUNCTION CalTidal(Mratio, Vratio)
   DOUBLE PRECISION, INTENT(IN) :: Mratio, Vratio
   DOUBLE PRECISION :: Cube

   CalTidal = min(1.d0, 0.2d0 * Mratio / Cube(Vratio))
END FUNCTION CalTidal
!!$ ========================================================================
!!$     === Return the Merging Timescale via Dynamical Friction
!!$         by Inputting Mhalo,cent/Mhalo,sat (=Mratio) and density [] ===
DOUBLE PRECISION FUNCTION CalTaumrg(Mratio, Mcent, density, zp1)
   use global_var
   DOUBLE PRECISION, INTENT(IN) :: Mratio, Mcent, density, zp1
   DOUBLE PRECISION :: eps    ! circularity of the satellite orbit
   DOUBLE PRECISION :: Mcrit  ! halo mass that collapsed at fixed redshift
   DOUBLE PRECISION :: Square ! function
   DOUBLE PRECISION :: z

   ! NY04 (Chandrasekar formula)
   !CalTaumrg = 1.17d0 * Mratio &
   !     / (log(1.d0 + Square(Mratio)) * sqrt(0.5d0 * density))

   ! circularity is determined by Mhalo and redshift
   ! the following derivation is obtained from Wetzel et al 2011
   z     = zp1 - 1.d0
   Mcrit = 10.d0**(-1.58d0 + z * (-1.56d0 + 0.038d0 * z)) / param%h ! [10^14 Msun]
   eps   = 1.05d0 / (1.05d0 + 0.242d0 * (1.d0 + 2.36d0 *(Mcent/Mcrit)**0.108d0))
   ! eps  = 0.5d0 ! mean value neglecting halo mass and redshift dependences

   ! Jiang et al. 2008, 2010
   CalTaumrg = (0.9d0 * eps**0.47d0 + 0.60) / 2.d0 / 0.43d0 &
               * Mratio / log(1.d0 + Mratio) / sqrt(0.5d0 * density)
END FUNCTION CalTaumrg
!!$ ========================================================================
!!$     === Return the Bulge Radius [kpc/h] by Inputting Mstar+Mcold (=Mtotal)
!!$         [param%munit] and Bulge Velocity [km/s] ===
DOUBLE PRECISION FUNCTION CalRbulgeMerger(id1, Mtotal, y1b, velocity)
   use global_var
   INTEGER, INTENT(IN) :: id1
   DOUBLE PRECISION, INTENT(IN) :: Mtotal, y1b, velocity
   DOUBLE PRECISION :: f_dm
   DOUBLE PRECISION :: Square ! function
   DOUBLE PRECISION :: FracDM
   f_dm = FracDM(id1)
   CalRbulgeMerger = const%Rb * (f_dm * y1b + Mtotal) / Square(velocity) ! kpc/h
   CalRbulgeMerger = 2.d0 * CalRbulgeMerger ! Shirakata (2017/Oct/20)
END FUNCTION CalRbulgeMerger
!!$ ========================================================================
!!$     === Return the Bulge Radius [kpc/h] by Inputting Mstar+Mcold (=Mtotal)
!!$         [param%munit] and Bulge Velocity [km/s] ===
DOUBLE PRECISION FUNCTION CalRbulgeDI(Mtotal, velocity)
   use global_var
   DOUBLE PRECISION, INTENT(IN) :: Mtotal, velocity
   DOUBLE PRECISION :: Square ! function

   CalRbulgeDI = const%Rb * Mtotal / Square(velocity) ! kpc/h
   CalRbulgeDI = 2.d0 * CalRbulgeDI ! Shirakata (2017/Oct/20)
END FUNCTION CalRbulgeDI
!!$ ========================================================================
!!$     === Return the Disk Radius [kpc/h] by Inputting Disk Circular
!!$         Velocity [km/s] and Density [] ===
DOUBLE PRECISION FUNCTION CalRdisk(velocity, density)
   DOUBLE PRECISION, INTENT(IN) :: velocity, density

   CalRdisk = 10.d0 * velocity * sqrt(2.d0 / density) ! kpc/h
END FUNCTION CalRdisk
!!$ ========================================================================
!!$     === Return the Bulge Velocity [km/s] by Inputting
!!$         the Total Energy of Each Galaxy (e1 and e2),
!!$         the Total Baryonic Mass of Each Galaxy (y1 and y2) [param%munit],
!!$         and Mstar+Mcold (=Mtotal) [param%munit]
DOUBLE PRECISION FUNCTION CalVbulgeMerger(id1, e1d, e1b, e1rem, e2d, e2b, y1b, y1, y2, fgas, massb)
   !!$     +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   !!$     E1 = -K1 = (M1b + Mbh1) * V1b**2 / 2 + (M1d + M1c) * V1d**2
   !!$     E2 = -K2 = (M2b + Mbh2) * V2b**2 / 2 + (M2d + M2c) * V2d**2
   !!$     Kint = M1*M2 * V**2 / (M1+M2) / 2,  M1 = M1b+M1d+M1c,  M2 = M2b+M2d+M2c
   !!$                                V: relative velocity
   !!$     E0 = -M0b * V0b**2 / 2 because of burst
   !!$     energy conservation: E0 = E1 + E2 + Kint
   !!$      IF V**2 = (Mlb * Vlb**2 + (Mld+Mlc) * Vlc**2) / Ml  (l: larger galaxy)
   !!$      THEN the below equation is O.K.
   !!$     +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   use global_var
   INTEGER, INTENT(IN) :: id1
   DOUBLE PRECISION, INTENT(IN) :: e1d, e1b, e1rem, e2d, e2b, y1b, y1, y2, fgas, massb 
   DOUBLE PRECISION :: xx, e_b, e_sum
   DOUBLE PRECISION :: e1, e2, mass1, mass2 
   DOUBLE PRECISION :: f_dm, fg
   DOUBLE PRECISION :: FracDM !function

   f_dm  = FracDM(id1)
   !fg    = fgas * (1e-3/MAX(gal(id1)%Mstard,1e-6))**0.2d0 
   fg = fgas
   mass1 = y1 + f_dm * y1b;  mass2 = y2 
   e1   = e1d + (1.d0 + f_dm) * e1b; e2 = e2d + e2b
   e_b = ((e1 + e1rem) * e2) / (mass2 * (e1 + e1rem)  / mass1 + mass1 * e2 / mass2)
   e_sum = e1 + e2 + e_b
   xx = (1.d0 + fg * param%fdiss) * e_sum / (massb + f_dm * y1b)
   CalVbulgeMerger = 1.d-4
   IF(xx > 1.d-8) CalVbulgeMerger = sqrt(xx)
END FUNCTION CalVbulgeMerger
!!$ ========================================================================
DOUBLE PRECISION FUNCTION CalVbulgeDI(ed, eb, yd, yb, fgas, massb)
   !!$     +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   !!$     From Virial Theorem,
   !!$       Ed = -Kd = (Md + Mc) * Vd**2
   !!$       Eb = -Kb = (Mb + Mbh) * Vb**2 
   !!$       E0 = -Mb * Vb**2  because of burst
   !!$     Energy Conservation: 
   !!$       E0 = Ed + Eb 
   !!$     Assuming V**2 = E0 / M
   !!$     +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   use global_var
   DOUBLE PRECISION, INTENT(IN) :: ed, eb, yd, yb, fgas, massb
   DOUBLE PRECISION :: xx, e_sum !, e_b

   !e_b = (ed * eb) / (yb * ed / yd + yd * eb / yb)
   e_sum = ed + eb !+ e_b
   xx = (1.d0+fgas*param%fdiss) * e_sum / massb
   CalVbulgeDI = 1.d-4
   IF(xx > 1.d-8) CalVbulgeDI = sqrt(xx)
END FUNCTION CalVbulgeDI
!!$ ========================================================================
!!$ Return the halo concentration parameter
DOUBLE PRECISION FUNCTION CalCST(zp1, Mhalo)
   use global_var
   DOUBLE PRECISION, INTENT(IN) :: zp1, Mhalo
   DOUBLE PRECISION :: Ez2, OmegaMz, x, Delta_c
   DOUBLE PRECISION :: Mhalo_tmp, c200
   INTEGER :: izp1, iM
   DOUBLE PRECISION :: Cube

   Ez2     = param%OMEGA_L + param%OMEGA0 * Cube(zp1)
   OmegaMz = 1.d0 - param%OMEGA_L / Ez2
   x       = OmegaMz - 1.d0
   Delta_c = 18.d0 * const%PI2 + (82.d0 - 39.d0 * x) * x
   CP%a = -1.119d0 * log10(Delta_c) + 3.537d0
   CP%b = -0.967d0 * log10(Delta_c) + 2.181d0

   Mhalo_tmp = Mhalo * param%h * param%munit ! Msun/h


   ! Bullock et al. 2001
   ! CalCST = (9.d0 / zp1) * (Mhalo_tmp/1.5d13)**(-0.13d0)

   ! Prada+ 2012
   izp1 = int((log10(zp1 ) -CP%zp1_base) / CP%zp1_step) + 1
   IF(izp1 < 1) izp1 = 1
   IF(izp1 > CP%NZ) izp1 = CP%NZ
   iM = int((log10(Mhalo_tmp) - CP%M_base) / CP%M_step) + 1
   IF(iM < 1) iM = 1
   IF(iM > CP%NM) iM = CP%NM
   c200 = CP%CST(izp1, iM)

   ! convert c200 to c100, see eq 16. of arXiv:1005.0411v1
   CalCST = CP%a * c200 + CP%b
END FUNCTION CalCST
!!$ ========================================================================
!!$ Return free fall radius related value
DOUBLE PRECISION FUNCTION FFR(x, delt, dens, cst, Vctmp, Mhalo) 
   use global_var
   DOUBLE PRECISION, INTENT(IN) :: x, delt, dens, cst, Vctmp, Mhalo
   DOUBLE PRECISION :: A, B ! coefficient
   DOUBLE PRECISION :: y, tmp, rs, Mdm
   DOUBLE PRECISION :: tff_s
   DOUBLE PRECISION :: Square, Cube, CalRdisk ! function

! --- x = r/r_s = 0.22 * (r/r_c)
   tmp   = 1.d0 + cst
   tmp   = log(tmp) + 1.d0 / tmp - 1.d0
   tmp   = 1.d0 / tmp
   rs    = CalRdisk(Vctmp, dens) / cst ! kpc/h
   Mdm   = Mhalo * (1.d0 - param%bar_rat)
   tff_s = delt * param%th_yr * 3.15d+7 ! [hubble time] -> [s]
   A     = 4.d-25 * Cube(param%h) * Square(tff_s) * Mdm / Cube(rs)
   B     = A * tmp
   y     = 1.d0 + x
   FFR   = cube(x) - B * (log(y) + 1.d0 / y - 1.d0)
! FFR(0) = 0 and only once FFR(x) becomes 0 other than x = 0
END FUNCTION FFR
!!$ ========================================================================
!!$ return DM fraction for galaxy mergers
DOUBLE PRECISION FUNCTION FracDM(id1)
   use global_var
   INTEGER, INTENT(IN):: id1
   FracDM = 1.d0 / param%bar_rat &
            * (gal(id1)%Morg/param%Mh0) ** param%alpha_halo
   IF(gal(id1)%flag_c == 0) FracDM = 0.d0
END FUNCTION FracDM
!!$ ========================================================================
