
MODULE global_var
  implicit none

  INTEGER :: nloopz
  DOUBLE PRECISION :: tlife, tlife_yr, tlife_sec


!!$     ==================================
!!$     ===   Structure named "mrgt"   ===
!!$     ==================================
  TYPE nbody1
     INTEGER :: num_tot
     INTEGER :: num_galarray ! the total number of galaxy array
     ! the following arrays are 0-offset (0:num_tot-1)
     INTEGER, ALLOCATABLE :: f_des(:) ! descendant halo ID at z=0
     INTEGER, ALLOCATABLE :: n_des(:) ! descendant halo ID at next time slice
     INTEGER, ALLOCATABLE :: f_prg(:) ! progenitor halo ID at previous time slice
     INTEGER, ALLOCATABLE :: numb(:)  ! number of progenitor halos
     INTEGER, ALLOCATABLE :: hori(:)
     INTEGER(KIND=8), ALLOCATABLE :: mpi(:) ! marker particle index of a halo
                                            !  to identify the 3D position
     DOUBLE PRECISION, ALLOCATABLE :: mhalo(:) ! halo mass
     INTEGER, ALLOCATABLE :: hstart(:)
     INTEGER, ALLOCATABLE :: c_halo(:), c2_halo(:)
                             ! mrgt%c_halo(me), mrgt%c2_halo(me): halo IDs of the most
                             !  and the 2nd most massive progenitors of the halo "me"
     INTEGER, ALLOCATABLE :: num_g(:)
  END TYPE nbody1
  TYPE(nbody1) :: mrgt
!!$     ==============================


!!$     ===================================
!!$     ===    Structure named "mrgp"   ===
!!$     ===================================
  TYPE nbody2
     INTEGER :: num_step
     ! the following arrays are 1-offset (1:num_step)
     DOUBLE PRECISION, ALLOCATABLE :: time(:) ! <-- times(:)
     DOUBLE PRECISION, ALLOCATABLE :: zp1ar(:)
     INTEGER, ALLOCATABLE :: num_now(:) ! num_now(itime): total # of halos
                                        !  at time snapshot "itime"
     INTEGER, ALLOCATABLE :: st_halo(:) ! st_halo(itime): the 1st halo ID
                                        !  at time snapshot "itime"
     INTEGER, ALLOCATABLE :: num_tot_gal(:)
  END TYPE nbody2
  TYPE(nbody2) :: mrgp
!!$     ===================================


!!$     ================================================================
!!$     ===   Structure for the model parameters named "parameters"  ===
!!$     ================================================================
  TYPE parameters
     INTEGER :: n_forest ! the number of N-body data files (forest of halo merger trees)
     LOGICAL :: run_all_ssp ! .true.: run all SSPs,
                            ! .false.: run a certain SSP
     LOGICAL :: off_collision ! .true.  : random collision does not occur
                              ! .false. : random collision occurs
     LOGICAL :: equal_mass_merger 
                ! .true.  : equal mass merger used in subroutine star
                ! .false. : non-equal mass merger
     LOGICAL :: dyn_resp_bulge ! logical valuable for dyn. resp. of bulge
     LOGICAL :: dyn_resp_disk  ! logical valuable for dyn. resp. of disk
     LOGICAL :: dyn_resp_halo  ! logical valuable for dyn. resp. of halo
     LOGICAL :: ML  ! logical valuable for calculation of ML as a function of M_B
     LOGICAL :: Mbh ! logical valuable for calculation of Mbh as a function of Mbulge
     LOGICAL :: LAE ! logical valuable for calculation of LAE
     LOGICAL :: SFH ! logical valuable for calculation of SFH
     LOGICAL :: SN  ! logical valuable for calculation of SN
     LOGICAL :: SB  ! logival valuable for calculation of LF with a given
                    !   surface brightness limit (added by MARK, 2015/Apr/11)
     LOGICAL :: DI  ! logical valuable for calculation disk instability
                    ! (added by Shirakata, 2016/Feb/25)
     INTEGER :: ML_AGN ! mass/lumi. of AGNs 1: L = eMc^2 2: Watarai+00
     DOUBLE PRECISION :: Edd_max ! for  ML_AGN = 3
     INTEGER :: BolCor ! Bolometric correction obtained by 1: Marconi+04 
     INTEGER :: ObsFrac ! observable fraction of AGNs  (Added by Shirakata 18/Mar/05)
                        ! 1: Shirakat+18 2:Hopkins+07

     LOGICAL :: delay ! logical valuable for delayed accretion model of AGN
                      ! (added by Shirakata, 2016/Aug/20)
     CHARACTER(LEN=30) :: file_nbody

     INTEGER :: run_type ! 1:run at a redshift
                         ! 2:run at all redshift
                         ! 3:MCMC
     INTEGER :: traceIDs ! 0: off 1: on
     INTEGER :: loopz, izout
     INTEGER :: idum
     INTEGER :: iy, iv(32), iset ! added by Makiya, for ran1 and gasdev (2014/Mar/18)

     INTEGER :: type_ssp ! 1:KA97, 2:PEGASE, 3:BC03
     INTEGER :: type_mag ! 1:Vega, 2:AB (added by MARK, 2014/Feb/02)
     INTEGER :: exttype  ! 1:MW, 2:LMC, 3:SMC (Pei92), 4: Calzetti-law (Calzetti00)
     INTEGER :: SFmodel  ! 1:CSF, 2:DSF 3: Makiya14(added by Makiya, 2015/Apr/04)

     INTEGER :: nwave ! # of filters used in calculation
     INTEGER, ALLOCATABLE :: iwave(:) ! iwave(1:param%nwave)
                                      ! this array provides the correspondence
                                      !  to integer of filter number for lumi
     INTEGER :: nwp1, tnw ! nwp1 = nwave + 1, tnw = 2 * nwave
     CHARACTER(LEN=500) :: line_filters 
                      ! characters which represents the filter names
                      !  used in calculation
     LOGICAL :: wdim_band ! .true. = bands w/ dimension are used
                          ! (added by MARK on 2018/Mar/02)

     ! corresponding integers to each band filter
     !  --- these are defined in SUBROUTINE "SubstFilterNumb" in main_nugc.f
     !       and within 1-param%nwave
     !  --- e.g., Bband: mag(iBband), gal(*):lumg(iBband)
     INTEGER :: iBband, iUband, iVband, iRcband, iRJband
     INTEGER :: iIcband, iIJband, izband, iJband, iHband
     INTEGER :: iKband, iKpband, iLband
     INTEGER :: iSuprime_Bband, iSuprime_gband, iSuprime_Vband, iSuprime_rband, &
          iSuprime_RLband, iSuprime_ipband, iSuprime_Iband, iSuprime_zpband, &
          iSuprime_zp_redband
     INTEGER :: iobs5um
     INTEGER :: i2MASS_Jband, i2MASS_Hband, i2MASS_Ksband
     INTEGER :: iACS_F775Wband, iACS_F850LPband
     INTEGER :: iCISCO_zband, iCISCO_Jband, iCISCO_Hband,&
          iCISCO_Kband, iCISCO_Kpband
     INTEGER :: iGALEX_FUVband, iGALEX_NUVband
     INTEGER :: iGOODS_RBband, iGOODS_RGband, iGOODS_RSband
     INTEGER :: iHST_F300wband, iHST_F450wband, iHST_F555wband,&
          iHST_F606wband, iHST_F702wband, iHST_F814wband
     INTEGER :: iSDSS_upband, iSDSS_gpband, iSDSS_rpband, &
          iSDSS_ipband, iSDSS_zpband
     INTEGER :: iNLyC, iL1216, iL1400, iL1500, iL1600, iL1700, iL2800, &
          iL4861, iL6563
     INTEGER :: iIRACch1band, iIRACch2band, iIRACch3band, iIRACch4band
     INTEGER :: iWFCAM_zband, iWFCAM_Yband, iWFCAM_Jband, iWFCAM_Hband, iWFCAM_Kband
     ! --- following filters are added since 2014/Oct/12
     INTEGER :: iVIRCAM_Yband, iVIRCAM_Jband, iVIRCAM_Hband, iVIRCAM_Kband
     INTEGER :: iCIBER_Iband, iCIBER_Hband
     INTEGER :: iNEWFIRM_J1band, iNEWFIRM_J2band, iNEWFIRM_J3band, iNEWFIRM_H1band, &
          iNEWFIRM_H2band, iNEWFIRM_Ksband
     INTEGER :: iAKARI_N2band, iAKARI_N3band, iAKARI_N4band, iAKARI_S7band, &
          iAKARI_S9Wband, iAKARI_S11band
     INTEGER :: iWISH_0band, iWISH_1band, iWISH_2band, iWISH_3band, &
          iWISH_4band, iWISH_5band
     INTEGER :: iACS_F435Wband, iACS_F475Wband, iACS_F606Wband, iACS_F814Wband, &
          iNICMOS_F160Wband, iWFPC2_F300Wband, iWFPC2_F450Wband, iWFPC2_F555Wband, &
          iWFC3_F125Wband, iWFC3_F140Wband, iWFC3_F160Wband
     INTEGER :: iHSC_gband, iHSC_rband, iHSC_iband, iHSC_zband, iHSC_Yband

     INTEGER :: iBband_r, iUband_r, iVband_r, iRcband_r, iRJband_r
     INTEGER :: iIcband_r, iIJband_r, izband_r, iJband_r, iHband_r
     INTEGER :: iKband_r, iKpband_r, iLband_r
     INTEGER :: iSuprime_Bband_r, iSuprime_gband_r, iSuprime_Vband_r, &
          iSuprime_rband_r, iSuprime_RLband_r, iSuprime_ipband_r, iSuprime_Iband_r, &
          iSuprime_zpband_r, iSuprime_zp_redband_r
     INTEGER :: iobs5um_r
     INTEGER :: i2MASS_Jband_r, i2MASS_Hband_r, i2MASS_Ksband_r
     INTEGER :: iACS_F775Wband_r, iACS_F850LPband_r
     INTEGER :: iCISCO_zband_r, iCISCO_Jband_r, iCISCO_Hband_r,&
          iCISCO_Kband_r, iCISCO_Kpband_r
     INTEGER :: iGALEX_FUVband_r, iGALEX_NUVband_r
     INTEGER :: iGOODS_RBband_r, iGOODS_RGband_r, iGOODS_RSband_r
     INTEGER :: iHST_F300wband_r, iHST_F450wband_r, iHST_F555wband_r,&
          iHST_F606wband_r, iHST_F702wband_r, iHST_F814wband_r
     INTEGER :: iSDSS_upband_r, iSDSS_gpband_r, iSDSS_rpband_r, &
          iSDSS_ipband_r, iSDSS_zpband_r
     INTEGER :: iNLyC_r, iL1216_r, iL1400_r, iL1500_r, iL1600_r, &
          iL1700_r, iL2800_r, iL4861_r, iL6563_r
     INTEGER :: iIRACch1band_r, iIRACch2band_r, iIRACch3band_r, iIRACch4band_r
     INTEGER :: iWFCAM_zband_r, iWFCAM_Yband_r, iWFCAM_Jband_r, iWFCAM_Hband_r, &
          iWFCAM_Kband_r
     ! --- following filters are added since 2014/Oct/12
     INTEGER :: iVIRCAM_Yband_r, iVIRCAM_Jband_r, iVIRCAM_Hband_r, iVIRCAM_Kband_r
     INTEGER :: iCIBER_Iband_r, iCIBER_Hband_r
     INTEGER :: iNEWFIRM_J1band_r, iNEWFIRM_J2band_r, iNEWFIRM_J3band_r, &
          iNEWFIRM_H1band_r, iNEWFIRM_H2band_r, iNEWFIRM_Ksband_r
     INTEGER :: iAKARI_N2band_r, iAKARI_N3band_r, iAKARI_N4band_r, iAKARI_S7band_r, &
          iAKARI_S9Wband_r, iAKARI_S11band_r
     INTEGER :: iWISH_0band_r, iWISH_1band_r, iWISH_2band_r, iWISH_3band_r, &
          iWISH_4band_r, iWISH_5band_r
     INTEGER :: iACS_F435Wband_r, iACS_F475Wband_r, iACS_F606Wband_r, iACS_F814Wband_r, &
          iNICMOS_F160Wband_r, iWFPC2_F300Wband_r, iWFPC2_F450Wband_r, &
          iWFPC2_F555Wband_r, iWFC3_F125Wband_r, iWFC3_F140Wband_r, iWFC3_F160Wband_r
     INTEGER :: iHSC_gband_r, iHSC_rband_r, iHSC_iband_r, iHSC_zband_r, iHSC_Yband_r

     ! --- cosmology related parameters
     DOUBLE PRECISION :: OMEGA0, OMEGA ! OMEGA0=Omega_M, OMEGA=Omega_b
     DOUBLE PRECISION :: OMEGA_L, OMEGA_rat ! OMEGA_L=Omega_L, OMEGA_rat=Omega_L/Omega_M
     DOUBLE PRECISION :: bar_rat ! = OMEGA / OMEGA0
     DOUBLE PRECISION :: h ! <-- hubble
     DOUBLE PRECISION :: munit, log10munit ! unit for mass = 10^14 Msun
     DOUBLE PRECISION :: th, th_yr ! Hubble time

     DOUBLE PRECISION :: zsp1
     DOUBLE PRECISION :: zsp1_input ! z+1 in input file (added by MARK, 2016/Jun/24)

     ! --- model parameters
     DOUBLE PRECISION :: alpst
     DOUBLE PRECISION :: tau0st
     DOUBLE PRECISION :: eps_SF ! for DSF mode (added by Makiya, 2015/Apr/04)
     DOUBLE PRECISION :: emax, emin ! for Makiya14 mode (added by Shirakata, 2017/11/02)
     DOUBLE PRECISION :: Zch, taud_th ! for Makiya14 mode (added by Shirakata, 2017/11/02)
     DOUBLE PRECISION :: fmerge
     DOUBLE PRECISION :: fmajor
     DOUBLE PRECISION :: Krem ! for remained energy fraction in the disk 
                              ! with bulge growth (Shirakata, 2017/Oct/20) 
     DOUBLE PRECISION :: freheat
     DOUBLE PRECISION :: fb
     DOUBLE PRECISION :: Vcut
     DOUBLE PRECISION :: Vlow
     DOUBLE PRECISION :: tauV0
     DOUBLE PRECISION :: fdm
     DOUBLE PRECISION :: fdiss   ! dissipation fraction during major merger
     DOUBLE PRECISION :: alpha_tau ! tau propto (1+z)^(-alpha_tau)
     DOUBLE PRECISION :: alpha_rad ! r_disk propto R_vir (1+z)^alpha_rad
     DOUBLE PRECISION :: alpha_burst ! tburst/taudyn
     DOUBLE PRECISION :: Vhot(2)
     DOUBLE PRECISION :: alphot(2)
     DOUBLE PRECISION :: Vst    ! added by Shirakata (2016/Oct/30)
     DOUBLE PRECISION :: alp_ret ! returned gas mass fraction from reservoir to
                                 ! hot halo (added by Shirakata 2016/Oct/31) 
     DOUBLE PRECISION :: f_ngal ! the number fraction of galaxies to DM halos
                                !  (added by MARK, 2016/Jun/29)

     ! --- SMBH related parameters
     INTEGER :: nwaveq
     DOUBLE PRECISION :: fbh(2)
     DOUBLE PRECISION :: Mbhseed
     DOUBLE PRECISION :: eps_agn
     DOUBLE PRECISION :: tacc_0, tacc_1
     DOUBLE PRECISION :: acc
     DOUBLE PRECISION :: alpha_ad, gamma_ad, tad_0

     ! --- AGN feedback related parameters
     INTEGER :: AGNFB_key ! 1:Vcut, 2:Croton, 3:Bower
     DOUBLE PRECISION :: kappa_croton, eta_croton
     DOUBLE PRECISION :: alp_bower, eps_bower
     DOUBLE PRECISION :: Lcool0, Ledd0

     ! --- Cooling Function added by Shirakata (2018/Aug/12)
     INTEGER :: CoolFN ! 1: Satherland and Dopita, 2: Okamoto cloudy 
     ! --- Disk instabiltiy
     DOUBLE PRECISION :: em
     DOUBLE PRECISION :: fdi

     ! --- Bulge size and velocity dispersion
     DOUBLE PRECISION :: Mh0, alpha_halo
     ! --- UV feedback related parameters
     INTEGER :: UVfb ! 1:Vlow, 2:Okamoto+08
     INTEGER :: iO08_p, iO08_n0, iO08_n1, iO08_n2, iO08_n3, iO08_n4
                ! integers for Okamoto+08 UV feedback
     DOUBLE PRECISION :: zreion

     DOUBLE PRECISION :: mu = 0.593d0 ! added by MARK (2015/Apr/07)
     DOUBLE PRECISION :: SBlimit ! surface brightness limit for calculation of LF
                                 ! only valid if param%SB = .true.
                                 ! added by MARK (2015/Apr/11)

     ! --- quantities used in subroutine cool1 (added by MARK, 2015/Apr/13)
     DOUBLE PRECISION :: Tvir_min = 1.d+4 ! [K]
     DOUBLE PRECISION :: Tvir_max = 10.d0**8.5d0 ! [K]
  END TYPE parameters
  TYPE(parameters) :: param
!!$     ================================================================


!!$     ================================================
!!$     ===   SSP related parameters and structure   ===
!!$     ================================================
  INTEGER, PARAMETER :: N_SSP = 3 ! 1:KA97, 2:PEGASE, 3:BC03
  INTEGER, PARAMETER :: NSFR_BASE(N_SSP)  = (/54, 516, 221/)
  INTEGER, PARAMETER :: NCHEM_BASE(N_SSP) = (/9, 7, 6/)
  INTEGER, PARAMETER :: NWAVE_ALL  = 107 ! # of filters of ssp
  INTEGER, PARAMETER :: TNWAVE_ALL = 214 ! = 2 * NWAVE_ALL
  INTEGER, PARAMETER :: N_SSP_TOT = 27 ! total number of ssp
                                       ! --- 1 KA97 + 2 BC03 (Salpeter + Chabrier IMFs)
                                       !     + 24 PEGASE (6 IMFs x 4 mass ranges)
  TYPE struct_ssp
     CHARACTER(LEN=20) :: bandname(TNWAVE_ALL)
     CHARACTER(LEN=100) :: filename(2) ! name of sspfile (1:burst, 2:quiescent)

     DOUBLE PRECISION, ALLOCATABLE :: time(:) ! time slices of ssp
     DOUBLE PRECISION, ALLOCATABLE :: chem(:) ! metallicity slices of ssp
     DOUBLE PRECISION, ALLOCATABLE :: lumi(:,:,:,:) 
     ! luminosity density [ergs/s/cm^2/Hz] at the distance of 10 pc at each
     !  time and metallicity slices separately given for burst and quiescent
     ! lumi(k,i,t,z) --- k: burst (1) or quiescent (2)
     !                   i: band filter (i=1-TNWAVE_ALL)
     !                   t: time slices (t=1-NSFR)
     !                   z: metallicity slices (z=1-NCHEM)

     DOUBLE PRECISION :: lam_c(TNWAVE_ALL)
                         ! central wavelength of each band [um]
     ! alp(i), R(i), p(i), y(i) --- i = 1:burst, 2:quies
     DOUBLE PRECISION :: alp(2) ! locked-up mass fraction
     DOUBLE PRECISION :: R(2)   ! returned mass fraction
     DOUBLE PRECISION :: y(2)   ! yeild
     DOUBLE PRECISION :: p(2)   ! p(:) = alp(:) * y
  END TYPE struct_ssp
  TYPE(struct_ssp) :: ssp
  INTEGER :: NSFR  ! # of time slices of ssp
  INTEGER :: NCHEM ! # of metallicity slices of ssp

  ! --- for the bands w/ dimension (added by MARK on 2018/Mar/02)
  INTEGER, PARAMETER :: NWAVE_wdim = 9 ! 1:LyC, 2:1216A, 3:1400A, 4:1500A, 5:1600A,
                                       ! 6:1700A, 7:2800A, 8:4861A, 9:6563A
!!$     ================================================


!!$     ==============================================
!!$     ===   Structure for galaxy named "galaxy"  ===
!!$     ==============================================
  TYPE galaxy
     INTEGER :: flag_burst  ! flag for starburst induced by mergers or DI
     INTEGER :: flag_di     ! flag for disk instability
     INTEGER :: flag_merger ! flag for merger
     INTEGER :: flag_seed   ! flag for putting a seed BH
     INTEGER :: flag_ccut   ! = CoolingCutoff_key
     INTEGER :: flag_dead
     INTEGER :: flag_cfalse ! satellite is larger than central
     INTEGER :: IDhost  ! host halo ID for mrgt: mrgt%mhalo(gal(id)%IDhost) provides
                        !  the mass of host halo
     INTEGER :: IDprog  ! progenitor halo ID (added by MARK on 2013/Sep/27)
     INTEGER :: id_cgal ! ID of central galaxy
     INTEGER :: flag_c
     INTEGER :: hori   ! <-- horiG
     INTEGER(KIND=8) :: mpi ! marker particle index to identify the 3D position
     INTEGER :: hstart ! <-- hstartG
     INTEGER :: hfinal ! <-- hfinalG
     INTEGER :: n_merge ! number of merging galaxies at this timestep 

     DOUBLE PRECISION :: z_col   ! first collapsed redshift
     DOUBLE PRECISION :: z_form1 ! halo formation redshift
     DOUBLE PRECISION :: Mhalo   ! halo mass in which the gal. is the cen. gal.
                                 ! --- cen.: gal(id)%Mhalo = mrgt%mhalo(gal(id)%host)
                                 ! --- sat.: gal(id)%Mhalo provides subhalo mass
     DOUBLE PRECISION :: Mhotout ! expelled gas from halo by superwind
     DOUBLE PRECISION :: Vc      ! halo circular velocity
     DOUBLE PRECISION :: Mreheat 
     DOUBLE PRECISION :: Telapse 
     DOUBLE PRECISION :: cst
     DOUBLE PRECISION :: Mratio  ! = Mhot / Mhalo
     DOUBLE PRECISION :: Morg    !
     DOUBLE PRECISION :: Mhotorg !
     DOUBLE PRECISION :: MZhorg  !
     DOUBLE PRECISION :: Mstard, Mstarb ! disk and bulge stellar masses
     DOUBLE PRECISION :: MZd, MZb       ! heavy element masses in disk and bulge
     DOUBLE PRECISION :: Vdisk, Vbulge  ! rotation velocity of disk and velocity
                                        !  dispersion of bulge
     DOUBLE PRECISION :: Vmax           ! for disk instability
     DOUBLE PRECISION :: rdisk, rbulge  ! effective sizes of disk and bulge
     DOUBLE PRECISION :: Mhot           ! gas mass in hot gas phase
     DOUBLE PRECISION :: Mcoold, Mcoolb ! gas mass in cold gas phase
                                        ! disk, bulge, respectively
                                        ! (Shirakata, 2016/Jul/20)
     DOUBLE PRECISION :: MZcd, MZcb, MZh ! heavy element mass in hot gas phase
     DOUBLE PRECISION :: growthrate
     DOUBLE PRECISION :: diskmas
     DOUBLE PRECISION :: density
     DOUBLE PRECISION :: Mtd, Mtb
     DOUBLE PRECISION :: clps
     DOUBLE PRECISION :: Vcent
     DOUBLE PRECISION :: MZc_rem
     DOUBLE PRECISION :: SFR
     DOUBLE PRECISION :: mSFR ! mean SFR during the past const%t0yr [yr]
                              ! (added by MARK on 2014/Nov/04)
     DOUBLE PRECISION :: taust ! taustar (added by MARK on 2017/Mar/16)
     DOUBLE PRECISION :: dMstar_burst ! stellar mass newly formed in starburst


     DOUBLE PRECISION :: mzg, mtg

     DOUBLE PRECISION, ALLOCATABLE :: lumg(:) ! 1:param%tnw
                                      ! --- lumg(1:param%nwave):        bulge
                                      !     lumg(param%nwp1:param%tnw): disk
     DOUBLE PRECISION :: Zmassb, Zmassd ! mass-weighted stellar metallicity
                                        ! (added by MARK on 2014/Nov/02)
     DOUBLE PRECISION :: BT             ! BT ratio (added by Shirakata on 2016/03/24)

     DOUBLE PRECISION :: taumrg       ! dynamical friction timescale (for satellites)
                                      ! renewed when most massive halo has been
                                      !  changed (added by Shirakata on 2017/02/27)

     DOUBLE PRECISION :: M2M1_av, M2M1_max ! merger mass ratio
                                           !(single-averaged and maximum values
                                           ! added by shirakata on 2017/08/24)
     ! --- related to BH ---
     DOUBLE PRECISION :: Mbh
     DOUBLE PRECISION :: Lpeak
     DOUBLE PRECISION, ALLOCATABLE :: lumq(:) ! QSO luminosity
     DOUBLE PRECISION, ALLOCATABLE :: agn(:) ! AGN properties  (Shirakata)
     ! --- begining of star formation redshift (by MN)
!!$         DOUBLE PRECISION :: z_sf
     INTEGER :: z_sf

     ! --- LAE related (by MARK) ---
     DOUBLE PRECISION :: Tmass, Tlum_b, Tlum_d
                         ! stellar mass and luminosity-weighted ages
     DOUBLE PRECISION :: beta ! SN feedback intensity

     ! --- SFH related (by MARK) ---
     DOUBLE PRECISION, ALLOCATABLE :: SFH(:,:) ! SFH(itime,iZ)

     ! --- ID trace (by Shirakata) --- 
     DOUBLE PRECISION, ALLOCATABLE :: mem(:)
!!$     DOUBLE PRECISION, ALLOCATABLE :: SFH0(:,:) ! SFH0(itime,iZ)
  END TYPE galaxy
  TYPE(galaxy), ALLOCATABLE :: gal(:), gal_prev(:), gal_next(:)
!!$     ==============================================

!!$     ==============================================
!!$     ===     For tracing some galaxies' IDs     ===
!!$     ==============================================
! ------ Added by Shirakata (2018/02/07)
TYPE targets
   INTEGER :: iforest, fdes, hstart 
END TYPE targets
TYPE(targets), ALLOCATABLE :: targ(:)

INTEGER :: num_tar
!!$     ==============================================

  INTEGER :: num_total ! tot. # of gals w/ non-zero baryon at a certain timestep
  INTEGER :: num_cum ! cumulative # of gals in halos from the halo "ihalo" = 1
                     !   at the time-step "ihalo"
  DOUBLE PRECISION :: taust, Mc0, Ms0, MZc0, Zc0, Zh0, Mbh0, Mrem
  DOUBLE PRECISION :: dMcold, dMstar, dMhot, Zc, dMZh, dMbh, dMbh_norm
  DOUBLE PRECISION :: beta
  DOUBLE PRECISION :: Tmass0, Tlum_b0, Tlum_d0 ! added by MARK
  DOUBLE PRECISION :: Zmassb0, Zmassd0 ! added by MARK (2014/Nov/02)
  DOUBLE PRECISION :: SFR0, mSFR0 ! added by MARK (2014/Nov/07)
  DOUBLE PRECISION, ALLOCATABLE :: mag(:), mag_d(:) ! 1:param%nwave
  DOUBLE PRECISION, ALLOCATABLE :: xi(:) ! 1:param%nwp1
  DOUBLE PRECISION, ALLOCATABLE :: mag_q(:) ! 1:param%nwaveq


!!$     ===================================================
!!$     ===   cooling function parameter and variables  ===
!!$     ===================================================
  INTEGER, PARAMETER :: NCOOLFN = 100
  DOUBLE PRECISION :: lambda(NCOOLFN), lambdas(NCOOLFN)
!!$     ===================================


!!$     ========================================================
!!$     ===   Structure for useful constants named "const"   ===
!!$     ========================================================
  TYPE constants
     DOUBLE PRECISION :: PI, PI2 ! PI2 = PI^2
     DOUBLE PRECISION :: Zsun ! solar metallicity
     DOUBLE PRECISION :: corr ! = -5logh
     DOUBLE PRECISION :: corr_wdim(NWAVE_wdim) ! added by MARK (2018/Mar/02)
     DOUBLE PRECISION :: nub  ! = frequency at the highest sensitivity in B-band
     DOUBLE PRECISION :: B2UV ! mag_UV = mag_B + B2UV for AGNs 
                              ! added by Shirakata (2018/Mar/05)
     DOUBLE PRECISION :: ML ! for Hydrogen mass to L_B
     DOUBLE PRECISION :: delc0 ! for the function of "delc"
     DOUBLE PRECISION :: V1, V2  ! for the function of "CircVel"
     DOUBLE PRECISION :: Rb ! = 2.175d+8*param%h [kpc/h] for the function of "CalRbulge"
     DOUBLE PRECISION :: fracH ! Hydrogen mass fraction in cold gas mass
     DOUBLE PRECISION :: tout, toutGyr ! the age of the universe at z = param%zsp1-1
                                       !  in the units of hubble time and Gyr
     DOUBLE PRECISION :: t0yr = 10.d+6, t0
                         ! used to calculate the mean SFR during the past t0yr
                         ! t0 = t0yr / param%th_yr
     ! --- for baryon fraction in Okamoto+08 UV feedback mode
     DOUBLE PRECISION :: alpha_UV, UV1, UV2 ! UV1 = 2^{alpha_UV/3} - 1, UV2 = -3/alpha_UV

     ! --- add by Makiya, for AGN feedback calc.
     DOUBLE PRECISION :: c = 2.99792458d+8 ! speed of light [m/s]
     DOUBLE PRECISION :: Msolar = 1.989d+30 ! solar mass [kg]
     DOUBLE PRECISION :: yr2sec = 3.1556926d+7 ! yr [sec]
     ! --- the followings are added by MARK (2015/Apr/11)
     DOUBLE PRECISION :: mp = 1.67262158d0 / 1.9884d0 ! proton mass [10^-57 Msun]
     DOUBLE PRECISION :: kB = 1.38d0 ! Boltzman const. [10^-16 erg/K]
     DOUBLE PRECISION :: t_dyn0, dMbh0 ! used only in the Bower AGN feedback mode
     ! --- the followings are added by Shirakata (2015/May/15)
     DOUBLE PRECISION :: G = 4.905d-15  ! gravitational const. [pc^3 Msun^-1 yr^-2]
     DOUBLE PRECISION :: T_si = 1.d0    ! silicate grain sublimation temperature [1500K]
  END TYPE constants
  TYPE(constants) :: const
!!$     ========================================================

  ! --- CSFR, by Makiya
  DOUBLE PRECISION, ALLOCATABLE :: CSFR(:,:), CSFR_all(:,:)

  ! --- SMF at all-z
  DOUBLE PRECISION, ALLOCATABLE :: SMF_z(:,:), SMF_z_all(:,:)
  DOUBLE PRECISION :: SMF_base
  INTEGER :: SMF_bin

  ! --- for parallel run, by Makiya
  INTEGER :: nnode, inode, istat

  ! --- LF and MF related constants
  INTEGER, PARAMETER :: NbinLF = 200
  DOUBLE PRECISION, PARAMETER :: stepLF = 0.5d0, baseLF = -40.25d0
  ! for mass functions
  INTEGER, PARAMETER :: NtypeMF = 11, NbinMF = 100
  DOUBLE PRECISION, PARAMETER :: stepMF = 0.2d0
  ! for M_HI/L_B
  INTEGER, PARAMETER :: NbinML = 50
  INTEGER :: NgalMax = 4000000, NgalMax2 = 20000000
  DOUBLE PRECISION, PARAMETER :: stepML = 1.d0, baseML = -30.d0
  ! for Mbh - Mbulge, by Makiya
  DOUBLE PRECISION, PARAMETER :: stepMbhMbulge = 0.2d0, baseMbhMbulge = 3.d0
  ! for Scaling relation of disk by Shirakata
  DOUBLE PRECISION, PARAMETER :: stepDiskScale = 0.2d0, baseDiskScale = 1.d0
  ! for Scaling relation of bulge by Shirakata
  DOUBLE PRECISION, PARAMETER :: stepSphScale = 0.5d0, baseSphScale = -30.d0

  ! --- for AGN feedback, by Makiya (2015/Apr/04)
  INTEGER :: CoolingCutoff_key

  ! --- for conversion of surface brightness from total to Petrosian,
  !       by Makiya (2015/Apr/04): modified by MARK (2014/Apr/11)
  TYPE SurfaceBrightness
     INTEGER :: N_Rratio = 500, N_BT = 500, iband = 3
     INTEGER :: N_Rratiop1, N_BTp1
     DOUBLE PRECISION, ALLOCATABLE :: Pet2TotRatio(:,:,:)
                                      ! Pet2TotRatio(k, 1:N_Rratio+1, 1:N_BT+1)
                                      ! k=1: -2.5log10(Lpet/Ltot)
                                      ! k=2: R50(Lpet)/Re(bulge)
                                      ! k=3: R50(Lpet)/Re(disk)
     DOUBLE PRECISION :: Rratio_base = -3.d0
     DOUBLE PRECISION :: Rratio_step = 6.d0 / 500.d0
                         ! Rratio is provided from 10^-3 to 10^3 divided equally
                         !   spaced in log-scale with N_Ratio
                         ! 6.d0 = log10(10^3) - log10(10^-3)
     DOUBLE PRECISION :: BT_base = 0.d0, BT_step = 2.d-3
                         ! BT ratio is provided from 0 to 1 divided equally
                         !   spaced in linear-scale with N_BT
     INTEGER :: iRratio, iBT
     DOUBLE PRECISION :: BT, mag_Pet, r_Pet
     DOUBLE PRECISION :: mu0 ! constant for surface brightness calculation
  END type SurfaceBrightness
  TYPE(SurfaceBrightness) :: SB

  !--- for parallel run, added by Makiya 2015/09/27
  TYPE DistributionFunction1
     INTEGER :: Nbin, N1, N2
     DOUBLE PRECISION :: step, invstep, base ! invstep = 1.d0 / step
     DOUBLE PRECISION, ALLOCATABLE :: bin(:)
     DOUBLE PRECISION, ALLOCATABLE :: n(:,:,:), n_brst(:,:,:)
                                      ! n(i,j,k) --- i:band,
                                      !              j:mor(1:E,2:S0,3:S,4:all),
                                      !              k:mag
  END type DistributionFunction1
  TYPE(DistributionFunction1), ALLOCATABLE :: lf(:), lf_d(:), lf_q(:) ! LFs
  TYPE(DistributionFunction1), ALLOCATABLE :: mf(:) ! MFs
  INTEGER :: nforest_this, iforest_start
  INTEGER :: N1, N2, N3, N1N2Nbin
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: lf_n_all
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: mf_n_all
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: lf_q_all
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:)     :: ML_bin
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: ML_n_all, ML_x_all
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: ML_med
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:)     :: MbhMbulge_bin
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: MbhMbulge_n_all, &
       MbhMbulge_xn_all, MbhMbulge_xxn_all
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:)     :: DiskScale_bin
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: DiskScale_n_all, &
       DiskScale_xn_all, DiskScale_xxn_all
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:)     :: SphScale_bin
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: SphScale_n_all, &
       SphScale_xn_all, SphScale_xxn_all

  ! --- for parallel run, by Makiya 2015/09/27
  CHARACTER(LEN=100) :: file_catalog, file_catalog_q, file_catalog_t, catalog_all
  ! CHARACTER(LEN=100) :: command
  CHARACTER(LEN=5000) :: command ! modified by MARK 2015/12/17
  INTEGER, PARAMETER :: NFILE = 8

  ! --- for concentration parameter, added bu Makiya 2015/10/07
  TYPE ConcentrationParameter
     DOUBLE PRECISION, ALLOCATABLE :: CST(:,:)
     INTEGER :: Nz = 50, NM = 100
     DOUBLE PRECISION :: zp1_base = 0.d0, zp1_step = log10(25.)/50.d0
     DOUBLE PRECISION :: M_base = 8.d0, M_step = 7.d0/100.d0
     DOUBLE PRECISION :: a, b
  END TYPE ConcentrationParameter
  TYPE(ConcentrationParameter) :: CP
!!$     ================================================================
!!$     ===   Structure for the MCMC parameters named "mcmcparams"  ===
!!$     ================================================================
  TYPE mcmcparams
     INTEGER :: iMCMC_max
     ! --- Star formation
     LOGICAL :: alpst, Vst, tau0st, eps_SF, emax, emin, Zch, taud_th
     DOUBLE PRECISION :: alpst_min, Vst_min, tau0st_min, eps_SF_min, &
                         emax_min, emin_min, Zch_min, taud_th_min
     DOUBLE PRECISION :: alpst_max, Vst_max, tau0st_max, eps_SF_max, &
                         emax_max, emin_max, Zch_max, taud_th_max
     ! --- SN feedback and returned gas  
     LOGICAL :: Vhot(2), alphot(2), alp_ret
     DOUBLE PRECISION ::Vhot_min(2), alphot_min(2), alp_ret_min
     DOUBLE PRECISION ::Vhot_max(2), alphot_max(2), alp_ret_max
     ! --- Mergers of galaxies 
     LOGICAL :: fmerge, fmajor, Krem, fdiss, Mh0, alpha_halo 
     DOUBLE PRECISION :: fmerge_min, fmajor_min, Krem_min, &
                         fdiss_min, Mh0_min, alpha_halo_min
     DOUBLE PRECISION :: fmerge_max, fmajor_max, Krem_max, &
                         fdiss_max, Mh0_max, alpha_halo_max
     ! --- AGN feedback 
     LOGICAL :: Vcut, kappa_croton, eta_croton, alp_bower, eps_bower
     DOUBLE PRECISION :: Vcut_min, kappa_croton_min, eta_croton_min, &
                         alp_bower_min, eps_bower_min
     DOUBLE PRECISION :: Vcut_max, kappa_croton_max, eta_croton_max, &
                         alp_bower_max, eps_bower_max
     ! --- Disc instabilities 
     LOGICAL :: em, fdi 
     DOUBLE PRECISION :: em_min, fdi_min
     DOUBLE PRECISION :: em_max, fdi_max
     ! --- SMBHs and AGNs
     LOGICAL :: fbh(2), Mbhseed, tacc_0, tacc_1, alpha_ad, gamma_ad, Edd_max 
     DOUBLE PRECISION :: fbh_min(2), Mbhseed_min, tacc_0_min, tacc_1_min, &
                         alpha_ad_min, gamma_ad_min, Edd_max_min
     DOUBLE PRECISION :: fbh_max(2), Mbhseed_max, tacc_0_max, tacc_1_max, &
                         alpha_ad_max, gamma_ad_max, Edd_max_max
  END TYPE mcmcparams
  TYPE(mcmcparams) :: mcmc
  INTEGER, PARAMETER :: NObsProp = 15 ! observational properties using MCMC 
END MODULE global_var
!===============================================================================
