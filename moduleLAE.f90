module LAErelated
! --- Schaerer (2003) related PARAMETERs ---
  TYPE parametersLAE
     INTEGER :: nfile = 16, base = 500 ! for output files
     INTEGER :: nloop = 10 ! # of loop for elapsed time since the onset of starburst
     INTEGER :: type_b = 2 ! burst type (1: instantaneous, 2: exponential)
     INTEGER :: exttype = 2 ! dust geometry (1:screen, 2:slab)
     DOUBLE PRECISION :: Rtburst = 10.d0  ! t_burst/t_dyn
     DOUBLE PRECISION :: fphaseLya = 4.d0
  END TYPE parametersLAE
  TYPE(parametersLAE) :: paramLAE
  TYPE constantsLAE
     DOUBLE PRECISION :: fLya = 1.10626d-11 ! fLya = ene_Lya * QH02NLya
!!$     DOUBLE PRECISION :: ene_Lya = 1.63406d-11
!!$                         ! the energy of a Lyman alpha photon [erg]
!!$                         !  = 10.2 [eV] * 1.602176e-12 [erg/eV]
!!$     DOUBLE PRECISION :: QH02NLya = 0.677d0
!!$                         ! conversion factor from Q(HI) [s^-1]
!!$                         !  to N(Lya) [s^-1] in the case B recomb.
!!$                         !  at  T = 10^4 [K]
     DOUBLE PRECISION :: Zsun, OptV, OptB, AV0, AB0, Ngas, Tdyn
     DOUBLE PRECISION :: convH2 ! conv. factor from [erg/s] to [erg/s/h^2]
  END TYPE constantsLAE
  TYPE(constantsLAE) :: constLAE
  TYPE ThresholdForLAE
     DOUBLE PRECISION :: time   = 5.d+7  ! [yr]
     DOUBLE PRECISION :: LLya   = 3.d+41 ! [erg/s/h^2]
     DOUBLE PRECISION :: Mstar  = 0.d0   ! [Msun]
     DOUBLE PRECISION :: magUV  = -15.d0 ! [ABmag]
  END TYPE ThresholdForLAE
  TYPE(ThresholdForLAE) :: thLAE

  ! --- Schaerer (2003) Population Synthesis model
  INTEGER, PARAMETER :: NT_S03 = 1001
  INTEGER, PARAMETER :: NZ_S03(2) = (/9, 4/)
  INTEGER, PARAMETER :: NWAVE_S03 = 6 ! 1: LyC, 2:1500A, 3:2800A, 4:6563A, 5:4861A,
                                      ! 6:1216A
  TYPE Schaerer03
     DOUBLE PRECISION :: tstep, invtstep, base, base_iS03
     DOUBLE PRECISION, DIMENSION(:),     ALLOCATABLE :: time, metal
     DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE :: Lis, Lcs
  END TYPE Schaerer03
  TYPE(Schaerer03) :: S03base(2) ! raw data of S03
  TYPE Schaerer03Related
     DOUBLE PRECISION :: lum(NWAVE_S03, NT_S03)
  END TYPE Schaerer03Related
  TYPE(Schaerer03Related) :: S03(2) ! metallicity corrected data of S03

  TYPE PhysicalQuantitiesForAllGalaxies
     INTEGER :: num_b
     INTEGER :: flag_c ! 0:satellite, 1:central
     DOUBLE PRECISION :: Zc, Zc_pre
     DOUBLE PRECISION :: Mc_pre, Ms_pre
     DOUBLE PRECISION :: rb_pre, Vb_pre
     DOUBLE PRECISION :: sfr(2), beta(2) ! 1:starburst, 2:quiescent
     DOUBLE PRECISION :: taueff(2, 2)    ! 1:starburst, 2:quiescent
     DOUBLE PRECISION, ALLOCATABLE :: lumg_pre(:)
  END TYPE PhysicalQuantitiesForAllGalaxies
  TYPE(PhysicalQuantitiesForAllGalaxies), ALLOCATABLE :: allgal(:), &
       allgal_prev(:), allgal_next(:)

  TYPE PhysicalQuantitiesForLAE
     INTEGER :: id
     INTEGER :: flag_c ! 0:satellite, 1:central
     INTEGER :: b_or_q ! 1:starburst, 2:quiescent
     INTEGER :: sftype ! 1:quiescent, 2:single major merger, 3:multiple-merger
     INTEGER, DIMENSION(:), ALLOCATABLE :: flag_mag
     DOUBLE PRECISION :: NLyC, LLya
     DOUBLE PRECISION, DIMENSION(NWAVE_S03) :: LS03, LS03d, LS03pre
     DOUBLE PRECISION, DIMENSION(NWAVE_S03) :: magS03, magS03d
     DOUBLE PRECISION :: weight, tweight, weightAll
     DOUBLE PRECISION :: SFR, mSFR
     DOUBLE PRECISION :: beta, ab ! alp + beta
     DOUBLE PRECISION :: Mgas, Mgas_pre, Mgas_q, McZc_q
     DOUBLE PRECISION :: Zc, Zc_pre, reff, Ngas, NZ
     DOUBLE PRECISION :: Mstar, Ms_pre, Ms_q, Ms_b
     DOUBLE PRECISION :: Mhost, Mprog
     DOUBLE PRECISION :: tel, taust_yr, twind, Tmass
     DOUBLE PRECISION :: Vc  ! the circ. vel. of the host halo [km/s]
     DOUBLE PRECISION :: Vcc ! the circ. vel. of the central halo [km/s]
     DOUBLE PRECISION :: Vel ! the bulge vel. disp. or disk circ. vel. [km/s]
     DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: lumg, lumg_d, mag, mag_d
  END TYPE PhysicalQuantitiesForLAE
  TYPE(PhysicalQuantitiesForLAE) :: LAE

  DOUBLE PRECISION, DIMENSION(NWAVE_S03) :: xi_S03, corrS03

  ! --- For Luminosity Functions and Distributions
  INTEGER, PARAMETER :: NLF1 = 200, NLF2 = 10
  TYPE DistributionForLAE
     DOUBLE PRECISION, DIMENSION(2,3,NLF1) :: count, wcount
     DOUBLE PRECISION :: base, step, invstep
  END TYPE DistributionForLAE
  INTEGER, PARAMETER :: n_distLAE = 8
  TYPE(DistributionForLAE) :: LyaLF
  TYPE(DistributionForLAE) :: distLAE(n_distLAE)
                              ! 1:Ngas, 2:Z, 3:NZ, 4:Mstar, 5:Mhost,
                              ! 6:reff, 7:SFR, 8:taueff
  INTEGER, PARAMETER :: n_distLya = 4
  TYPE Lya_x_Distribution
     DOUBLE PRECISION :: count(2, NLF2, NLF1)
     DOUBLE PRECISION :: base, step, invstep
  END TYPE Lya_x_Distribution
  TYPE(Lya_x_Distribution) :: distLya(n_distLya)
                              ! 1:Lya-Ngas, 2:Lya-Z, 3:Lya-Mstar, 4:Lya-NZ
  TYPE LuminosityFunction
     ! LFs w/o and w/ dust extinction
     DOUBLE PRECISION, DIMENSION(3,3,NLF1) :: count, countd
  END TYPE LuminosityFunction
  TYPE(LuminosityFunction) :: LF_S03(NWAVE_S03)

  TYPE IntegersForLAE
     INTEGER :: LyaLF, contLF
     INTEGER :: distLya(n_distLya), dist(n_distLAE)
     INTEGER :: eachLAE, each
     INTEGER :: numLya, numUV
     INTEGER :: Makiya
  END TYPE IntegersForLAE
  TYPE(IntegersForLAE) :: iLAE
  CHARACTER(LEN=100), ALLOCATABLE :: fnameLAE(:)
  CHARACTER(LEN=100) :: fnameLAE_catalog_LyA, fnameLAE_catalog_UV
                        ! for parallel run

  TYPE DistributionFunction_LAE
     INTEGER :: Nbin, N1, N2, N3
     DOUBLE PRECISION :: step, invstep, base ! invstep = 1.d0 / step
     DOUBLE PRECISION, ALLOCATABLE :: bin(:) ! 1-Nbin
     DOUBLE PRECISION, ALLOCATABLE :: n(:,:,:,:) ! n(i,j,k,l)
                                      ! --- i:b_or_q = 1:starburst, 2:quiescent
                                      !     j:mor    = 1:E, 2:S0, 3:S, 4:all
                                      !     k:band   = 1-param%nwave
                                      !     l:bin    = 1-Nbin
  END type DistributionFunction_LAE
  TYPE(DistributionFunction_LAE), ALLOCATABLE :: LF_LAE(:), LFd_LAE(:)
                                                 ! luminosity functions

  DOUBLE PRECISION, ALLOCATABLE :: CSFRD(:,:,:) ! CSFRD(i,j,k)
                                   ! --- i:iz     = 1-nz_end
                                   !     j:b_or_q = 1:starburst, 2:quiescent
                                   !     k:mor    = 1:E, 2:S0, 3:S, 4:all
END module LAErelated
