module SNrelated
  TYPE ParametersForSN
     INTEGER :: nfile = 2, base = 150 ! for output files
     INTEGER :: Nt, Nz, Ntp1, Nzp1 ! Nt = paramSFH%Nt, Nz = paramSFH%Nz
     INTEGER :: N_DTD = 18 ! 1-5: Greggio & Renzini83(Z/Zsun=1/50,1/5,1/2.5,1,2.5),
                           ! 6-8: Delta Function (0.5, 1, 3 Gyr),
                           ! 9-13: power-law (t^-1),
                           !        from 0.05, 0.1, 0.5, 1, 3 Gyr
                           ! 14-17: Gaussian (exp[-(log10(t)-log10(t_age))^2/sigma^2],
                           !         sigma=0.3dex: t_age=0.1, 0.5, 1, 3 Gyr
                           ! 18: flat from 1 Myr--15 Gyr
     INTEGER :: iband
     DOUBLE PRECISION :: thMstar ! threshold for stellar mass
  END TYPE ParametersForSN
  TYPE(ParametersForSN) :: paramSN


  INTEGER, ALLOCATABLE :: i_SN(:) ! i_SN(i) for i=1-paramSN%nfile
  CHARACTER*100, ALLOCATABLE :: fnameSN(:) ! fnameSN(i) for i=1-paramSN%nfile


  DOUBLE PRECISION, ALLOCATABLE :: R_SN(:) ! SN rate [yr^-1]
  DOUBLE PRECISION, ALLOCATABLE :: CSNR(:) ! cosmic SN rate density [h^3/Mpc^3/yr]


  TYPE BaseQuantitiesForSN
     DOUBLE PRECISION, ALLOCATABLE :: Z(:) ! for 1-paramSN%Zp1
     DOUBLE PRECISION, ALLOCATABLE :: t_univ(:), t_lk(:) ! for 1-paramSN%Ntp1
     DOUBLE PRECISION :: tl, tu, tstep_yr
  END type BaseQuantitiesForSN
  TYPE(BaseQuantitiesForSN) :: SNbase


  TYPE DelayTimeDistribution
     CHARACTER :: name*100
     INTEGER, ALLOCATABLE :: flag(:) ! flag(i) for i=1-paramSN%Ntp1
                                     ! --- 0: r=0, 1:non-zero r
     DOUBLE PRECISION, ALLOCATABLE :: r(:) ! r(i) for i=1-paramSN%Ntp1
  END type DelayTimeDistribution
  TYPE(DelayTimeDistribution), ALLOCATABLE :: DTD(:) ! N_DTD


  TYPE ConstantsForSNrelated
     DOUBLE PRECISION :: hinv, hCube, corr
  END type ConstantsForSNrelated
  TYPE(ConstantsForSNrelated) :: constSN
END module SNrelated

