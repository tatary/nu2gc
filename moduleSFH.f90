module SFHrelated
  TYPE ParametersForSFH
     INTEGER :: div_type = 2 ! 1:equal # of bin, 2:equal timestep
     INTEGER :: nfile = 6, base = 100 ! for output files
     INTEGER :: loop
     INTEGER :: Nt = 100 ! 100 (high res.), 50 (low res.): used only for paramSFH%div_type=1
     INTEGER :: Nz = 6  ! BC03: 6
     INTEGER :: Ntp1, Nzp1
     INTEGER :: N_AV = 10 ! 50 (high res.), 30 (low res.)
     DOUBLE PRECISION :: stepAV = 0.2d0
     DOUBLE PRECISION :: SFRunit
  END TYPE ParametersForSFH
  TYPE(ParametersForSFH) :: paramSFH
  INTEGER, ALLOCATABLE :: i_SFH(:) ! i_SFH(i) for i=1-paramSFH%nfile
  CHARACTER(LEN=100), ALLOCATABLE :: fnameSFH(:) ! fnameSFH(i) for i=1-paramSFH%nfile


  TYPE BaseQuantitiesForSFH
     DOUBLE PRECISION, ALLOCATABLE :: Z(:), t(:), t_yr(:), tlk_yr(:)
     DOUBLE PRECISION :: zp1l, zp1u, tl, tu, tstep, tstep_yr
  END TYPE BaseQuantitiesForSFH
  TYPE(BaseQuantitiesForSFH) :: SFHbase


  DOUBLE PRECISION, ALLOCATABLE :: mSFH(:,:)
  DOUBLE PRECISION, ALLOCATABLE :: CSFH(:,:,:)

  DOUBLE PRECISION, ALLOCATABLE :: reff(:), tauV_SFH(:)
END module SFHrelated

