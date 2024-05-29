!!$  ----------------- function z2t -------------------------------
DOUBLE PRECISION FUNCTION z2t(zp1)
  use global_var
  DOUBLE PRECISION, INTENT(IN) :: zp1
  DOUBLE PRECISION :: w
  DOUBLE PRECISION :: Cube ! function

  w   = param%OMEGA_rat / Cube(zp1)
  z2t = log(1.d0 + 2.d0 * w + 2.d0 * sqrt(w * (1.d0+w))) &
        / (3.d0 * sqrt(1.d0 - param%OMEGA0))
END FUNCTION z2t
!!$  ----------------- function dens ------------------------------
DOUBLE PRECISION FUNCTION dens(zp1)
  ! ref. Bryan & Norman (1998), Andreu Benson PhD-Thesis
  use global_var
  DOUBLE PRECISION, INTENT(IN) :: zp1
  DOUBLE PRECISION :: Ez2, OmegaMz, x, Delta_c
  DOUBLE PRECISION :: Cube ! function

  Ez2     = param%OMEGA_L + param%OMEGA0 * Cube(zp1)
  OmegaMz = 1.d0 - param%OMEGA_L / Ez2
  x       = OmegaMz - 1.d0
  Delta_c = 18.d0 * const%PI2 + (82.d0 - 39.d0 * x) * x
  dens = Delta_c * Ez2
END FUNCTION dens
!!$  -------------- function taustar -----------
DOUBLE PRECISION FUNCTION taustar(Vctmp, Rdisk, zplus1, Mcold, McoldZ, b_or_q)
  use global_var
  DOUBLE PRECISION, INTENT(IN) :: Vctmp, Rdisk, zplus1, Mcold, McoldZ
  INTEGER, INTENT(IN) :: b_or_q
  DOUBLE PRECISION :: dens
  DOUBLE PRECISION :: Zcold, taudust, em, re
  DOUBLE PRECISION :: CalDiskReff ! function

  IF(param%SFmodel == 1) THEN ! CSF
    taustar = param%tau0st &
         * (1.d0 + (Vctmp / param%Vst)**param%alpst) ! current
         ! * (Vctmp/param%Vhot(b_or_q))**param%alpst ! N05
  ELSEIF(param%SFmodel == 2) THEN ! DSF
    taustar = (Rdisk / Vctmp / param%h) * 0.9784d0 & ! [Gyr]
         / param%eps_SF &
         * (1.d0 + (Vctmp / param%Vst)**param%alpst) & ! dyn. time of galaxy
         / param%th
    ! taustar = param%tau0st &
    !     * (1.d0 + (Vctmp / param%Vhot(b_or_q))**param%alpst) &
    !     * sqrt(dens(0.d0) / dens(zplus1)) ! dyn. time of DM halo
  ELSEIF(param%SFmodel == 3) THEN ! Makiya14
    Zcold = 0.d0
    IF(Mcold > 0.d0) Zcold = McoldZ / Mcold / const%Zsun 
    em = param%emin * exp(-Zcold/param%Zch)

    re = CalDiskReff(Rdisk, 0.d0) * 1.d3 / param%h ! pc 
    taudust = 2.d-3 * (Mcold * param%munit * Zcold / (re * re)) ! Msun Zsun pc^-2 


    IF(Zcold <= 0.d0) THEN
       taustar = 1.d0 / em 
    ELSE
       taustar = 1.d0 / (param%emax * exp(-param%taud_th/taudust) + em) ! [Gyr]
    ENDIF
    taustar = min(taustar, 1.d+10)
  ENDIF
END FUNCTION taustar
!!$  ----------------- function delc ------------------------------
DOUBLE PRECISION FUNCTION delc(zp1)
  use global_var
  DOUBLE PRECISION, INTENT(IN) :: zp1
  DOUBLE PRECISION :: ww, xwv
  DOUBLE PRECISION :: Cube ! function

  ww  = param%OMEGA_rat / Cube(zp1)
  xwv = 1.d0 + 1.d0 / (5.1066d0 * ww**(-1.0812d0) + 2.0215d0 * ww**(-0.35396d0))
  delc = const%delc0 * zp1 * xwv
END FUNCTION delc
