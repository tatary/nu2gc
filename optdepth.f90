!!$    estimation of optical depths relative to V-band optical depth
!!$    xi(i)=tau(i)/tau(V)
!!$    This function is taken from Pei (1992).
SUBROUTINE optdepth(zsp1, ionum, iforest, nz)
  use global_var
  implicit none

  DOUBLE PRECISION, INTENT(IN) :: zsp1 ! = param%zsp1
  INTEGER, INTENT(IN) :: ionum, iforest, nz
  INTEGER :: iband
  DOUBLE PRECISION :: x, xi_Vrest
  DOUBLE PRECISION, ALLOCATABLE :: xp(:)
  DOUBLE PRECISION :: Pei92, Calzetti00 ! functions


  ! allocate(xp(param%nwp1))
  allocate(xp(param%nwave))
  DO iband = 1, param%nwave
     xp(iband) = ssp%lam_c(param%iwave(iband))
  ENDDO
  ! xp(param%nwp1) = 0.55121d0 ! (central wavelength of V-band [um])

  IF(iforest == 0 .and. inode == 0 .and. nz == 1 .and. param%run_type /= 3) THEN
     IF(param%exttype == 1) THEN
        print '(A)', '# MW type extinction curve is adopted'
        write(ionum, '(A)') '# MW type extinction curve is adopted'
     ELSEIF(param%exttype == 2) THEN
        print '(A)', '# LMC type extinction curve is adopted'
        write(ionum, '(A)') '# LMC type extinction curve is adopted'
     ELSEIF(param%exttype == 3) THEN
        print '(A)', '# SMC type extinction curve is adopted'
        write(ionum, '(A)') '# SMC type extinction curve is adopted'
     ELSEIF(param%exttype == 4) THEN
        print '(A)', '# Calzetti-law is adopted'
        write(ionum, '(A)') '# Calzetti-law is adopted'
     ENDIF
  ENDIF

  DO iband = 1, param%nwave
     x = xp(iband) ! rest-frame
     IF(param%iwave(iband) <= NWAVE_ALL) THEN ! observer-frame
        IF(param%iwave(iband) < 46 .or. param%iwave(iband) > 54) x = xp(iband) / zsp1
        ! --- param%iwave(iband) = 46-54 is omitted here since they are always
        !      luminosities at rest-frame wavelength
        !      (46:NLyC, 47:L1216, 48:L1400, 49:L1500, 50:L1600, 51:L1700,
        !       52:L2800, 53:L4861, 54:L6563)
     ENDIF

     IF(param%exttype == 4) THEN ! Calzetti-law
        xi(iband) = Calzetti00(x)
     ELSE ! Pei92 extinction curve
        xi(iband) = Pei92(x)
     ENDIF

     IF(xi(iband) < 0.d0) xi(iband) = 0.d0
  ENDDO
  xi_Vrest = xi(2)

!!$  x = xp(param%nwp1)
!!$  IF(param%exttype == 4) THEN ! Calzetti-law
!!$     xi(param%nwp1) = Calzetti00(x)
!!$  ELSE ! Pei92 extinction curve
!!$     xi(param%nwp1) = Pei92(x)
!!$  ENDIF

!!$  print '(A, F8.5, A, F10.4)', '# lamV = ', xp(param%nwp1), &
!!$       ' [um], xi(lamV) = ', xi(param%nwp1)
!!$  write(ionum, '(A, F8.5, A, F10.4)') '# lamV = ', xp(param%nwp1), &
!!$       ' [um], xi(lamV) = ', xi(param%nwp1)
  IF(iforest == 0 .and. inode == 0 .and. nz == 1 .and. param%run_type /= 3) THEN
     print '(A, F8.5, A, F10.5)', '# lamV = ', xp(2), ' [um], xi(lamV) = ', xi_Vrest
     write(ionum, '(A, F8.5, A, F10.5)') '# lamV = ', xp(2), ' [um], xi(lamV) = ', xi_Vrest
  ENDIF
  DO iband = 1, param%nwave
     ! xi(iband) = xi(iband) / xi(param%nwp1)
     xi(iband) = xi(iband) / xi_Vrest
     IF(iforest == 0) THEN
        x = xp(iband) ! rest-frame
        IF(param%iwave(iband) <= NWAVE_ALL) THEN ! observer-frame
           IF(param%iwave(iband) < 46 .or. param%iwave(iband) > 54) x = xp(iband) / zsp1
        ENDIF
        IF(inode == 0 .and. nz == 1 .and. param%run_type /= 3)THEN
           print '(A, F8.5, A, F10.5, A)', '   --- xi(', x, '[um])/xi(lamV) = ', xi(iband), &
                ' ('//trim(ssp%bandname(param%iwave(iband)))//'-band)'
           write(ionum, '(A, F8.5, A, F10.5, A)') '# --- xi(', x, '[um])/xi(lamV) = ', xi(iband), &
                ' ('//trim(ssp%bandname(param%iwave(iband)))//'-band)'
        ENDIF
     ENDIF
  ENDDO

  deallocate(xp)
END SUBROUTINE optdepth
!!$=============================================================================
DOUBLE PRECISION FUNCTION Pei92(lam)
  use global_var
  implicit none
  DOUBLE PRECISION, INTENT(IN) :: lam
  INTEGER :: i
  DOUBLE PRECISION :: tmp
  TYPE dustext
     DOUBLE PRECISION :: a(6), lam(6), b(6), n(6)
  END TYPE dustext
  TYPE(dustext) :: dust
  ! --- MW type extinction curve parameters
  DOUBLE PRECISION :: a_MW(6)   = (/165.d0,14.d0,0.045d0,0.002d0,0.002d0,0.012d0/)
  DOUBLE PRECISION :: lam_MW(6) = (/0.047d0,0.08d0,0.22d0,9.7d0,18.d0,25.d0/)
  DOUBLE PRECISION :: b_MW(6)   = (/90.d0,4.00d0,-1.95d0,-1.95d0,-1.8d0,0.d0/)
  DOUBLE PRECISION :: n_MW(6)   = (/2.d0,6.5d0,2.d0,2.d0,2.d0,2.d0/)
  ! --- LMC type extinction curve parameters
  DOUBLE PRECISION :: a_LMC(6)   = (/175.d0,19.d0,0.023d0,0.005d0,0.006d0,0.020d0/)
  DOUBLE PRECISION :: lam_LMC(6) = (/0.046d0,0.08d0,0.22d0,9.7d0,18.d0,25.d0/)
  DOUBLE PRECISION :: b_LMC(6)   = (/90.d0,5.50d0,-1.95d0,-1.95d0,-1.8d0,0.d0/)
  DOUBLE PRECISION :: n_LMC(6)   = (/2.d0,4.5d0,2.d0,2.d0,2.d0,2.d0/)
  ! --- SMC type extinction curve parameters
  DOUBLE PRECISION :: a_SMC(6)   = (/185.d0,27.d0,0.005d0,0.010d0,0.012d0,0.030d0/)
  DOUBLE PRECISION :: lam_SMC(6) = (/0.042d0,0.08d0,0.22d0,9.7d0,18.d0,25.d0/)
  DOUBLE PRECISION :: b_SMC(6)   = (/90.d0,5.50d0,-1.95d0,-1.95d0,-1.8d0,0.d0/)
  DOUBLE PRECISION :: n_SMC(6)   = (/2.d0,4.0d0,2.d0,2.d0,2.d0,2.d0/)

  IF(param%exttype == 1) THEN
     dust%a(:) = a_MW(:); dust%lam(:) = lam_MW(:)
     dust%b(:) = b_MW(:); dust%n(:)   = n_MW(:)
  ELSEIF(param%exttype == 2) THEN
     dust%a(:) = a_LMC(:); dust%lam(:) = lam_LMC(:)
     dust%b(:) = b_LMC(:); dust%n(:)   = n_LMC(:)
  ELSEIF(param%exttype == 3) THEN
     dust%a(:) = a_SMC(:); dust%lam(:) = lam_SMC(:)
     dust%b(:) = b_SMC(:); dust%n(:)   = n_SMC(:)
  ENDIF

  Pei92 = 0.d0
  DO i = 1, 6
     tmp = (lam / dust%lam(i))**dust%n(i)
     Pei92 = Pei92 + dust%a(i) / (tmp + 1.d0 / tmp + dust%b(i))
  ENDDO
END FUNCTION Pei92
!!$=============================================================================
DOUBLE PRECISION FUNCTION Calzetti00(lam)
  DOUBLE PRECISION, INTENT(IN) :: lam
  DOUBLE PRECISION, PARAMETER :: RV = 4.05d0

!!$           IF(lam < 0.12d0 .or. lam > 2.20d0) THEN
!!$              print '(A, F6.3, A)', ' # Corresponding attenuation law '//
!!$    $              'cannot be given at wavelength of ', lam, '[micron]'
!!$              Calzetti = 0.d0
!!$           ELSEIF(lam >= 0.12d0 .and. lam <= 0.63d0) THEN
  IF(lam <= 0.63d0) THEN
     Calzetti00 = RV + 2.659d0 * (-2.156d0 + (1.509d0 &
          + (-0.198d0 + 0.011d0 / lam) / lam) / lam)
!!$           ELSEIF(lam <= 2.20d0) THEN
  ELSE
     Calzetti00 = RV + 2.659d0 * (-1.857d0 + 1.040d0 / lam)
  ENDIF
END FUNCTION Calzetti00
!!$=============================================================================
