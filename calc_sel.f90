module global_var
  implicit none
  integer k
  integer :: n_prof
  integer :: num_step
  double precision :: zsp1, xmM
  double precision :: OMEGA0, hubble, lambda0, Hubble0
  double precision :: mag(5,2)
  double precision :: dprop(5000),dlum(5000),dad(5000),diz(5000)
  double precision,parameter :: KE_lam=0.d0
  double precision,parameter :: sFWHM=0.98d0,theta_min=0.255d0
  double precision,parameter :: Sth(5)=(/29.24d0,28.67d0,28.73d0,28.52d0,&
       27.61d0/)
  ! for Subaru Deep Field given by Kashikawa-san
  double precision th(5)           !isophotal radius
  double precision, parameter :: PI=acos(-1.d0)

  double precision :: A_HI(5) ! intervenning HI absorption
  integer,parameter:: Nh1=1000
  double precision :: zh1(Nh1),tauh1(5,Nh1)
end module global_var

program selection_effects
  use global_var
  implicit none
  integer :: i,j,Ngal
  integer :: loop
  integer :: pos
  real :: x
  character :: model*5,nbody*8,ci*3
  character :: filename*30, fileout*30
  integer, parameter :: N=600000
  integer :: morph(N)
  double precision :: xtmp(20)
  double precision :: mag0(5,N),size(N),zsp1ar(N)

!!$     write(*,*) 'model (5 characters) ='
!!$     read(*,'(A)') model
  call getarg(1, model)

  nbody='w390141a'

  open(20, file = '../'//nbody//'/param', status = 'old')
  read(20, *) x             ! simulation box size [Mpc/h]
  read(20, *) OMEGA0
  read(20, *) hubble
  read(20, *) x             ! Omega_b
  close(20)
  lambda0 = 1.d0 - OMEGA0

  open(1,file='distance.dat',status='old')
  read(1,*)
  read(1,*)
  do i=1,5000
     read(1,*)diz(i),dprop(i),dlum(i),dad(i)
  enddo
  close(1)

  call read_tauh1

  num_step = 55
  filename='../data/'//trim(model)//'1_001.dat'

  do loop = 1, num_step - 1
     write(ci,'(I3.3)') loop
     pos = index(filename,'.dat') - 3
     filename(pos:pos + 2) = ci

     fileout = filename
     pos = index(fileout, '.dat') - 5
     fileout(pos:pos) = 's'

     open(1,file=filename,status='old')
     read(1,*)
     read(1,*)
     i=1
10   read(1,*,end=20) morph(i),(xtmp(j),j=1,15)
     zsp1ar(i) = xtmp(1)    ! 1+z
     mag0(1:5,i) = xtmp(2:6) + 5. * log10(hubble)
     ! must be SCam-B,V,R,i,z
     size(i) = xtmp(8)      ! effective radius [kpc/h]
     i=i+1
     goto 10
20   continue
     close(1)
     Ngal=i-1

     zsp1 = zsp1ar(1)

     call interHI ! making grid for HI absorption at z
     do i = 1, Ngal
        mag0(1:5,i) = mag0(1:5,i) + A_HI(1:5)
     enddo

     print *, 'loop=',loop,'Ngal=', Ngal !, (A_HI(k),k=1,5)

     size(:)=size(:)/hubble ! [kpc/h] -> [kpc]

     open(2,file=fileout,status='unknown')
     write(2,'(A)') '# selection effects for '//model
     write(2,'(A)') '# (1-5)isophotal mag of SCam-BVRiz '//&
          '(no 5log(h)), '//&
          '(6-10) isophotal size for BVRiz, (11)m-M'
     do i=1,Ngal
        if(morph(i) == 3) then
           n_prof = 1       ! spiral
        else
           n_prof = 4       ! E/S0
        endif
        mag(1:5,1)=mag0(1:5,i)
        call sel_lcdm(size(i)) ! return: th, mag, xmM

        write(2,'(11E13.5)')&
             (mag(k,2),k=1,5), ! SCam-B,V,R,i,z&
        (th(k),k=1,5), ! isophotal size for the mags&
        xmM            ! m-M at z
     enddo
     close(2)
  enddo
end program selection_effects
!!$    ---------------------------------------------------
subroutine read_tauh1
  use global_var
  integer i

!!$     open(1,file='../opth1/tau_miyRcBVIcK.dat',status='old')
  open(1,file='../../opth1/tau_SuprimeCam5.dat',status='old')
  do i=1,Nh1
     read(1,*)zh1(i),(tauh1(j,i),j=1,5)
  enddo
  close(1)   
  tauh1(:,:)=2.5d0*tauh1(:,:)*log10(exp(1.d0))

!!$     open(1,file='../opth1/tau_SuprimeCam5.dat',status='old')
!!$     open(1,file='../opth1/tau_F775W.dat',status='old')
!!$     open(1,file='../opth1/tau_CISCO_J.dat',status='old')
!!$     open(1,file='../opth1/tau_F850LP.dat',status='old')
!!$     open(1,file='../opth1/tau_CISCO_Kp.dat',status='old')
!!$     open(1,file='../opth1/data/tau_f814w.dat',status='old')
  return
end subroutine read_tauh1
!!$    =================== subroutine interHI ================
subroutine interHI
  use global_var
  implicit none
  integer iz

  iz=max(1,int((zsp1-1.d0)*1.d+2))
  if(iz.lt.Nh1)then
     A_HI(:)=(tauh1(:,iz+1)-tauh1(:,iz))*&
          (zsp1-1.d0-zh1(iz))/(zh1(iz+1)-zh1(iz))+tauh1(:,iz)
  else
     A_HI(:)=tauh1(:,Nh1)
  endif
  return
end subroutine interHI

subroutine sel_lcdm(r_e)
  use global_var
  implicit none
  integer :: i,j
  double precision :: M_lam,z,r_e
  double precision :: S_L,theta,apm_lam
  double precision :: Mc,Zcold
  double precision :: A,Splusm
  double precision :: gauss,sig,theta1,theta2,theta3,S_1,S_2,S_FWHM

!!$     Omega0s=OMEGA0

!!$     KE_lam = 0d0
  z = zsp1-1.d0
  Hubble0 = 100.d0*hubble

  do i=1,5                  !i: B,V,R,i',z'
     call theta_cal(mag(i,1), z, r_e, Sth(i), th(i))
     if(th(i).gt.theta_min)then
        call apm_lam_cal(mag(i,1), z, r_e, th(i), mag(i,2))
     else
        mag(i,2)=100.d0
     endif
  enddo

  return
end subroutine sel_lcdm
!!$---------------------------------------------------------------------------
!!$---------------------------------------------------------------------------
Subroutine theta_cal(M_lam, z, r_e, S_L, theta)
  use global_var
  Implicit None
!!$---------------------------------------------------------------------------
!!$---------------------------------------------------------------------------
!!$input:
  double precision M_lam,z,r_e
!!$     double precision Hubble0, Omega0s, lambda0 ! Hubble0 in km/s/Mpc
  double precision S_L
  double precision theta  ! viewing angle of a galaxy above S_L, in arcsec
  double precision theta_upf
  parameter (theta_upf = 5d0) 
  real(4) theta_acc4
  parameter (theta_acc4 = 0.001d0) ! arcsec
  double precision theta_up,theta_e,S_lam
  double precision d_prop, d_curv, d_lum, d_ad ! Mp!!$
  integer i
  real(4) rtbis
  real(4) theta_solve4
  external theta_solve4
  double precision M_lam_com, KE_lam_com, z_com, r_e_com
  common /theta_com2/  M_lam_com, KE_lam_com, z_com, r_e_com
  double precision S_L_com, sFWHM_com
  common /obs_cond/ S_L_com, sFWHM_com
  double precision H_com, O_com, L_com
  common /cosmo_para/ H_com, O_com, L_com
  double precision c
  parameter (!!$= 2.9979d5) 
  real t_rt,t_rt0,t_rt1
  common /trt/t_rt,t_rt0,t_rt1

  call s_lam_cal(M_lam, z, r_e, 0d0, S_lam)

  If (S_lam .ge. S_L) then
     theta = 0d0
     return
  End if
!!$      call cos_d_cal(Hubble0, Omega0, lambda0, z, 
!!$     &     d_prop, d_curv, d_lum, d_ad)
  i=int(z*1.d+2)
  d_ad=(c/Hubble0)*&
       ((dad(i+1)-dad(i))*(z-diz(i))/(diz(i+1)-diz(i))+dad(i))

  theta_e = r_e / (d_ad*1d3) * 2.06264d5 ! arcsec
!!$     theta_up = theta_upf * theta_e ! arcsec

!!$find the upper theta braketing the solution:
!!$     Do i = 1, 5
  do i = 1, 12
     theta_up = theta_e * theta_upf** dble(i)
     call s_lam_cal(M_lam, z, r_e, theta_up, S_lam)
     if (S_lam .gt. (S_L + 1d-4)) then
        goto 1
     endif
  enddo
  print *, 'S_lam', S_lam
  print *, 'S_L', S_L
  print *, 'r_e', r_e
  print *, 'theta_e', theta_e
  write(*,*) 'S(50 r_e) should be fainter than S_L';stop
1 continue

  M_lam_com = M_lam
  KE_lam_com = KE_lam
  z_com = z
  r_e_com = r_e

!!$     H_com = Hubble0
!!$     O_com = Omega0s
!!$     L_com = lambda0

  S_L_com = S_L
  sFWHM_com = sFWHM

  theta = dble(rtbis(theta_solve4, 0e0, real(theta_up), theta_acc4) )
  return
end Subroutine theta_cal
!!$-------------------------
real(4) Function theta_solve4(theta4)
  implicit none
  real(4) theta4
  double precision S_lam ! mag arcsec^-2
  double precision M_lam_com, KE_lam_com, z_com, r_e_com
  common /theta_com2/ M_lam_com, KE_lam_com, z_com, r_e_com
  double precision S_L_com, sFWHM_com
  common /obs_cond/ S_L_com, sFWHM_com
  double precision H_com, O_com, L_com
  common /cosmo_para/ H_com, O_com, L_com

  call s_lam_cal(M_lam_com, z_com, r_e_com, dble(theta4), S_lam)
!!$     print *, 'S_lam 3', S_lam, S_L_com

  theta_solve4 = real( S_lam - S_L_com )

  return
end Function theta_solve4
!!$---------------------------------------------------------------------------
Subroutine s_lam_cal(M_lam, z, r_e, theta, S_lam)
  use global_var
  Implicit None
!!$     double precision M_lam,KE_lam,z,r_e
  double precision M_lam,z,r_e
!!$     double precision Hubble0, Omega0, lambda0 ! Hubble0 in km/s/Mpc
  double precision theta,S_lam
  double precision d_ad,d_prop, d_curv, d_lum
  double precision beta, sigma_t
  double precision sm_g_tilde
  double precision g_inf,vG_inf

  call cos_d_cal(z, d_prop, d_curv, d_lum, d_ad)

  beta = (theta*4.848d-6) * (d_ad*1d3) / r_e
  sigma_t = (sFWHM*4.848d-6) / 2.3548d0 * (d_ad*1d3) / r_e

  call sm_g_tilde_cal(beta,sigma_t,sm_g_tilde)

  vG_inf = 2d0 * PI / 1.679**2d0
  S_lam = M_lam + KE_lam + 5d0 * dlog10( r_e * 1d3 * (1d0+z)**2 / 10d0 ) &
       + 26.5721d0 - 2.5d0 * dlog10( sm_g_tilde / vG_inf)
!!$     print *, 'S_lam 0', S_lam, M_lam, KE_lam, r_e, sm_g_tilde, vG_inf
  return
end Subroutine s_lam_cal
!!$    === dop ===
Subroutine dop(M_lam, z, r_e, theta, S_lam)
  use global_var
  Implicit None
!!$     double precision M_lam,KE_lam,z,r_e
  double precision M_lam,z,r_e
!!$     double precision Hubble0, Omega0, lambda0 ! Hubble0 in km/s/Mpc
  double precision theta,S_lam,d_ad,d_prop, d_curv, d_lum
  double precision beta,sigma_t,sm_g_tilde,sm_g_tilde2
  double precision g_inf,vG_inf

  call cos_d_cal(z, d_prop, d_curv, d_lum, d_ad)

  beta = (theta*4.848d-6) * (d_ad*1d3) / r_e
  sigma_t = (sFWHM*4.848d-6) / 2.3548d0 * (d_ad*1d3) / r_e
  call sm_g_tilde_cal(beta,sigma_t, sm_g_tilde)
  call sm_g_tilde_cal(0.d0,sigma_t, sm_g_tilde2)

  vG_inf = 2d0 * PI / 1.679**2d0
  S_lam=M_lam+KE_lam + 5d0*dlog10(r_e*1d+3*(1d0+z)**2/10d0) + 26.5721d0 &
       - 2.5d0*dlog10(0.5*(sm_g_tilde+sm_g_tilde2)/vG_inf)
  return
end Subroutine dop
!!$-------------------------
double precision Function G_inf
  use global_var
  Implicit None

  if (n_prof .eq. 1) then
     G_inf = 2d0 * PI / 1.679**2d0
  elseif (n_prof .eq. 4) then
     G_inf = 40320d0 * PI / 7.670**8d0
  else
     write(*,*) 'n_prof must be 1 or 4 in G_inf';stop
  endif
  return
End Function G_inf
!!$---------------------------------------------------------------------------
Subroutine apm_lam_cal(M_lam, z, r_e, theta, apm_lam)
  use global_var
  Implicit None
  double precision M_lam,z,r_e
  double precision theta,apm_lam
  double precision d_ad,d_lum,d_prop,d_curv
  double precision beta0, sigma_t
  double precision lg_g_tilde  ! large G-tilder(beta, sigma_t)
  double precision G_inf,vG_inf

  if (theta .le. 0d0) then
     apm_lam = 1d4
     return
  endif

  call cos_d_cal(z, d_prop, d_curv, d_lum, d_ad)

  beta0 = (theta*4.848d-6) * (d_ad*1d3) / r_e
  sigma_t = (sFWHM*4.848d-6) / 2.3548d0 * (d_ad*1d3) / r_e

  vG_inf = 2d0 * PI / 1.679**2d0
  if (theta .lt. 1d10) then
     call lg_g_tilde_cal(beta0, sigma_t, lg_g_tilde)
  Else 
     lg_g_tilde = vG_inf
  End if

  apm_lam=M_lam-2.5d0*dlog10(lg_g_tilde/vG_inf)
!!$     apm_lam=M_lam-2.5d0*dlog10(lg_g_tilde/G_inf(n_prof))
!!$    &     + KE_lam 
!!$    &     + 5d0 * dlog10( d_lum*1d6 / 10d0 )
  xmM=5d0 * dlog10( d_lum*1d6 / 10d0 )
  return
end Subroutine apm_lam_cal
!!$---------------------------------------------------------------------------
Subroutine sm_g_tilde_cal(beta, sigma_t, sm_g_tilde)
  use global_var
  Implicit None
  double precision beta,sigma_t
  double precision sm_g_tilde  ! small g tilder
  double precision up_limit_sigma
!!$     parameter (up_limit_sigma = 4d0)
  parameter (up_limit_sigma = 3.5d0)
  double precision a_n(1:4)
  real(4) g_int_acc4
  parameter (g_int_acc4 = 1e-3)
  double precision up_limit_s,r,s,s_low_tmp
  real(4) s_low,s_up,g4
  real(4) sm_g_tilde_int
  external sm_g_tilde_int
  double precision beta_com, sigma_t_com
  common /sm_g_com2/ beta_com, sigma_t_com
  common /sm_int8/r

  a_n(1) = 1.679d0; a_n(4) = 7.670d0

  beta_com = beta; sigma_t_com = sigma_t
  r=beta/sigma_t
!!$---- when sigma_t is much smaller than beta:
  if ( r .gt. 1d5 ) then
     sm_g_tilde = dexp( -1.679d0*beta )
     return
  endif
!!$---- integration range determination:

!!$     up_limit_s=(up_limit_sigma**2/a_n(n_prof))**n_prof/sigma_t
!!$     up_limit_s=(up_limit_sigma**2/1.679d0)/sigma_t
  up_limit_s=(40.d0/1.679d0)/sigma_t
  s_low_tmp=r-sqrt(2.d0)*up_limit_sigma
  if(s_low_tmp.gt.up_limit_s)then
     sm_g_tilde=1.d-36
     return
  endif
  s_low=real(max(0.d0,s_low_tmp))
  s_up=real(min(r+sqrt(2.d0)*up_limit_sigma,up_limit_s))

  call myqromb(sm_g_tilde_int, s_low, s_up, g4, g_int_acc4)

  sm_g_tilde = dble( g4 )

  return
end Subroutine sm_g_tilde_cal
!!$------------
real(4) function sm_g_tilde_int(xi_prime4)
  Implicit None
  real(4) xi_prime4
  double precision x,bess_part
  double precision beta_com, sigma_t_com
  common /sm_g_com2/ beta_com, sigma_t_com
  real(4) bessi0mod
  double precision r,s
  common /sm_int8/r

  s = dble(xi_prime4)
  x = r*s
  bess_part = dble( bessi0mod(real(x)))

  sm_g_tilde_int = real(s*dexp( -1.679d0*(sigma_t_com*s) ) * bess_part &
       * dexp(-0.5d0*(r-s)**2))

  return
end function sm_g_tilde_int
!!$---------------------------------------------------------------------------
!!$---------------------------------------------------------------------------
Subroutine lg_g_tilde_cal(beta, sigma_t, lg_g_tilde)
  use global_var
  Implicit None
!!$---------------------------------------------------------------------------
!!$    large g_tilder calculation
!!$---------------------------------------------------------------------------
!!$input:
  double precision beta
  double precision sigma_t
!!$output:
  double precision lg_g_tilde  ! large g tilder
!!$const:
  real(4) g_int_acc4
  parameter (g_int_acc4 = 1e-3)
!!$local:
  real(4) g4
  real(4) xi_prime_up4
!!$function:
  real(4) lg_g_tilde_int4
  external lg_g_tilde_int4
!!$common:
  double precision sigma_t_com
  common /lg_g_com2/ sigma_t_com
!!$begin:
  sigma_t_com = sigma_t

  xi_prime_up4 = real( beta )
  call myqromb2(lg_g_tilde_int4, 0e0, xi_prime_up4, g4, g_int_acc4)

  lg_g_tilde = dble(g4)
  return
end Subroutine lg_g_tilde_cal
!!$------------
real(4) function lg_g_tilde_int4(xi_prime4)
  use global_var
  Implicit None
!!$input:
  real(4) xi_prime4  ! xi_prime = xi**(1/n_prof)
!!$local:
  real(4) xi4
  double precision sm_g_tilde
  double precision lg_g8
!!$common:
  double precision sigma_t_com
  common /lg_g_com2/ sigma_t_com
!!$function:
!!$begin:
  xi4 = xi_prime4
  call sm_g_tilde_cal(dble(xi4), sigma_t_com, sm_g_tilde)
  lg_g8 = 2d0 * PI * sm_g_tilde * dble(xi4)
  If (lg_g8 .le. 1d-30) then
     lg_g_tilde_int4 = 1e-30
  Else
     lg_g_tilde_int4 = real( lg_g8 )
  End if

  return
end function lg_g_tilde_int4
!!$---------------------------------------------------------------------------
!!$---------------------------------------------------------------------------
Subroutine cos_d_cal(z, d_prop, d_curv, d_lum, d_ad)
  use global_var
  Implicit None
!!$---------------------------------------------------------------------------
!!$    calculation of proper, luminosity, and angular diameter distance [Mpc]
!!$    from cosmological parameters and redshift
!!$---------------------------------------------------------------------------
!!$input:
  double precision z                  ! redshift
  double precision d_prop  ! Proper Distance [Mpc]     
  real(4) d_prop4 ! [se!!$Mp!!$/ km]
  double precision d_curv ! = R_0 * sin(chi)       if k = 1
  double precision d_lum ! Luminosity Distance [Mpc]     
  double precision d_ad ! Angular Diameter Distance [Mpc]     
  double precision delta_chi  !  disrance in comoving coodinate
  double precision R_0  ! present scale factor [Mpc]
  double precision !!$ ! speed of light [km/sec]
  parameter (!!$= 2.9979d5)
  integer i

  i=int(z*1.d+2)
  d_prop=(c/Hubble0)*&
       ((dprop(i+1)-dprop(i))*(z-diz(i))/(diz(i+1)-diz(i))+dprop(i))
!!$flat:
  d_curv = d_prop
  d_lum = d_curv * (1 + z); d_ad = d_curv / (1 + z)

  return
end Subroutine cos_d_cal
!!$    ---
SUBROUTINE myqromb(func,a,b,ss, EPS)
  INTEGER JMAX,JMAXP,K,KM
  REAL a,b,func,ss,EPS
  EXTERNAL func
!!$     PARAMETER (EPS=1.e-6, JMAX=20, JMAXP=JMAX+1, K=5, KM=K-1)
  PARAMETER (JMAX=20, JMAXP=JMAX+1, K=5, KM=K-1)
!!$U    USES polint,trapzd
  INTEGER j
  REAL dss,h(JMAXP),s(JMAXP)

  h(1)=1.
  do j=1,JMAX
     call trapzd(func,a,b,s(j),j)
     if (j.ge.K) then
        call polint(h(j-KM),s(j-KM),K,0.,ss,dss)
        if (abs(dss).le.EPS*abs(ss))then
           t_qrom=t_qrom+t_qrom1-t_qrom0
           return
        endif
     endif
     s(j+1)=s(j)
     h(j+1)=0.25*h(j)
  enddo
  write(*,*) 'too many steps in myqromb';stop
END SUBROUTINE myqromb
!!$ (C) Copr. 1986-92 Numerical Recipes Software z!0(0.
!!$    ---
SUBROUTINE myqromb2(func,a,b,ss, EPS)
  INTEGER JMAX,JMAXP,K,KM
  REAL a,b,func,ss,EPS
  EXTERNAL func
!!$     PARAMETER (EPS=1.e-6, JMAX=20, JMAXP=JMAX+1, K=5, KM=K-1)
  PARAMETER (JMAX=20, JMAXP=JMAX+1, K=5, KM=K-1)
!!$U    USES polint,trapzd
  INTEGER j
  REAL dss,h(JMAXP),s(JMAXP)

  h(1)=1.
  do j=1,JMAX
     call trapzd2(func,a,b,s(j),j)
     if (j.ge.K) then
        call polint(h(j-KM),s(j-KM),K,0.,ss,dss)
        if (abs(dss).le.EPS*abs(ss))then
           t_qrom=t_qrom+t_qrom1-t_qrom0
           return
        endif
     endif
     s(j+1)=s(j)
     h(j+1)=0.25*h(j)
  enddo
  write(*,*) 'too many steps in myqromb2';stop
END SUBROUTINE myqromb2
!!$ (C) Copr. 1986-92 Numerical Recipes Software z!0(0.
!!$    ---
SUBROUTINE trapzd2(func,a,b,s,n)
  INTEGER n
  REAL a,b,s,func
  EXTERNAL func
  INTEGER it,j
  REAL del,tot,tnm,x
  if (n.eq.1) then
     s=0.5*(b-a)*(func(a)+func(b))
  else
     it=2**(n-2)
     tnm=it
     del=(b-a)/tnm
     x=a+0.5*del
     tot=0.
     do j=1,it
        tot=tot+func(x)
        x=x+del
     enddo
     s=0.5*(s+(b-a)*tot/tnm)
  endif
  return
END SUBROUTINE trapzd2
!!$ (C) Copr. 1986-92 Numerical Recipes Software z!0(0.

!!$followings are the Numerical Recipes routines.
!!$    ---
FUNCTION bessi0(x)
  REAL bessi0,x
  REAL ax
  DOUBLE PRECISION p1,p2,p3,p4,p5,p6,p7,q1,q2,q3,q4,q5,q6,q7,q8,q9,y
  SAVE p1,p2,p3,p4,p5,p6,p7,q1,q2,q3,q4,q5,q6,q7,q8,q9
  DATA p1,p2,p3,p4,p5,p6,p7/1.0d0,3.5156229d0,3.0899424d0,1.2067492d0,&
       0.2659732d0,0.360768d-1,0.45813d-2/
  DATA q1,q2,q3,q4,q5,q6,q7,q8,q9/0.39894228d0,0.1328592d-1,0.225319d-2,&
       -0.157565d-2,0.916281d-2,-0.2057706d-1,0.2635537d-1,-0.1647633d-1,&
       0.392377d-2/
  if (abs(x).lt.3.75) then
     y=(x/3.75)**2
     bessi0=p1+y*(p2+y*(p3+y*(p4+y*(p5+y*(p6+y*p7)))))
  else
     ax=abs(x)
     y=3.75/ax
     bessi0=(exp(ax)/sqrt(ax))*(q1+y*(q2+y*(q3+y*(q4+y*(q5+y*(q6+y*(q7&
          +y*(q8+y*q9))))))))
  endif
  return
END FUNCTION bessi0
!!$    ---
FUNCTION rtbis(func,x1,x2,xacc)
  INTEGER JMAX
  REAL rtbis,x1,x2,xacc,func
  EXTERNAL func
  PARAMETER (JMAX=40)
  INTEGER j
  REAL dx,f,fmid,xmid

  fmid=func(x2)
  f=func(x1)

  If(f*fmid.ge.0.) then
     write (*, *) x1, f
     write (*, *) x2, fmid
     write(*,*) 'root must be bracketed in rtbis';stop
  End if
  if(f.lt.0.)then
     rtbis=x1
     dx=x2-x1
  else
     rtbis=x2
     dx=x1-x2
  endif
  do j=1,JMAX
     dx=dx*.5
     xmid=rtbis+dx
     fmid=func(xmid)
     if(fmid.le.0.)rtbis=xmid
     if(abs(dx).lt.xac!!$.or. fmid.eq.0.)return
  enddo
  write(*,*) 'too many bisections in rtbis';stop
END FUNCTION rtbis
!!$ (C) Copr. 1986-92 Numerical Recipes Software z!0(0.
!!$    ---
SUBROUTINE trapzd(func,a,b,s,n)
  INTEGER n
  REAL a,b,s,func
  EXTERNAL func
  INTEGER it,j
  REAL del,tot,tnm,x

  if (n.eq.1) then
     s=0.5*(b-a)*(func(a)+func(b))
  else
     it=2**(n-2)
     tnm=it
     del=(b-a)/tnm
     x=a+0.5*del
     tot=0.
     do j=1,it
        tot=tot+func(x)
        x=x+del
     enddo
     s=0.5*(s+(b-a)*tot/tnm)
  endif
  return
END SUBROUTINE trapzd
!!$ (C) Copr. 1986-92 Numerical Recipes Software z!0(0.
!!$    ---
SUBROUTINE polint(xa,ya,n,x,y,dy)
  INTEGER n,NMAX
  REAL dy,x,y,xa(n),ya(n)
  PARAMETER (NMAX=10)
  INTEGER i,m,ns
  REAL den,dif,dift,ho,hp,w,c(NMAX),d(NMAX)

  ns=1
  dif=abs(x-xa(1))
  do i=1,n
     dift=abs(x-xa(i))
     if (dift.lt.dif) then
        ns=i
        dif=dift
     endif
     c(i)=ya(i)
     d(i)=ya(i)
  enddo
  y=ya(ns)
  ns=ns-1
  do m=1,n-1
     do i=1,n-m
        ho=xa(i)-x
        hp=xa(i+m)-x
        w=c(i+1)-d(i)
        den=ho-hp
        if(den.eq.0.)then
           write(*,*)'failure in polint';stop
        endif
        den=w/den
        d(i)=hp*den
        c(i)=ho*den
     enddo
     if (2*ns.lt.n-m)then
        dy=c(ns+1)
     else
        dy=d(ns)
        ns=ns-1
     endif
     y=y+dy
  enddo
  return
END SUBROUTINE polint
!!$ (C) Copr. 1986-92 Numerical Recipes Software z!0(0.
!!$    --------------------------- function gauss(x) -------------------
function gauss(x)
  double precision myerf,xx
  double precision gauss,x,gus0,gus1,gus2
  if(x.lt.0.)then
     xx=(abs(x))
     if(x.gt.-.5)then
        gauss=.5d0+gus0(xx)
     else if(x.gt.-1.5)then
        gauss=.5d0+gus1(xx)
     else if(x.gt.-2.5)then
        gauss=.5d0+gus2(xx)
     else 
        gauss=.5d0+(myerf(xx*.7071067690849304D+00)/2.D0)
     endif
  else 
     xx=x
     if(x.lt..5)then
        gauss=.5d0-gus0(xx)
     else if(x.lt.1.5)then
        gauss=.5d0-gus1(xx)
     else if(x.lt.2.5)then
        gauss=.5d0-gus2(xx)
     else
        gauss=.5d0-(myerf(xx*.7071067690849304D+00)/2.D0)
     endif
  endif
  return
end function gauss
!!$    ---------------------------- function gus0,1,2 -----------------
function gus0(t)
  double precision t,f0,f1
  double precision gus0

  f1=t**2*(0.2976190531626344D-02-t**2*(0.2893518540076911D-03-t**2&
       *0.2367424167459831D-04))
  f0=t*(1.d0-t**2*(0.1666666716337204D+00-t**2*(0.2500000037252903D-01&
       -f1)))*0.3989422804014327D+00
  gus0=f0
  return
end function gus0

function gus1(t)
  double precision h0,h1,h2,h3,h4,t
  double precision gus1

  h0=0.3989422804014327d0
  h1=0.855624391892149d0+(t-1.d0)*(0.6065306597126334d0&
       -0.3032653298563167d0*(t-1.d0))
  h2=(t-1.d0)**4*(0.05054422164271945d0-0.01010884432854389d0*(t-1.d0)&
       -0.005054422164271946d0*(t-1.d0)**2)
  h3=(t-1.d0)**7*(0.001925494157817884d0+0.0003008584621590444d0*(t-1.d0))
  h4=(t-1.d0)**9*(-0.0002206295389166325d0-4.680020522474028d-6*(t-1.d0)&
       +0.00001847696414067667d0*(t-1.d0)**2)
  gus1=(h0*(h1+h2+h3+h4))
  return 
end function gus1

function gus2(t)
  double precision a0,a1,a2,a3,a4,t
  double precision gus2

  a0=0.3989422804014327d0
  a1=1.196288013322608d0+(t-2.d0)*(0.1353352832366127d0&
       -0.1353352832366127d0*(t-2.d0))
  a2=(t-2.d0)**3*(0.06766764161830636d0-0.01127794026971772d0*(t-2.d0))
  a3=(t-2.d0)**5*(-0.005638970134858862d0+0.003383382080915318d0*(t-2.d0)&
       -0.0002953746261116547d0*(t-2.d0)**2)
  a4=(t-2.d0)**8*(-0.0002886615664272989d0+0.0000928639923002551d0*(t-2.d0)&
       +7.086007444597778d-6*(t-2.d0)**2)-8.88632799631137d-6*(t-2.d0)**11
  gus2=(a0*(a1+a2+a3+a4))
  return
end function gus2

function erfcc(x)
  double precision erfcc,x
  double precision t,z

  z=dabs(x)
  t=1.d0/(1.d0+0.5d0*z)
  erfcc=t*dexp(-z*z-1.26551223d0+t*(1.00002368d0+t*(.37409196d0&
       +t*(.09678418d0+t*(-.18628806d0+t*(.27886807d0+t*(-1.13520398d0&
       +t*(1.48851587d0+t*(-.82215223d0+t*.17087277d0)))))))))
  if(x.lt.0.)erfcc=2.d0-erfcc
  return
end function erfcc

function myerf(x)
  double precision myerf,x
  double precision gammp

  if(x.lt.0.d0)then
     myerf=-gammp(.5d0,x**2)
  else
     myerf=gammp(.5d0,x**2)
  endif
  return
end function myerf

function gammp(a,x)
  double precision a,gammp,x
  double precision gammcf,gamser,gln
  if(x.lt.0.d0.or.a.le.0.d0)then
     write(*,*)'bad arguments in gammp';stop
  endif
  if(x.lt.a+1.d0)then
     call gser(gamser,a,x,gln)
     gammp=gamser
  else
     call gcf(gammcf,a,x,gln)
     gammp=1.d0-gammcf
  endif
  return
end function gammp

function gammq(a,x)
  double precision a,gammq,x
  double precision gammcf,gamser,gln

  if(x.lt.0.d0.or.a.le.0.d0)then
     write(*,*)'bad arguments in gammq';stop
  endif
  if(x.lt.a+1.d0)then
     call gser(gamser,a,x,gln)
     gammq=1.d0-gamser
  else
     call gcf(gammcf,a,x,gln)
     gammq=gammcf
  endif
  return
end function gammq

subroutine gser(gamser,a,x,gln)
  integer ITMAX
  double precision a,gamser,gln,x,eps
  parameter(ITMAX=100,eps=3.d-7)
  integer n
  double precision ap,del,tot,gammln
  gln=gammln(a)
  if(x.le.0.d0)then
     if(x.lt.0.d0)then
        write(*,*)'x<0 in gser';stop
     endif
     gamser=0.d0
     return
  endif
  ap=a
  tot=1.d0/a
  del=tot
  do n=1,ITMAX
     ap=ap+1.d0
     del=del*x/ap
     tot=tot+del
     if(dabs(del).lt.dabs(tot)*eps)goto 1
  enddo
  write(*,*)'a too large, ITMAX too small in gser';stop
1 gamser=tot*dexp(-x+a*dlog(x)-gln)
  return
end subroutine gser

subroutine gcf(gammcf,a,x,gln)
  integer ITMAX
  double precision a,gammcf,gln,x,eps,FPMIN
  parameter(ITMAX=100,eps=3.d-7,FPMIN=1.d-30)
  integer i
  double precision an,b,c,d,del,h,gammln

  gln=gammln(a)
  b=x+1.d0-a
  c=1.d0/FPMIN
  d=1.d0/b
  h=d
  do i=1,ITMAX
     an=-dble(i)*(dble(i)-a)
     b=b+2.d0
     d=an*d+b
     if(dabs(d).lt.FPMIN)d=FPMIN
     c=b+an/c
     if(dabs(c).lt.FPMIN)c=FPMIN
     d=1.d0/d
     del=d*c
     h=h*del
     if(dabs(del-1.d0).lt.eps)goto 1
  enddo
  write(*,*)'a too large, ITMAX too small in gcf';stop
1 gammcf=dexp(-x+a*dlog(x)-gln)*h
  return
end subroutine gcf

function gammln(xx)
  double precision gammln,xx
  integer j
  double precision ser,stp,tmp,x,y,cof(6)
  save cof,stp
  data cof,stp/76.18009172947146d0,-86.50532032941677d0,24.01409824083091d0,&
       -1.231739572450155d0,.1208650973866179d-2,-.5395239384953d-5,&
       2.5066282746310005d0/

  x=xx
  y=x
  tmp=x+5.5d0
  tmp=(x+0.5d0)*dlog(tmp)-tmp
  ser=1.000000000190015d0
  do j=1,6 
     y=y+1.d0
     ser=ser+cof(j)/y
  enddo
  gammln=tmp+dlog(stp*ser/x)
  return
end function gammln

FUNCTION bessi0mod(x)
  REAL bessi0mod,x
  REAL ax
  DOUBLE PRECISION p1,p2,p3,p4,p5,p6,p7,q1,q2,q3,q4,q5,q6,q7,q8,q9,y
  SAVE p1,p2,p3,p4,p5,p6,p7,q1,q2,q3,q4,q5,q6,q7,q8,q9
  DATA p1,p2,p3,p4,p5,p6,p7/1.0d0,3.5156229d0,3.0899424d0,1.2067492d0,&
       0.2659732d0,0.360768d-1,0.45813d-2/
  DATA q1,q2,q3,q4,q5,q6,q7,q8,q9/0.39894228d0,0.1328592d-1,0.225319d-2,&
       -0.157565d-2,0.916281d-2,-0.2057706d-1,0.2635537d-1,-0.1647633d-1,&
       0.392377d-2/
  if (abs(x).lt.3.75) then
     y=(x/3.75)**2
     bessi0mod=exp(-x)*(p1+y*(p2+y*(p3+y*(p4+y*(p5+y*(p6+y*p7))))))
!!$       bessi0E=p1+y*(p2+y*(p3+y*(p4+y*(p5+y*(p6+y*p7)))))
  else
     ax=abs(x)
     y=3.75/ax
!!$       bessi0E=(exp(ax)/sqrt(ax))*(q1+y*(q2+y*(q3+y*(q4+y*(q5+y*(q6+y*
     bessi0mod=(1./sqrt(ax))*(q1+y*(q2+y*(q3+y*(q4+y*(q5+y*(q6+y*&
          (q7+y*(q8+y*q9))))))))
  endif
  return
END FUNCTION bessi0mod

