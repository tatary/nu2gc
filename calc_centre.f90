module global_var
  integer num
!!$     integer,parameter:: N=300
  integer,parameter:: N=3000
  real lambda(N),response(N)
end module global_var

program main
  use global_var
  implicit none
  integer i,j
  real est_pII,p_ira,t0,t1,e_ira,est_pII2,e_ira_tmp,p_ira2
  external est_pII,est_pII2

  i=1
!!$     open(1,file='../SuprimeCam_z.dat',status='old')
!!$     open(1,file='../GALEX/GALEX-FUV_EA-BandPass-DS-200309.txt',
!!$     open(1,file='../GALEX/GALEX-NUV_EA-BandPass-DS-200309.txt',
!!$     open(1,file='../2MASS/2MASS_Kstot.dat',
!!$     open(1,file='../filters/filter12m.dat',
!!$     open(1,file='../filters/filter14m.dat',
!!$     open(1,file='../ACS/F775W.dat',
!!$     open(1,file='../ACS/F850LP.dat',
!!$     open(1,file='../CISCO/filter_Kp.dat',
  open(1,file='../CISCO/filter_K.dat', status='old')
  read(1,*)
10 read(1,*,end=20)lambda(i),response(i)
!!$10   read(1,*,end=20)i,lambda(i),response(i)
  i=i+1
  goto 10
20 continue
  close(1)
  num=i-1
!!$     num=i

  p_ira=0.
  t0=lambda(1)
  t1=lambda(num)
  call qsimp_SNIa(est_pII,t0,t1,p_ira)
  call qsimp_SNIa(est_pII2,t0,t1,p_ira2)
  write(*,*)p_ira2/p_ira,p_ira2,p_ira  
end program main

!!$    ======================== function est_pII =============
function est_pII(t)
  use global_var
  real est_pII,t

  call locate(lambda,num,t,k)
  est_pII=(response(k+1)-response(k))*(t-lambda(k))/&
       (lambda(k+1)-lambda(k))+response(k)

  return
end function est_pII
!!$    ======================== function est_pII =============
function est_pII2(t)
  use global_var
  real est_pII2,t

  call locate(lambda,num,t,k)
  est_pII2=(response(k+1)-response(k))*(t-lambda(k))/&
       (lambda(k+1)-lambda(k))+response(k)
  est_pII2=est_pII2*t

  return
end function est_pII2
!!$    === taken from Numerical Recipes
SUBROUTINE qsimp_SNIa(func,a,b,s)
  INTEGER JMAX
  REAL a,b,func,s,EPS
  EXTERNAL func
!!$     PARAMETER (EPS=1.e-1, JMAX=25)
  PARAMETER (EPS=1.e-4, JMAX=25)
!!$U    USES trapzd_SNIa
  INTEGER j
  REAL os,ost,st

  ost=-1.e30
  os= -1.e30
  do j=1,JMAX
     call trapzd_SNIa(func,a,b,st,j)
     s=(4.*st-ost)/3.
     if (abs(s-os).lt.EPS*abs(os)) return
     os=s
     ost=st
  enddo
  write(*,*) 'too many steps in qsimp_SNIa'
  stop
END SUBROUTINE qsimp_SNIa

SUBROUTINE trapzd_SNIa(func,a,b,s,n)
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
     do 11 j=1,it
        tot=tot+func(x)
        x=x+del
     enddo
     s=0.5*(s+(b-a)*tot/tnm)
  endif
  return
END SUBROUTINE trapzd_SNIa

subroutine locate(xx,n,x,j)
  integer j,n
  real x,xx(n)
  integer jl,jm,ju

  jl=0
  ju=n+1
10 if(ju-jl.gt.1)then
     jm=(ju+jl)/2
     if((xx(n).ge.xx(1)).eqv.(x.ge.xx(jm)))then
        jl=jm
     else
        ju=jm
     endif
     goto 10
  endif
  if(x.eq.xx(1))then
     j=1
  elseif(x.eq.xx(n))then
     j=n-1
  else
     j=jl
  endif
  return
end subroutine locate
