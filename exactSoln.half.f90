!!!------------------------------------------------------------------!!!
!!Module for solving the self-similar flow equation                   !!
!!   V_yyt = (a'/a)(V_yy +yV_yyy)+(1/Ra^2)V_yyyy -VV_yyy +V_yV_yy)    !!
!!(check notes). The self-similar flow is given by the equations      !!
!!        u=-(xb)V_y, v=aV, p=(x^2 b^2)P(t) +P_0 (y,t)                !!
!!The code is self contained. the solution is given by vectv, where   !!
!!vectv(0) is V(-1), vectv(m) is V(0). boundary conditions have to be !!
!!changed for vector f within advance_exact_dt if whole channel       !!
!!simulations are required, and also the value of hy. vector f        !!
!!corresponds to a transformation of V:                               !!
!!        V(y,t)=-(f(y,t)+a')/a.                                      !!
!!Contains subroutine penta for solving penta-diagonal linear systems !!
!!with standard gaussian elimination.                                 !!
!!                                                                    !!
!!Author: Leonardo Espin.   Modified: 1/14/2011                       !! 
!!!------------------------------------------------------------------!!!
module exact_solution
  implicit none
  private
  public :: vectv
  public :: initialize_exsolution,advance_dt_exact,oscillations,finalize_exsolution,penta,kick_start,recover_ss
  integer, parameter :: double=kind(1.0d0)
  real(double), parameter :: pi=3.1415926535897931
  real(double) :: hy
  real(double), dimension(:), allocatable  :: f,vectv,etaL
contains
  subroutine initialize_exsolution(m)
    integer :: m,i 
    hy=1.d0/dble(m)
    if(.not.allocated(vectv)) then
       allocate(f(0:m+2))
       allocate(vectv(0:m))
       allocate(etaL(m+1))
    endif
    f=0.d0
    vectv=0.d0
    etaL(1:m+1)=(/ (-1.d0+dble(i)*hy, i=0,m) /)
  end subroutine initialize_exsolution

  subroutine advance_dt_exact(m,t,dt,del1,del2,rm,phase)
    integer, intent(in) :: m
    integer :: n
    real(double), intent(in) :: t,dt,del1,del2,rm,phase
    real(double), dimension(m-1) :: a,b,c,d,e,abar,bbar,cbar,dbar,ebar,tnonlin
    real(double), dimension(m-1) :: r,x,trm1,trm1b
    real(double), dimension(3) :: sw
    real(double) :: uwp,trm2,trm2b
    real(double) :: fna,fnad,fnar,fnao,fnaod,fnaor
       
    n=m-1   !for matching # of interior points in the n-s code
    sw=oscillations(t,del1)
    fnao=sw(1)
    fnaod=sw(2) 
    fnaor=sw(3) 
    sw=oscillations(t+dt,del1)
    fna=sw(1) 
    fnad=sw(2) 
    fnar=sw(3) 
    trm1(1:n)=(1.d0+etaL(2:n+1))*dt*fnar/hy  
    trm2=dt/(hy*hy*rm*fna*fna)
    trm1b(1:n)=(1.d0+etaL(2:n+1))*dt*fnaor/hy 
    trm2b=dt/(hy*hy*rm*fnao*fnao)
    a(1:n)=.25d0*trm1(1:n)-.5d0*trm2
    b(1:n)=1.d0-dt*fnar-.5d0*trm1(1:n)+2.d0*trm2
    c(1:n)=-2.d0+2.d0*dt*fnar-3.d0*trm2
    d(1:n)=1.d0-dt*fnar+.5d0*trm1(1:n)+2.d0*trm2
    e(1:n)=-.25d0*trm1(1:n)-.5d0*trm2
    abar(1:n)=-.25d0*trm1b(1:n)+.5d0*trm2b
    bbar(1:n)=1.d0+dt*fnaor+.5d0*trm1b(1:n)-2.d0*trm2b
    cbar(1:n)=-2.d0-2.d0*dt*fnaor+3.d0*trm2b
    dbar(1:n)=1.d0+dt*fnaor-.5d0*trm1b(1:n)-2.d0*trm2b
    ebar(1:n)=.25d0*trm1b(1:n)+.5d0*trm2b
    tnonlin(1:n)=.5d0*f(2:n+1)*(f(4:n+3)-2.d0*f(3:n+2)+2.d0*f(1:n)-f(0:n-1))&
         -.5d0*(f(3:n+2)-f(1:n))*(f(3:n+2)-2.d0*f(2:n+1)+f(1:n))
    tnonlin(1:n)=dt*tnonlin(1:n)/fnao/hy
    r(1:n)=abar(1:n)*f(0:n-1)+bbar(1:n)*f(1:n)+cbar(1:n)*f(2:n+1)&
         +dbar(1:n)*f(3:n+2)+ebar(1:n)*f(4:n+3)+tnonlin(1:n)
    c(1)=c(1)+a(1)
    c(n)=c(n)-e(n)
    sw=oscillations(t+dt+phase,del2)
    uwp=2.d0*hy*fna*sw(3)!horizontal oscillation
    r(1)=r(1)+uwp*a(1)
    r(n-1)=r(n-1)+fnad*e(n-1)
    r(n)=r(n)+fnad*(d(n)+2.d0*e(n))
    call penta(x,a,b,c,d,e,r,n)
    f(2:n+1)=x(1:n)
    f(0)=f(2)-uwp
    f(1)=0.d0
    f(n+2)=-fnad
    f(n+3)=-f(n+1)-2.d0*fnad
    vectv(0:n+1)=-(f(1:n+2)+fnad)/fna
  end subroutine advance_dt_exact

  real(double) function kick_start(m,del1,del2,rm,phase) result(r)
    integer, intent(in) :: m
    real(double), intent(in) :: del1,del2,rm,phase
    real(double) :: dt 
    r=0.d0
    dt=2.d0*pi/60000.d0
    do 
       call advance_dt_exact(m,r,dt,del1,del2,rm,phase)
       r=r+dt
       if(r.gt.(300.d0*pi)) exit  
    enddo
  end function kick_start
  
  subroutine recover_ss(m,r,d1,d2,ph,exsol)
    integer, intent(in) :: m
    real(double), intent(in) :: r,d1,d2,ph
    real(double) :: fna,hy,uwp
    real(double), dimension(3) :: sw
    real(double), dimension(0:), intent(in) :: exsol 
    hy=1.d0/dble(m)
    sw=oscillations(r,d1)
    vectv(0:m)=exsol(0:m)
    f(1:m+1)=-sw(1)*exsol(0:m)-sw(2)
    f(m+2)=-f(m)-2.d0*sw(2)
    fna=sw(1)
    sw=oscillations(r+ph,d2)
    uwp=2.d0*hy*fna*sw(3) 
    f(0)=f(2)-uwp
  end subroutine recover_ss

  subroutine penta(x,a,b,c,d,e,h,n)
    integer, intent(in) :: n
    integer :: i
    real(double), dimension(:), intent(inout) :: x,a,b,c,d,e,h
    d(1)=d(1)/c(1)
    e(1)=e(1)/c(1)
    h(1)=h(1)/c(1)
    c(2)=c(2)-b(2)*d(1)
    d(2)=(d(2)-b(2)*e(1))/c(2)
    h(2)=(h(2)-b(2)*h(1))/c(2)
    e(2)=e(2)/c(2)
    do i=3,n
       b(i)=b(i)-a(i)*d(i-2)
       c(i)=c(i)-a(i)*e(i-2)
       h(i)=h(i)-a(i)*h(i-2)
       c(i)=c(i)-b(i)*d(i-1)
       d(i)=(d(i)-b(i)*e(i-1))/c(i)
       e(i)=e(i)/c(i)
       h(i)=(h(i)-b(i)*h(i-1))/c(i)
    enddo
    x(n)=h(n)
    x(n-1)=h(n-1)-d(n-1)*x(n)
    do i=n-2,1,-1
       x(i)=h(i)-d(i)*x(i+1)-e(i)*x(i+2)
    enddo
  end subroutine penta
  
  function oscillations(t,del)
     real(double), dimension(3) :: oscillations
     real(double) :: t,del
     oscillations(1)=1.d0+del*sin(t)  
     oscillations(2)=del*cos(t)
     oscillations(3)=oscillations(2)/oscillations(1)
  end function oscillations
  
  subroutine finalize_exsolution
    deallocate(f,vectv,etaL)
  end subroutine finalize_exsolution
end module exact_solution
