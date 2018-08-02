!!!------------------------------------------------------------------!!!
!!Computes Flow in a squeeze-stretch channel. code assumes viscous    !!
!!scaling of the Navier-Stokes equations and a staggered grid         !!
!!discretization (check notes)                                        !!
!!                                                                    !!
!!Main subroutine is advance_dt_ns which was modified so that dt is   !!
!!is computed only during first step of ADI and remains equal during  !!
!!the second step.                                                    !!
!!Requires files pSolver.f90 (the Poisson solver) and                 !!
!!exactSoln.half.f90                                                  !!
!!Code assumes half-channel computations. Have to adjust hy, vectors  !!
!!eta and xi in subroutine initialize_ns plus boundary condition in   !!
!!subroutine bcs (and module exact_solution as well) if whole-channel !!
!!simulations are required.                                           !!
!!Subroutines compute_FG and prbcs might have to be adjusted          !!
!!depending on the boundary conditions used with the system           !!
!!(e.g computations in periodic domain)                               !!
!!                                                                    !!
!!  *** contains new self-similar compatible boundary conditions ***  !!
!!                                                                    !!
!!Author: Leonardo Espin.   Modified: 3/23/2011                       !! 
!!!------------------------------------------------------------------!!!
module ns_solver
  use poisson_solver 
  use exact_solution
  implicit none
  private
  public :: u,v,p,t,dt 
  public :: initialize_ns,self_similar_ns,recover_ns,advance_dt_ns,finalize_ns 
  integer, parameter :: double=kind(1.0d0)
  real(double), parameter :: pi=3.1415926535897931
  integer :: nx,ny,i,j 
  real(double)  :: rm,hx,hy,gmma,eps,t,dt 
  real(double)  :: del1,del2,phase
  real(double)  :: fna,fnad,fnar,fnb,fnbd,fnbr
  real(double), dimension (3) :: sw    
  real(double), dimension (:,:), allocatable :: F,G,u,v,p,rhs  
  real(double), dimension (:), allocatable :: eta,xi 
contains
  subroutine initialize_ns(m,n,L,re,d1,d2,ph,gma)
    integer :: m,n
    real(double)  :: re,L,d1,d2,ph,gma
    ny=m
    nx=n
    eps=L
    rm=re
    phase=ph 
    gmma=gma
    del1=d1
    del2=d2
    hx=1.d0/dble(nx)       
    hy=1.d0/dble(ny)
    if(.not.allocated(u)) then
       allocate(u(0:ny+1,0:nx+1))
       allocate(v(0:ny+1,0:nx+1))
       allocate(p(0:ny+1,0:nx+1))
       allocate(rhs(0:ny+1,0:nx+1))
       allocate(F(0:ny+1,0:nx+1))
       allocate(G(0:ny+1,0:nx+1))
       allocate(eta(0:ny))
       allocate(xi(0:nx))
    end if
    p=0.d0
    u=0.d0
    v=0.d0
    F=0.d0
    G=0.d0
    rhs=0.d0
    eta(0:ny)=(/ (-1.d0+ dble(i)*hy, i=0,ny) /)
    xi(0:nx)=(/ (dble(i)*hx, i=0,nx) /)
    call initialize_hockney(ny)
  end subroutine initialize_ns

  subroutine self_similar_ns 
    sw=oscillations(t,del1) 
    fna=sw(1)
    fnad=sw(2)
    sw=oscillations(t+phase,del2)
    fnb=sw(1)
    fnbd=sw(2)
    do j=0,nx+1  !sets boundary conditions for v as well 
       v(0:ny,j)=fna*vectv(0:ny)
    enddo
    do j=0,nx
       u(1:ny,j)=-xi(j)*fnb*(vectv(1:ny)-vectv(0:ny-1))/hy
    enddo 
    !-----exterior points--------- 
    u(0,1:nx)=2.d0*fnbd*xi(1:nx)-u(1,1:nx) !check bcs subroutine to confirm right boundary contions!
    u(ny+1,1:nx)=u(ny,1:nx) 
  end subroutine self_similar_ns

  subroutine recover_ns(eu,ev,ep)
    real(double), dimension(0:,0:), intent(in) :: eu,ev,ep
    u(0:ny+1,0:nx+1)=eu(0:ny+1,0:nx+1)
    v(0:ny+1,0:nx+1)=ev(0:ny+1,0:nx+1)
    p(0:ny+1,0:nx+1)=ep(0:ny+1,0:nx+1)
  end subroutine recover_ns

  real(double) function compute_dt() result(r)
    real(double) :: umax,vmax,dtu,dtv,dtR   
    umax=0.0d0
    vmax=0.0d0
    dtu=10.0d0
    dtv=10.0d0
    dtR=0.5d0*rm/(1.d0/hx/hx +1.d0/hy/hy)
    umax=maxval(abs(u))
    vmax=maxval(abs(v))
    if(umax.gt.0.d0) dtu=hx/umax 
    if(vmax.gt.0.d0) dtv=hy/vmax 
    r=2.0d0*min(dtR,dtv,dtu)!0.3=fast 0.08=slow 0.01=vslow
  end function compute_dt

  subroutine advance_dt_ns(first_half)
    logical, intent(in) :: first_half
    if(first_half) then
       dt=compute_dt()
       dt=dt/2.d0
    end if
    call compute_FG
    call compute_rhs
    call prbcs(1.d0/eps/eps/fna/fna,1.d0/fnb/fnb)  
    call hockney(ny,nx,hy,hx,1.d0/eps/eps/fna/fna,1.d0/fnb/fnb,p,rhs) 
    call update_vel
    call bcs  
    t=t+dt
  end subroutine advance_dt_ns

  subroutine compute_FG  
    real(double)  :: g2v,du2dx,dv2dy,duvdy,duvdx,g2u
    sw=oscillations(t,del1)
    fna=sw(1)
    fnad=sw(2)
    fnar=sw(3)
    sw=oscillations(t+phase,del2)
    fnb=sw(1)
    fnbd=sw(2)
    fnbr=sw(3)
!$omp parallel num_threads(8) default(shared) private(du2dx,dv2dy,duvdy,duvdx,g2u,g2v)
!$omp do 
    do i=1,nx-1
       do j=1,ny
          du2dx=0.25d0*(((u(j,i)+u(j,i+1))/fnb -(xi(i+1)+xi(i))*fnbr)*(u(j,i)+u(j,i+1))&
               -((u(j,i-1)+u(j,i))/fnb -(xi(i)+xi(i-1))*fnbr)*(u(j,i-1)+u(j,i))&
               +gmma*(abs((u(j,i)+u(j,i+1))/fnb -(xi(i+1)+xi(i))*fnbr)*(u(j,i) -u(j,i+1))&
               -abs((u(j,i-1)+u(j,i))/fnb -(xi(i)+xi(i-1))*fnbr)*(u(j,i-1)-u(j,i))))/hx
          duvdy=0.25d0*(((v(j,i)+v(j,i+1))/fna -eta(j)*fnar)*(u(j,i)+u(j+1,i))&
               -((v(j-1,i)+v(j-1,i+1))/fna -eta(j-1)*fnar)*(u(j-1,i)+u(j,i)) &
               +gmma*(abs((v(j,i)+v(j,i+1))/fna -eta(j)*fnar)*(u(j,i)-u(j+1,i))&
               -abs((v(j-1,i)+v(j-1,i+1))/fna -eta(j-1)*fnar)*(u(j-1,i)-u(j,i))))/hy
          g2u=eps*eps*(u(j,i+1)-2.*u(j,i)+u(j,i-1))/fnb/fnb/hx/hx +(u(j+1,i)-2.*u(j,i)+u(j-1,i))/fna/fna/hy/hy
          F(j,i)=u(j,i)+dt*(g2u/rm -du2dx -duvdy-(fnar+fnbr)*u(j,i))
       enddo
    enddo
!$omp end do
!$omp single
    F(1:ny,0)=u(1:ny,0)   !BOUNDARY VALUE of u(1:ny,0) at tn+1
    F(1:ny,nx)=u(1:ny,nx)
!$omp end single
!$omp do 
    do i=1,nx
       do j=1,ny-1
          duvdx=0.25d0*(((u(j,i)+u(j+1,i))/fnb -xi(i)*fnbr)*(v(j,i)+v(j,i+1))&
               -((u(j,i-1)+u(j+1,i-1))/fnb -xi(i-1)*fnbr)*(v(j,i-1)+v(j,i))&
               +gmma*(abs((u(j,i)+u(j+1,i))/fnb -xi(i)*fnbr)*(v(j,i)-v(j,i+1))&
               -abs((u(j,i-1)+u(j+1,i-1))/fnb -xi(i-1)*fnbr)*(v(j,i-1)-v(j,i))))/hx
          dv2dy=0.25d0*(((v(j,i)+v(j+1,i))/fna -(eta(j+1)+eta(j))*fnar)*(v(j,i)+v(j+1,i))&
               -((v(j-1,i)+v(j,i))/fna -(eta(j)+eta(j-1))*fnar)*(v(j-1,i)+v(j,i))&
               +gmma*(abs((v(j,i)+v(j+1,i))/fna -(eta(j+1)+eta(j))*fnar)*(v(j,i)-v(j+1,i))&
               -abs((v(j-1,i)+v(j,i))/fna -(eta(j)+eta(j-1))*fnar)*(v(j-1,i)-v(j,i))))/hy
          g2v=eps*eps*(v(j,i+1)-2.*v(j,i)+v(j,i-1))/fnb/fnb/hx/hx +(v(j+1,i)-2.*v(j,i)+v(j-1,i))/fna/fna/hy/hy
          G(j,i)=v(j,i)+dt*(g2v/rm -duvdx -dv2dy-(fnar+fnbr)*v(j,i))
       enddo
    enddo
!$omp end do
!$omp end parallel
    G(0,1:nx)=v(0,1:nx)   !BOUNDARY VALUE of v(0,1:nx) at tn+1
    G(ny,1:nx)=v(ny,1:nx)
  end subroutine compute_FG

  subroutine compute_rhs
    sw=oscillations(t+dt,del1)
    fna=sw(1)
    fnad=sw(2)
    fnar=sw(3)
    sw=oscillations(t+dt+phase,del2)
    fnb=sw(1)
    fnbd=sw(2)
    fnbr=sw(3)
!$omp parallel do num_threads(8) default(shared) 
    do i=1,nx
       do j=1,ny
          rhs(j,i)=rm*((F(j,i)-F(j,i-1))/fnb/hx +(G(j,i)-G(j-1,i))/fna/hy)/dt!boundary values of F and G required
       enddo
    enddo
!$omp end parallel do
  end subroutine compute_rhs

  subroutine prbcs(cy,cx)
    !relaxation of pressure BC's
    !modifies the rhs according to dirichlet conditions for Poisson solver
    real(double)  cy,cx
    p(1:ny,0)=p(1:ny,1)
    p(1:ny,nx+1)=p(1:ny,nx)
    p(0,1:nx)=p(1,1:nx)
    p(ny+1,1:nx)=p(ny,1:nx)
    rhs(1:ny,1)=rhs(1:ny,1)-cx*p(1:ny,0)/hx/hx
    rhs(1:ny,nx)=rhs(1:ny,nx)-cx*p(1:ny,nx+1)/hx/hx
    rhs(1,1:nx)=rhs(1,1:nx)-cy*p(0,1:nx)/hy/hy
    rhs(ny,1:nx)=rhs(ny,1:nx)-cy*p(ny+1,1:nx)/hy/hy
  end subroutine prbcs

  subroutine update_vel
    !a,a' etc were computed at compute_rhs and are at t_{n+1}
    do j=1,ny
       u(j,1:nx-1)=F(j,1:nx-1)-dt*(p(j,2:nx)-p(j,1:nx-1))/fnb/hx/rm!boundary values of pressure required
    enddo
    do i=1,nx
       v(1:ny-1,i)=G(1:ny-1,i)-dt*(p(2:ny,i)-p(1:ny-1,i))/fna/hy/rm/eps/eps
    enddo
  end subroutine update_vel

  subroutine bcs
    !fna, etc were computed at t+dt in compute_rhs() above
    call advance_dt_exact(ny,t,dt,del1,del2,rm,phase)
    !----boundary values--------
    u(1:ny,0)=0.d0 
    u(1:ny,nx)=2.d0*u(1:ny,nx-1)-u(1:ny,nx-2)!-xi(nx)*fnb*(vectv(1:ny)-vectv(0:ny-1))/hy 
    v(0,1:nx)=-fnad
    v(ny,1:nx)=0.d0!fnad
    !-----exterior points---------
    u(0,1:nx)=2.d0*fnbd*xi(1:nx)-u(1,1:nx) 
    u(ny+1,1:nx)=u(ny,1:nx)
    v(1:ny,0)=v(1:ny,1)
    v(1:ny,nx+1)=v(1:ny,nx)!2.d0*fna*vectv(1:ny)-v(1:ny,nx)
  end subroutine bcs

  subroutine finalize_ns
    deallocate(u,v,p)
    deallocate(F,G,rhs)
    deallocate(eta,xi)
    call finalize_hockney
  end subroutine finalize_ns
end module ns_solver
