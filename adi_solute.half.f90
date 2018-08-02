!!!------------------------------------------------------------------!!!
!!Module for solving the oxygen transport and consumption model       !!
!!     C_t +(1/b)(u-xb')C_x +(1/a)(v-ya')C_y =                        !!
!!          (1/Pe){(del^2/b^2)C_xx +(1/a^2)C_yy}                      !!
!!     O_t -x(b'/b)O_x -y(a'/a)O_y=                                   !!
!!          (1/Pe){(del^2/b^2)O_xx +(1/a^2)O_yy -(lambda)O}           !!
!!C is the channel concentration, O (theta) the concentration in the  !!
!!channel walls.                                                      !!
!!                                                                    !!
!!Requires the files exactSoln.half.f90 plust the files nsSolver.f90  !!
!!pSolver.f90 if (u,v) are computed with the Navier-Stokes equations  !!
!!instead of assuming a self-similar velocity field.                  !!
!!Contains subroutine penta for solving tri-diagonal linear systems   !!
!!with standard gaussian elimination.                                 !!
!!                                                                    !!
!!Author: Leonardo Espin.   Modified: 1/14/2011                       !!
!!!------------------------------------------------------------------!!!
module adi_solute
  use ns_solver
  use exact_solution
  implicit none
  private
  public :: s
  public :: initialize_adi,finalize_adi,tri,advance_dt_adi_ss,advance_dt_adi
  integer, parameter :: double=kind(1.0d0)
  real(double), parameter :: pi=3.1415926535897931
  integer :: nx,ny,i,j
  real(double)  :: hx,hy,rm,lambda,Pe,delta,kap,cta,ctb,rate
  real(double)  :: del1,del2,phase
  real(double)  :: fna,fnad,fnar,fnb,fnbd,fnbr,fnao,fnado,fnaro,fnbo,fnbdo,fnbro
  real(double), dimension (3) :: sw
  real(double), dimension (:,:), allocatable :: s,stmp,vold,uold  
  real(double), dimension (:), allocatable :: eta,xi,exold
contains
  subroutine initialize_adi(m,n,del,re,d1,d2,ph,lamb,peclet,kappa)
    integer, intent(in) :: m,n
    real(double), intent(in) :: re,d1,d2,ph,lamb,peclet,kappa,del
    ny=m
    nx=n
    rm=re
    delta=del
    phase=ph 
    del1=d1
    del2=d2
    lambda=lamb
    Pe=peclet
    kap=kappa
    hx=1.d0/dble(nx)       
    hy=1.d0/dble(ny)
    rate=sqrt(lambda)/(1+kap*sqrt(lambda))
    if(.not.allocated(s)) then
       allocate(s(0:3*ny+3,0:nx+1))
       allocate(stmp(0:3*ny+3,0:nx+1))
       allocate(vold(0:ny+1,0:nx+1))
       allocate(uold(0:ny+1,0:nx+1))
       allocate(exold(0:ny))
       allocate(eta(0:3*ny+3))
       allocate(xi(0:nx+1))
    end if
    s=0.d0
!!$    do i=2*ny+2,3*ny+3
!!$       s(i,0:nx+1)=1.d0
!!$    enddo
    stmp=0.d0
    xi(0:nx+1)=(/ (dble(i-0.5)*hx, i=0,nx+1) /)
    eta(0:2*ny+1)=(/ (-3.d0+ dble(i-0.5)*hy, i=0,2*ny+1) /)
    eta(2*ny+2:3*ny+3)=(/ (-1.d0+ dble(i-2*ny-2-0.5)*hy, i=2*ny+2,3*ny+3) /)
  end subroutine initialize_adi

  subroutine advance_dt_adi(r)
    real(double), intent(in) :: r
    real(double) :: c_t                 !fortran passes arguments by reference. need current 
    c_t=r                               !value of t, since t changes after advance_dt_ns
    vold(0:ny+1,0:nx+1)=v(0:ny+1,0:nx+1)  
    call advance_dt_ns(.true.)          !this line is not in a parallel region to avoid nesting
    call adi_i(2.d0*dt,c_t)             !after previous line t=c_t+dt/2, c_t stores previous value of time
    uold(0:ny+1,0:nx+1)=u(0:ny+1,0:nx+1)!2.d0*dt is necessary since advance_dt_ns divides dt by 2
    call advance_dt_ns(.false.)
    call adi_ii(2.d0*dt,c_t)            !after previous line t=c_t+dt/2+dt/2
  end subroutine advance_dt_adi

  subroutine advance_dt_adi_ss(dr,r)
    real(double), intent(in) :: r,dr
    call adi_i_ss(dr,r)  !assuming self-similar velocities
    call adi_ii_ss(dr,r) 
  end subroutine advance_dt_adi_ss
     
  subroutine adi_i(dr,r)
    real(double), intent(in) :: r,dr
    real(double) :: sfp,sfm
    real(double) :: ai(nx),bi(nx),ci(nx),ri(nx),x1(nx),gpi(nx),gmi(nx),fpi(nx),fmi(nx)
    sw=oscillations(r+dr/2.d0,del1)
    fna=sw(1)
    fnad=sw(2)
    fnar=sw(3)
    sw=oscillations(r+dr/2.d0,del2)
    fnb=sw(1)
    fnbd=sw(2)
    fnbr=sw(3)
    sw=oscillations(r,del1)
    fnao=sw(1)
    fnado=sw(2)
    fnaro=sw(3)
    sw=oscillations(r,del2)
    fnbo=sw(1)
    fnbdo=sw(2)
    fnbro=sw(3)
    cta=0.5d0*dr*delta*delta/Pe/fnb/fnb/hx/hx
    ctb=0.5d0*dr/Pe/fnao/fnao/hy/hy
    gpi(1:nx)=0.5d0*dr*xi(1:nx)*max(fnbr,0.d0)/hx
    gmi(1:nx)=0.5d0*dr*xi(1:nx)*min(fnbr,0.d0)/hx 
!$omp parallel num_threads(8) default(shared) private(sfp,sfm,ai,bi,ci,ri,fpi,fmi,x1) firstprivate(gpi,gmi)
!$omp do 
    do i=1,2*ny
       sfp=0.5d0*dr*max(eta(i)*fnaro,0.d0)/hy
       sfm=0.5d0*dr*min(eta(i)*fnaro,0.d0)/hy
       ai(1:nx)=gmi(1:nx) -cta
       bi(1:nx)=1.d0 +gpi(1:nx)-gmi(1:nx) +2.d0*cta +dr*lambda/Pe
       ci(1:nx)=-gpi(1:nx) -cta
       ri(1:nx)=(-sfm+ctb)*s(i-1,1:nx)
       ri(1:nx)=ri(1:nx)+(1.d0 -sfp+sfm-2.d0*ctb)*s(i,1:nx)
       ri(1:nx)=ri(1:nx)+(sfp+ctb)*s(i+1,1:nx)
       bi(1)=bi(1)+ai(1)    
       bi(nx)=bi(nx)+ci(nx)
       call tri(x1,ai,bi,ci,ri,nx)
       stmp(i,1:nx)=x1(1:nx)
       stmp(i,0)=stmp(i,1)
       stmp(i,nx+1)=stmp(i,nx)
    end do
!$omp end do
!-------------------------solute concentration inside the channel----------------------------- 
!$omp do 
    do i=2*ny+3,3*ny+2       
       fpi(1:nx)=0.5d0*dr*max(eta(i)*fnaro-0.5d0*(vold(i-2*ny-2,1:nx)+vold(i-2*ny-3,1:nx))/fnao,0.d0)/hy
       fmi(1:nx)=0.5d0*dr*min(eta(i)*fnaro-0.5d0*(vold(i-2*ny-2,1:nx)+vold(i-2*ny-3,1:nx))/fnao,0.d0)/hy 
       gpi(1:nx)=0.5d0*dr*max(fnbr*xi(1:nx)-0.5d0*(u(i-2*ny-2,1:nx)+u(i-2*ny-2,0:nx-1))/fnb,0.d0)/hx
       gmi(1:nx)=0.5d0*dr*min(fnbr*xi(1:nx)-0.5d0*(u(i-2*ny-2,1:nx)+u(i-2*ny-2,0:nx-1))/fnb,0.d0)/hx 
       ai(1:nx)=gmi(1:nx) -cta
       bi(1:nx)=1.d0 +gpi(1:nx)-gmi(1:nx) +2.d0*cta 
       ci(1:nx)=-gpi(1:nx) -cta
       ri(1:nx)=(-fmi(1:nx)+ctb)*s(i-1,1:nx)
       ri(1:nx)=ri(1:nx)+(1.d0 -fpi(1:nx)+fmi(1:nx)-2.d0*ctb)*s(i,1:nx)
       ri(1:nx)=ri(1:nx)+(fpi(1:nx)+ctb)*s(i+1,1:nx) 
       bi(1)=bi(1)+ai(1)                          
       bi(nx)=bi(nx)-ci(nx)               
       ri(nx)=ri(nx)-ci(nx)*2.d0
       call tri(x1,ai,bi,ci,ri,nx)
       stmp(i,1:nx)=x1(1:nx)
       stmp(i,nx+1)=2.d0-stmp(i,nx)
       stmp(i,0)=stmp(i,1)
    end do
!$omp end do
!$omp end parallel
  end subroutine adi_i

  subroutine adi_ii(dr,r)
    real(double), intent(in) :: r,dr
    real(double) :: aii(3*ny+2),bii(3*ny+2),cii(3*ny+2),dii(3*ny+2),eii(3*ny+2),rii(3*ny+2)
    real(double) :: x2(3*ny+2),gpii(3*ny+2),gmii(3*ny+2),fpii(3*ny+2),fmii(3*ny+2)
    sw=oscillations(r+dr,del1)
    fna=sw(1)
    fnad=sw(2)
    fnar=sw(3)
    sw=oscillations(r+dr,del2)
    fnb=sw(1)
    fnbd=sw(2)
    fnbr=sw(3)
    sw=oscillations(r+dr/2.d0,del1)
    fnao=sw(1)
    fnado=sw(2)
    fnaro=sw(3)
    sw=oscillations(r+dr/2.d0,del2)
    fnbo=sw(1)
    fnbdo=sw(2)
    fnbro=sw(3)
    cta=0.5d0*dr/Pe/fna/fna/hy/hy
    ctb=0.5d0*dr*delta*delta/Pe/fnbo/fnbo/hx/hx
    fpii(1:3*ny+2)=0.5d0*dr*max(eta(1:3*ny+2)*fnar,0.d0)/hy  
    fmii(1:3*ny+2)=0.5d0*dr*min(eta(1:3*ny+2)*fnar,0.d0)/hy
!$omp parallel do num_threads(8) default(shared) private(aii,bii,cii,dii,eii,rii,gpii,gmii,x2) firstprivate(fpii,fmii)
    do j=1,nx
       fpii(2*ny+3:3*ny+2)=0.5d0*dr*max(eta(2*ny+3:3*ny+2)*fnar -0.5d0*(v(1:ny,j)+v(0:ny-1,j))/fna,0.d0)/hy
       fmii(2*ny+3:3*ny+2)=0.5d0*dr*min(eta(2*ny+3:3*ny+2)*fnar -0.5d0*(v(1:ny,j)+v(0:ny-1,j))/fna,0.d0)/hy 
       gpii(1:2*ny+2)=0.5d0*dr*xi(j)*max(fnbro,0.d0)/hx    
       gmii(1:2*ny+2)=0.5d0*dr*xi(j)*min(fnbro,0.d0)/hx 
       gpii(2*ny+3:3*ny+2)=0.5d0*dr*max(fnbro*xi(j)-0.5d0*(uold(1:ny,j)+uold(1:ny,j-1))/fnbo,0.d0)/hx
       gmii(2*ny+3:3*ny+2)=0.5d0*dr*min(fnbro*xi(j)-0.5d0*(uold(1:ny,j)+uold(1:ny,j-1))/fnbo,0.d0)/hx 
       aii(1:3*ny+2)=0.d0
       bii(1:3*ny+2)=fmii(1:3*ny+2)-cta
       cii(1:3*ny+2)=1.d0 +fpii(1:3*ny+2)-fmii(1:3*ny+2) +2.d0*cta
       dii(1:3*ny+2)=-fpii(1:3*ny+2)-cta
       eii(1:3*ny+2)=0.d0
       rii(1:3*ny+2)=(-gmii(1:3*ny+2)+ctb)*stmp(1:3*ny+2,j-1)
       rii(1:3*ny+2)=rii(1:3*ny+2)+(1.d0 -gpii(1:3*ny+2)+gmii(1:3*ny+2)-2.d0*ctb)*stmp(1:3*ny+2,j)
       rii(1:3*ny+2)=rii(1:3*ny+2)+(gpii(1:3*ny+2)+ctb)*stmp(1:3*ny+2,j+1)
       bii(2*ny+1)=0.5d0*fna/kap -1.d0/hy 
       cii(2*ny+1)=0.5d0*fna/kap +1.d0/hy 
       dii(2*ny+1)=-0.5d0*fna/kap 
       eii(2*ny+1)=-0.5d0*fna/kap 
       rii(2*ny+1)=0.d0
       aii(2*ny+2)=-0.5d0*fna/kap 
       bii(2*ny+2)=-0.5d0*fna/kap
       cii(2*ny+2)=0.5d0*fna/kap +1.d0/hy
       dii(2*ny+2)=0.5d0*fna/kap -1.d0/hy 
       rii(2*ny+2)=0.d0
       cii(3*ny+2)=cii(3*ny+2)+dii(3*ny+2)  
       !cii(1)=cii(1)+bii(1)*(2.d0-rate*hy)/(2.d0+rate*hy)               !exp decreasing
       !cii(1)=cii(1)+bii(1)                 !h_newmann
       cii(1)=cii(1)-bii(1)                 !0.d0 dirchlet condition
       call penta(x2,aii,bii,cii,dii,eii,rii,3*ny+2)
       s(1:3*ny+2,j)=x2(1:3*ny+2)
       s(3*ny+3,j)=s(3*ny+2,j) 
       !s(0,j)=(2.d0-rate*hy)*s(1,j)/(2.d0+rate*hy)                !exp decreasing
       !s(0,j)=s(1,j)                  !h_newmann       
       s(0,j)=-s(1,j)                 !0.d0 dirchlet cond  
    end do
!$omp end parallel do
    s(1:3*ny+2,0)=s(1:3*ny+2,1)                   !how the interfacial values 
    s(1:3*ny+2,nx+1)=s(1:3*ny+2,nx)               !are filled is important for
    s(2*ny+2:3*ny+3,nx+1)=2.d0-s(2*ny+2:3*ny+3,nx)!nusselt # computation!!
  end subroutine adi_ii

  subroutine adi_i_ss(dr,r)
    real(double), intent(in) :: r,dr
    real(double) :: ai(nx),bi(nx),ci(nx),ri(nx),x1(nx),gpi(nx),gmi(nx),sfp,sfm
!$omp parallel num_threads(4) default(shared) private(ai,bi,ci,ri,gpi,gmi,sfp,sfm,x1)
!$omp single
    sw=oscillations(r+dr/2.d0,del1)
    fna=sw(1)
    fnad=sw(2)
    fnar=sw(3)
    sw=oscillations(r+dr/2.d0,del2)
    fnb=sw(1)
    fnbd=sw(2)
    fnbr=sw(3)
    sw=oscillations(r,del1)
    fnao=sw(1)
    fnado=sw(2)
    fnaro=sw(3)
    sw=oscillations(r,del2)
    fnbo=sw(1)
    fnbdo=sw(2)
    fnbro=sw(3)
    cta=0.5d0*dr*delta*delta/Pe/fnb/fnb/hx/hx
    ctb=0.5d0*dr/Pe/fnao/fnao/hy/hy
!$omp end single
    gpi(1:nx)=0.5d0*dr*xi(1:nx)*max(fnbr,0.d0)/hx !these are not shared variables, 
    gmi(1:nx)=0.5d0*dr*xi(1:nx)*min(fnbr,0.d0)/hx !line has to be executed by each processor.
!$omp do 
    do i=1,2*ny
       sfp=0.5d0*dr*max(eta(i)*fnaro,0.d0)/hy
       sfm=0.5d0*dr*min(eta(i)*fnaro,0.d0)/hy
       ai(1:nx)=gmi(1:nx) -cta
       bi(1:nx)=1.d0 +gpi(1:nx)-gmi(1:nx) +2.d0*cta +dr*lambda/Pe
       ci(1:nx)=-gpi(1:nx) -cta
       ri(1:nx)=(-sfm+ctb)*s(i-1,1:nx)
       ri(1:nx)=ri(1:nx)+(1.d0 -sfp+sfm-2.d0*ctb)*s(i,1:nx)
       ri(1:nx)=ri(1:nx)+(sfp+ctb)*s(i+1,1:nx)
       bi(1)=bi(1)+ai(1)    
       bi(nx)=bi(nx)+ci(nx)
       call tri(x1,ai,bi,ci,ri,nx)
       stmp(i,1:nx)=x1(1:nx)
       stmp(i,0)=stmp(i,1)
       stmp(i,nx+1)=stmp(i,nx)
    end do
!$omp end do
!-------------------------solute concentration inside the channel----------------------------- 
!$omp single
    exold(0:ny)=vectv(0:ny)
    call advance_dt_exact(ny,r,dr/2.d0,del1,del2,rm,phase)!computes V(eta,t+dt/2)
!$omp end single
!$omp do 
    do i=2*ny+3,3*ny+2       
       sfp=0.5d0*dr*max(eta(i)*fnaro-0.5d0*(exold(i-2*ny-2)+exold(i-2*ny-3)),0.d0)/hy
       sfm=0.5d0*dr*min(eta(i)*fnaro-0.5d0*(exold(i-2*ny-2)+exold(i-2*ny-3)),0.d0)/hy
       gpi(1:nx)=0.5d0*dr*xi(1:nx)*max(fnbr+(vectv(i-2*ny-2)-vectv(i-2*ny-3))/hy,0.d0)/hx
       gmi(1:nx)=0.5d0*dr*xi(1:nx)*min(fnbr+(vectv(i-2*ny-2)-vectv(i-2*ny-3))/hy,0.d0)/hx 
       ai(1:nx)=gmi(1:nx) -cta
       bi(1:nx)=1.d0 +gpi(1:nx)-gmi(1:nx) +2.d0*cta 
       ci(1:nx)=-gpi(1:nx) -cta
       ri(1:nx)=(-sfm+ctb)*s(i-1,1:nx)
       ri(1:nx)=ri(1:nx)+(1.d0 -sfp+sfm-2.d0*ctb)*s(i,1:nx)
       ri(1:nx)=ri(1:nx)+(sfp+ctb)*s(i+1,1:nx) 
       bi(1)=bi(1)+ai(1)                          
       bi(nx)=bi(nx)-ci(nx)               
       ri(nx)=ri(nx)-ci(nx)*2.d0
       call tri(x1,ai,bi,ci,ri,nx)
       stmp(i,1:nx)=x1(1:nx)
       stmp(i,nx+1)=2.d0-stmp(i,nx)
       stmp(i,0)=stmp(i,1)
    end do
!$omp end do
!$omp end parallel
  end subroutine adi_i_ss

  subroutine adi_ii_ss(dr,r)
    real(double), intent(in) :: r,dr
    real(double) :: aii(3*ny+2),bii(3*ny+2),cii(3*ny+2),dii(3*ny+2),eii(3*ny+2),rii(3*ny+2)
    real(double) :: x2(3*ny+2),gpii(3*ny+2),gmii(3*ny+2),fpii(3*ny+2),fmii(3*ny+2)
!$omp parallel num_threads(4)  default(shared) private(aii,bii,cii,dii,eii,rii,gpii,gmii,x2)
!$omp single
    sw=oscillations(r+dr,del1)
    fna=sw(1)
    fnad=sw(2)
    fnar=sw(3)
    sw=oscillations(r+dr,del2)
    fnb=sw(1)
    fnbd=sw(2)
    fnbr=sw(3)
    sw=oscillations(r+dr/2.d0,del1)
    fnao=sw(1)
    fnado=sw(2)
    fnaro=sw(3)
    sw=oscillations(r+dr/2.d0,del2)
    fnbo=sw(1)
    fnbdo=sw(2)
    fnbro=sw(3)
    cta=0.5d0*dr/Pe/fna/fna/hy/hy
    ctb=0.5d0*dr*delta*delta/Pe/fnbo/fnbo/hx/hx
    exold(0:ny)=vectv(0:ny)
    call advance_dt_exact(ny,r+dr/2.d0,dr/2.d0,del1,del2,rm,phase)!computes V(eta,t+dt)
    fpii(1:2*ny+2)=0.5d0*dr*max(eta(1:2*ny+2)*fnar,0.d0)/hy  
    fmii(1:2*ny+2)=0.5d0*dr*min(eta(1:2*ny+2)*fnar,0.d0)/hy
    fpii(2*ny+3:3*ny+2)=0.5d0*dr*max(eta(2*ny+3:3*ny+2)*fnar -0.5d0*(vectv(1:ny)+vectv(0:ny-1)),0.d0)/hy
    fmii(2*ny+3:3*ny+2)=0.5d0*dr*min(eta(2*ny+3:3*ny+2)*fnar -0.5d0*(vectv(1:ny)+vectv(0:ny-1)),0.d0)/hy 
!$omp end single
!$omp do 
    do j=1,nx
       gpii(1:2*ny+2)=0.5d0*dr*xi(j)*max(fnbro,0.d0)/hx    
       gmii(1:2*ny+2)=0.5d0*dr*xi(j)*min(fnbro,0.d0)/hx 
       gpii(2*ny+3:3*ny+2)=0.5d0*dr*xi(j)*max(fnbro+(exold(1:ny)-exold(0:ny-1))/hy,0.d0)/hx
       gmii(2*ny+3:3*ny+2)=0.5d0*dr*xi(j)*min(fnbro+(exold(1:ny)-exold(0:ny-1))/hy,0.d0)/hx 
       aii(1:3*ny+2)=0.d0
       bii(1:3*ny+2)=fmii(1:3*ny+2)-cta
       cii(1:3*ny+2)=1.d0 +fpii(1:3*ny+2)-fmii(1:3*ny+2) +2.d0*cta
       dii(1:3*ny+2)=-fpii(1:3*ny+2)-cta
       eii(1:3*ny+2)=0.d0
       rii(1:3*ny+2)=(-gmii(1:3*ny+2)+ctb)*stmp(1:3*ny+2,j-1)
       rii(1:3*ny+2)=rii(1:3*ny+2)+(1.d0 -gpii(1:3*ny+2)+gmii(1:3*ny+2)-2.d0*ctb)*stmp(1:3*ny+2,j)
       rii(1:3*ny+2)=rii(1:3*ny+2)+(gpii(1:3*ny+2)+ctb)*stmp(1:3*ny+2,j+1)
       bii(2*ny+1)=0.5d0*fna/kap -1.d0/hy 
       cii(2*ny+1)=0.5d0*fna/kap +1.d0/hy 
       dii(2*ny+1)=-0.5d0*fna/kap 
       eii(2*ny+1)=-0.5d0*fna/kap 
       rii(2*ny+1)=0.d0
       aii(2*ny+2)=-0.5d0*fna/kap 
       bii(2*ny+2)=-0.5d0*fna/kap
       cii(2*ny+2)=0.5d0*fna/kap +1.d0/hy
       dii(2*ny+2)=0.5d0*fna/kap -1.d0/hy 
       rii(2*ny+2)=0.d0
       cii(3*ny+2)=cii(3*ny+2)+dii(3*ny+2)  
       !cii(1)=cii(1)+bii(1)*(2.d0-rate*hy)/(2.d0+rate*hy)               !exp decreasing
       !cii(1)=cii(1)+bii(1)                 !h_newmann
       cii(1)=cii(1)-bii(1)                 !0.d0 dirchlet condition
       call penta(x2,aii,bii,cii,dii,eii,rii,3*ny+2)
       s(1:3*ny+2,j)=x2(1:3*ny+2)
       s(3*ny+3,j)=s(3*ny+2,j) 
       !s(0,j)=(2.d0-rate*hy)*s(1,j)/(2.d0+rate*hy)                !exp decreasing
       !s(0,j)=s(1,j)                  !h_newmann       
       s(0,j)=-s(1,j)                 !0.d0 dirchlet cond  
    end do
!$omp end do
!$omp end parallel
    s(1:3*ny+2,0)=s(1:3*ny+2,1)                   !how the interfacial values 
    s(1:3*ny+2,nx+1)=s(1:3*ny+2,nx)               !are filled is important for
    s(2*ny+2:3*ny+3,nx+1)=2.d0-s(2*ny+2:3*ny+3,nx)!nusselt # computation!!
  end subroutine adi_ii_ss

  subroutine tri(x,a,b,c,h,n)
    integer, intent(in) :: n
    integer :: i
    real(double), dimension(:), intent(inout) :: x,a,b,c,h
    c(1)=c(1)/b(1)
    h(1)=h(1)/b(1)
    do i=2,n
       b(i)=b(i)-a(i)*c(i-1)
       c(i)=c(i)/b(i)
       h(i)=(h(i)-a(i)*h(i-1))/b(i)
    enddo
    x(n)=h(n)
    do i=n-1,1,-1
       x(i)=h(i)-c(i)*x(i+1)
    enddo
  end subroutine tri

  subroutine finalize_adi
    deallocate(s,stmp,uold,vold)
    deallocate(xi,eta,exold)
  end subroutine finalize_adi
end module adi_solute
