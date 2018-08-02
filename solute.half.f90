!!!------------------------------------------------------------------!!!
!!    ****   this code is for HALF-CHANNEL computations   ****        !!
!!                                                                    !!
!!Takes care of storing temporary results, making movie frames and    !!
!!computing quantities related to the oxygen transport and            !!
!!consumption model.                                                  !!
!!                                                                    !!
!!Parameters required for computation are read from input file.       !!
!!Main solute code is in file adi_solute.half.f90, which contains the !!
!!adi solute solver. It also requires the files exactSoln.half.f90    !!
!!plus the files nsSolver.f90 and pSolver.f90 for when velocity       !!
!!field  (u,v) is computed with the Navier-Stokes equations instead   !! 
!!of assuming a self-similar velocity field.                          !!
!!                                                                    !!
!!Choice between self-similar computations/Navier-Stokes computations !!
!!is controlled by parameter NAVIER.                                  !!
!!                                                                    !!
!!Whole channel computations require changing variable hy from        !!
!!  hy=1.d0/dble(ny) to hy=2.d0/dble(ny) in subroutines read_vars,    !! 
!!initialize_adi, initialize_ns, initialize_exsolution and            !!
!!recover_ss. This last subroutine is in file exactSoln.half.f90.     !!
!!Also boundary conditions have to be adjusted in subroutines         !!
!!adi_i, adi_ii (or adi_i_ss, adi_ii_ss if self-similar computations  !!
!!are being done),  subroutine bcs in file nsSolver.f90 (see also     !! 
!!the initial comments in that module) and subroutines                !!
!!advance_dt_exact and recover_ss.                                    !!
!!                                                                    !!
!!Author: Leonardo Espin.   Modified: 1/14/2011                       !!
!!!------------------------------------------------------------------!!!
program solute
  use ns_solver 
  use exact_solution, only: vectv,oscillations,initialize_exsolution,finalize_exsolution,kick_start,recover_ss
  use adi_solute
  implicit none
  integer, parameter :: double=kind(1.0d0)
  real(double), parameter :: pi=3.1415926535897931
  logical, parameter :: navier=.true.
  logical :: integrating 
  integer :: nx,ny,iteration,length,j,frame,next_frame
  real(double) :: hx,hy,t_a,tmax,d1,d2,next_time,phase
  !real(double) :: t,dt   !comment if navier-stokes solver is used. these variables are provided by the solver. 
  real(double), dimension (3) :: sw
  real(double), allocatable :: vectf(:),vectfo(:),ave(:),intgr(:)
  
  call read_vars
 !----------------------main loop----------------------------------------
  iteration=1
  do
     if(navier) then
        call advance_dt_adi(t)
     else
        call advance_dt_adi_ss(dt,t) !if navier-stokes solver is used, that library takes care of t and dt. 
        t=t+dt                       !these variables have to be controled by main program otherwise.
     endif
     if (t.gt.tmax) exit
     if(integrating) then
        do j=0,nx+1
           intgr(0:ny+1)=s(2*ny+2:3*ny+3,j)
           vectf(j)=spatial_average(hy,ny,intgr)
        enddo
        sw=oscillations(t,d1)
        vectf(0:nx+1)=sw(1)*vectf(0:nx+1)
        if(navier) then
           call time_integrate(2.d0*dt,nx,vectfo,vectf,ave) !recall that advance_dt_ns divides dt by 2?
        else
           call time_integrate(dt,nx,vectfo,vectf,ave)
        endif
        vectfo(0:nx+1)=vectf(0:nx+1)
        if(navier) then
           if(t.ge.next_time) then
              call save_results
              next_time=make_frame2(t,frame)
           endif
        else
           if(iteration.eq.next_frame) then
              call save_results
              next_frame=make_frame(iteration,frame)
           endif
        endif
     endif
     if((t.ge.t_a).and.(.not.(integrating))) then !start computing averages over a time period
        do j=0,nx+1
           intgr(0:ny+1)=s(2*ny+2:3*ny+3,j)
           vectfo(j)=spatial_average(hy,ny,intgr)
        enddo
        sw=oscillations(t,d1)
        vectfo(0:nx+1)=sw(1)*vectfo(0:nx+1)
        integrating=.true.
        call save_results
        if(navier) then
           next_time=make_frame2(t,frame)
        else
           next_frame=make_frame(iteration,frame)
        endif
     endif
     if(mod(iteration,1000).eq.0) then
        call compute_nusselt
        if(navier) call write_trace
     endif
     if((mod(iteration,50000).eq.0).and.(t.lt.t_a))then  
        call save_results
        if(navier) call save_velocities
     endif
     iteration=iteration+1
  enddo
!--------end main loop--------------------------------------------
call save_results
call save_averages(nx,ave)!used for saving time and y-averaged channel concentration. 
call compute_nusselt
if(navier)then
   if(t.ge.next_time) next_time=make_frame2(t,frame)
   call save_velocities
else
   if(iteration.eq.next_frame) next_frame=make_frame(iteration,frame)
endif
call finalize 
stop

contains
  subroutine read_vars    
    integer :: icont
    real(double)  :: sigma,gmma,Pec
    real(double) :: lambda,re,del,kappa
    open(1,file='input',status='old')
    read(1,*) nx
    read(1,*) ny
    read(1,*) re
    read(1,*) d1
    read(1,*) d2
    read(1,*) phase
    read(1,*) sigma
    read(1,*) lambda
    read(1,*) del
    read(1,*) kappa
    read(1,*) tmax
    read(1,*) icont 
    read(1,*) gmma 
    hx=1.d0/dble(nx)       
    hy=1.d0/dble(ny)
    phase=phase*pi
    Pec=sigma*re    
    t_a=tmax-2.d0*pi
    integrating=.false.
    frame=1
    call initialize_exsolution(ny)
    call initialize_adi(ny,nx,del,re,d1,d2,phase,lambda,Pec,kappa)
    call initialize_vects
    if(navier) then
       call initialize_ns(ny,nx,del,re,d1,d2,phase,gmma)
       open(9,file='trace.dat',status='replace') 
       inquire(iolength=length)  s(0:ny+1,1)
       open(11,file='u.dat',status='unknown',form='unformatted',access='direct',recl=length)
       open(12,file='v.dat',status='unknown',form='unformatted',access='direct',recl=length)
       open(13,file='p.dat',status='unknown',form='unformatted',access='direct',recl=length)
    endif
    open(3,file='final',status='unknown')   
    open(8,file='nusselt.dat',status='unknown',position='append') 
    inquire(iolength=length)  hx
    open(4,file='ffin.dat',status='unknown',form='unformatted',access='direct',recl=length)
    open(10,file='average.dat',status='unknown',form='unformatted',access='direct',recl=length)
    inquire(iolength=length)  s(0:3*ny+3,1)
    open(7,file='solute.dat',status='unknown',form='unformatted',access='direct',recl=length)
    open(14,file='frame_times.dat',status='replace')
    if(icont.eq.1) then
       call read_results
       if(navier) then
          !call self_similar_ns !if starting ns computations from self-similar state
          call read_velocities !if resuming ns computations
          call write_trace
       endif
    else
       t=kick_start(ny,d1,d2,re,phase)
       rewind(8)
       if(navier) then 
          call self_similar_ns 
          call write_trace
       endif
    endif
    dt=2.d0*pi/3000.d0
  end subroutine read_vars

  subroutine initialize_vects
    if(.not.allocated(vectf)) then
       allocate(vectf(0:nx+1))
       allocate(vectfo(0:nx+1))
       allocate(ave(0:nx+1))
       allocate(intgr(0:ny+1))
    end if
    vectf=0.d0
    vectfo=0.d0
    ave=0.d0
    intgr=0.d0
  end subroutine initialize_vects

  subroutine time_integrate(dh,n,a,b,t_int)
    integer, intent(in) :: n
    real(double), intent(in)  :: dh 
    real(double), dimension(0:),  intent(in) :: a,b 
    real(double), dimension(0:), intent(inout) :: t_int 
    t_int(0:n+1)=t_int(0:n+1) +0.5d0*dh*(a(0:n+1)+b(0:n+1))
  end subroutine time_integrate

  real(double) function spatial_average(dh,n,f) result(r)
    integer, intent(in) :: n
    real(double), intent(in)  :: dh 
    real(double), dimension(0:),  intent(in)  :: f !indices 0...n+1
    r=0.5d0*dh*(f(1)+f(n)+2.d0*sum(f(2:n-1)))      !code asumes a staggered
    r=r+ dh*(3.d0*f(n)+f(n+1))/8.d0                !grid is used
    r=r+ dh*(3.d0*f(1)+f(0))/8.d0
  end function spatial_average

  subroutine compute_nusselt                          !since this calculation involves the whole solid-fluid interface, 
    real(double) :: nnt                               !is important to fill the residual boundary points s(2*ny+1,nx+1)
    real(double), dimension(0:nx+1) :: intg           !s(2*ny,nx+1) in a consistent way. I use the no-flux condition also
    intg(0:nx+1)=(s(2*ny+1,0:nx+1)-s(2*ny,0:nx+1))/hy !with this points.
    !sw=oscillations(t,d2)                             
    nnt=spatial_average(hx,nx,intg)/2.d0/pi  !!correct only if d1=d2, and q0=0
    write(8,*) t,nnt
    call flush(8)
  end subroutine compute_nusselt

  subroutine save_results
    integer :: i
    real(double) :: tmp         !for some reason the executable is not storing
    do i=0,ny                   !the ss_solution to file ffin.dat. the solution
       write(4,rec=i+1) vectv(i)!I found is to add the artificial line below
       read(4,rec=i+1) tmp      !which forces the executable to store the result
    enddo
    do i=0,nx+1  
       write(7,rec=i+1) s(0:3*ny+3,i)
    enddo
    rewind(3)
    write(3,*) t
    call flush(3)
  end subroutine save_results

  subroutine read_results
    integer :: i
    real(double), dimension(0:ny) :: exs
    read(3,*) t
    do i=0,ny  
       read(4,rec=i+1) exs(i)
    enddo
    do i=0,nx+1  
       read(7,rec=i+1) s(0:3*ny+3,i)
    enddo
    call recover_ss(ny,t,d1,d2,phase,exs)
  end subroutine read_results

  subroutine save_velocities
    integer :: i
    do i=0,nx+1  
       write(11,rec=i+1) u(0:ny+1,i)
    enddo
    do i=0,nx+1  
       write(12,rec=i+1) v(0:ny+1,i)
    enddo
    do i=0,nx+1  
       write(13,rec=i+1) p(0:ny+1,i)
    enddo
  end subroutine save_velocities
  
  subroutine read_velocities
    integer :: i
    real(double), dimension(0:ny+1,0:nx+1) :: lu,lv,lp
    do i=0,nx+1
       read(11,rec=i+1) lu(0:ny+1,i)
    enddo
    do i=0,nx+1 
       read(12,rec=i+1) lv(0:ny+1,i)
    enddo
    do i=0,nx+1 
       read(13,rec=i+1) lp(0:ny+1,i)
    enddo
    call recover_ns(lu,lv,lp)
  end subroutine read_velocities
  
  subroutine save_averages(n,f)
    integer :: i
    integer, intent(in) :: n
    real(double), intent(inout), dimension (0:) :: f
    f(0:n+1)=f(0:n+1)/2.d0/pi
    do i=0,n+1
       write(10,rec=i+1) f(i)
    enddo
    f=0.d0
  end subroutine save_averages

  integer function make_frame(current,fcount) result(nextframe)
    integer, intent(in) :: current    !current global iteration 
    integer, intent(inout) :: fcount  !current frame number
    integer :: m,niter
    character*30 cmand
    m=nint(2*pi/dt) !dt variable is public in library ADI
    niter=m/25      !25=# of frames desired
    if(fcount.lt.10) write(cmand,555) fcount!cmand contains command for bking up file solute.dat now
    if((fcount.ge.10).and.(fcount.lt.100)) write(cmand,556) fcount
    if((fcount.ge.100).and.(fcount.lt.1000)) write(cmand,557) fcount
    call system(cmand)                      !executes the system shell command containd in cmand
    write(14,*) fcount,t 
    fcount=fcount+1
    nextframe=current+niter                 

555 format('cp solute.dat slt', 1i1,'.dat')
556 format('cp solute.dat slt', 1i2,'.dat')
557 format('cp solute.dat slt', 1i3,'.dat')
  end function make_frame

  real(double) function make_frame2(current,fcount) result(nexttime)
    real(double), intent(in) :: current    !current time
    real(double) :: t_interval
    integer, intent(inout) :: fcount  !current frame number
    character*30 cmand
    t_interval=2*pi/25
    if(fcount.lt.10) write(cmand,555) fcount!cmand contains command for bking up file solute.dat now
    if((fcount.ge.10).and.(fcount.lt.100)) write(cmand,556) fcount
    if((fcount.ge.100).and.(fcount.lt.1000)) write(cmand,557) fcount
    call system(cmand)                      !executes the system shell command containd in cmand
    write(14,*) fcount,t
    fcount=fcount+1
    nexttime=current+t_interval                 

555 format('cp solute.dat slt', 1i1,'.dat')
556 format('cp solute.dat slt', 1i2,'.dat')
557 format('cp solute.dat slt', 1i3,'.dat')
  end function make_frame2

  subroutine write_trace
    write(9,*) t,energy(u,v),ss_energy() 
    call flush(9)
  end subroutine write_trace

  real(double) function energy(bu,bv) result(r)
    !this function calculates the integral in the whole
    !domain. so, bu(0:ny+1,0:nx) and bv(0:ny,0:nx+1) are required
    real(double), dimension(0:,0:), intent(in) :: bu,bv
    real(double), dimension(0:ny+1,0:nx+1) :: lu,lv
    real(double), dimension(0:nx) :: avu
    real(double), dimension(0:ny) :: avv
    real(double)  eu
    
    lu=bu*bu
    lv=bv*bv
    avu=(lu(0,0:nx)+lu(ny+1,0:nx)+7.d0*lu(1,0:nx)+7.d0*lu(ny,0:nx))/8.d0
    do j=2,ny-1
       avu=avu+lu(j,0:nx)
    enddo
    avu=hy*avu
    eu=0.5d0*(avu(0)+avu(nx))
    eu=hx*(eu+sum(avu(1:nx-1)))
    avv=(lv(0:ny,0)+lv(0:ny,nx+1)+7.d0*lv(0:ny,1)+7.d0*lv(0:ny,nx))/8.d0
    do j=2,nx-1
       avv=avv+lv(0:ny,j)
    enddo
    avv=hx*avv
    r=0.5d0*(avv(0)+avv(ny))
    r=hy*(r+sum(avv(1:ny-1)))
    r=r+eu
  end function energy
  
  real(double) function ss_energy() result(er)
    real(double), dimension(0:ny+1,0:nx+1) :: lu,lv
    call self_similar(lu,lv)
    er=energy(lu,lv) 
  end function ss_energy
  
  subroutine self_similar(lu,lv)
    real(double) :: fna,fnad,fnb,fnbd
    real(double), dimension (3) :: sw    
    real(double), dimension (0:nx) :: xi
    real(double), dimension(0:,0:), intent(inout) :: lu,lv
    xi(0:nx)=(/ (dble(j)*hx, j=0,nx) /)
    lu=0.d0
    lv=0.d0
    sw=oscillations(t,d1) 
    fna=sw(1)
    fnad=sw(2)
    sw=oscillations(t+phase,d2)
    fnb=sw(1)
    fnbd=sw(2)
    do j=0,nx+1 
       lv(0:ny,j)=fna*vectv(0:ny)
    enddo
    do j=0,nx
       lu(1:ny,j)=-xi(j)*fnb*(vectv(1:ny)-vectv(0:ny-1))/hy
    enddo
    !-----exterior points---------
    lu(0,1:nx)=2.d0*fnbd*xi(1:nx)-lu(1,1:nx) 
    lu(ny+1,1:nx)=lu(ny,1:nx) 
  end subroutine self_similar

  subroutine finalize
    deallocate(vectf,vectfo,ave,intgr)
    call finalize_adi
    call finalize_exsolution
    if(navier)  then
       call finalize_ns
       close(9)
       close(11) 
       close(12) 
       close(13) 
    endif
    close(3)
    close(8)
    close(4)
    close(10)
    close(7)
    close(14)
  end subroutine finalize
end program solute
           
