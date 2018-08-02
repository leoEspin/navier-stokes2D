!!!------------------------------------------------------------------!!!
!!module for solving the elliptic equation c_x p_xx +c_y p_yy = rhs   !!
!!with hockney's method. Requires the FFTW library, version 2.1.5     !!
!!Code assumes that system to be solved has a                         !!
!!Toeplitz-Symmetric-Tridiagonal (TST) structure                      !!
!!                                                                    !!  
!!Contains thomas subroutine for solving diagonally dominant, tri-    !!
!!diagonal systems                                                    !!
!!                                                                    !!  
!!Author: Leonardo Espin.   Modified: 1/14/2011                       !! 
!!!------------------------------------------------------------------!!!
module poisson_solver
  implicit none
  private
  public :: hockney,initialize_hockney,finalize_hockney
  integer, parameter :: double=kind(1.0d0)
  integer fftw_forward,fftw_backward
  parameter (fftw_forward=-1,fftw_backward=1)
  integer fftw_estimate,fftw_measure
  parameter (fftw_estimate=0,fftw_measure=1)
  integer fftw_threadsafe
  parameter (fftw_threadsafe=128)
  real(double), parameter :: pi=3.1415926535897931!4.d0*atan(1.d0)!
  integer(kind=8) plan
contains
  subroutine initialize_hockney(m)
    integer, intent(in) :: m
    !integer :: status
    call fftw_f77_create_plan(plan,2*m+2,fftw_backward,fftw_estimate+fftw_threadsafe)

    !!!!!if P has to many rows (large columns) it would be convenient to use 
    !!!!!a parallel fourier transform. initialize parallel fftw with line below
    !call fftw_f77_threads_init(status) 
    !returned value in status should be 0 if everything is ok, uncomment thread fftw
    !transforms in hockney subroutine and comment serial transforms
  end subroutine initialize_hockney

  subroutine hockney(m,n,dy,dx,cy,cx,lp,lrhs)
    integer, intent(in) :: m,n
    integer :: i,j
    real(double), dimension(0:,0:), intent(out) :: lp 
    real(double), dimension(0:,0:), intent(in) :: lrhs
    real(double), intent(in) ::  dx,dy,cy,cx
    real(double) :: kt
    complex(double) :: coeff(2*m+2),tmp(2*m+2)  
    real(double) :: a(n),b(n),c(n),x(n),lambs(m),nrhs(m,n)

    kt=sqrt(2.d0/dble(m+1))
    lambs=0.d0
    nrhs=0.d0
    coeff=(0.d0,0.d0)
    tmp=(0.d0,0.d0)
    lambs(1:m)=(/ (i, i=1,m) /)
    lambs(1:m)=-2.d0*cx/dx/dx -2.d0*cy*(1.d0 -cos(pi*lambs(1:m)/dble(m+1)))/dy/dy  
    a(1:n)=cx/dx/dx 
    c(1:n)=a(1:n)
!$omp parallel num_threads(8) default(shared) firstprivate(coeff,tmp) private(b,x)
!$omp do
    do i=1,n
       coeff(2:m+1)=cmplx(lrhs(1:m,i),0.d0,double)
       call fftw_f77_one(plan,coeff,tmp)
       !call fftw_f77_threads_one(4,plan,coeff,tmp)
       nrhs(1:m,i)=kt*aimag(tmp(2:m+1))
    enddo
!$omp end do
!$omp do 
    do j=1,m
       b(1:n)=lambs(j)
       x(1:n)=nrhs(j,1:n)
       call thomas(n,x,a,b,c) 
       lp(j,1:n)=x(1:n)
    enddo
!$omp end do
!$omp do 
    do i=1,n
       coeff(2:m+1)=cmplx(lp(1:m,i),0.d0,double)
       call fftw_f77_one(plan,coeff,tmp)
       !call fftw_f77_threads_one(4,plan,coeff,tmp)
       lp(1:m,i)=kt*aimag(tmp(2:m+1))
    enddo
!$omp end do
!$omp end parallel
  end subroutine hockney

  subroutine thomas(n,x,a,b,c)
    integer, intent(in) :: n
    integer k,i   
    real(double), dimension(:), intent(inout) :: x  !at input x is acting as rhs of the system
    real(double), dimension(:), intent(in) :: a,b,c !as output is the solution
    real(double) :: p(n+1),q(n+1)                   
    real(double) :: denom
    p=0.0d0
    q=0.0d0
    do i=1,n
       denom=a(i)*p(i)+b(i);
       p(i+1)=-c(i)/denom;
       q(i+1)=(x(i)-a(i)*q(i))/denom;
    enddo
    x(n)=q(n+1) !bdry conditions of tridg system are x(0)=0,x(n+1)=0
    do k=n-1,1,-1
       x(k)=p(k+1)*x(k+1)+q(k+1);
    enddo
  end subroutine thomas
  
  subroutine finalize_hockney
    call fftw_f77_destroy_plan(plan)
  end subroutine finalize_hockney
end module poisson_solver
