!---------------------------------------------------------------------
! T
!	- module containing subroutines for kinetic matrix elements
!---------------------------------------------------------------------
module T
  use key
  implicit none

contains
!---------------------------------------------------------------------
! T_calc
!	- calculates kinetic energy elements for trimmed Hamiltonian
!---------------------------------------------------------------------
! ndim		: int, number of dimensions
! npoints	: 1D int, number of points for each dim
! sys		: int, system coordinates
! delx		: 1D real*8, delta x for each dim
! lb		: 1D real*8, lower bound for each dim
! ub		: 1D real*8, upper bound for each dim
! g		: 2D real*8, diagonal g-matrix elements for each point
! om		: 1D real*8, omega for each dim
! Np		: int, size of trimmed Hamiltonian
! id_vec	: 1D int, id vector of trimmed Hamiltonian
! H		: 2D real*8, trimmed Hamiltonian 

subroutine T_calc(ndim,npoints,sys,delx,lb,ub,g,om,coord,Np,id_vec,H,error)
  implicit none
  real(kind=8), dimension(0:,0:), intent(inout) :: H
  real(kind=8), dimension(0:,0:), intent(in) :: g
  real(kind=8), dimension(0:), intent(in) :: delx,lb,ub,om
  integer, dimension(0:), intent(in) :: npoints,id_vec,coord
  integer, intent(inout) :: error
  integer, intent(in) :: ndim,Np,sys
  real(kind=8), dimension(0:ndim-1) :: dum2
  integer, dimension(0:ndim-1) :: arry_i,arry_j,key
  real(kind=8) :: T
  integer :: fid,dum1
  integer :: i,j
  error = 0
  write(*,*) "Calculating kinetic elements"

  !loop scheme...
  call key_make(ndim,npoints,key)
  do j=0,Np-1
    call key_eval(ndim,key,npoints,coord,id_vec(j),arry_j)
    do i=0,j
      call key_eval(ndim,key,npoints,coord,id_vec(i),arry_i)
      !cartesian coord system
      if (sys .eq. 1) then
        T = T_CART(ndim,npoints,delx,lb,ub,g(j,0:ndim-1),coord,arry_i,arry_j)
      else if (sys .eq. 2) then
        T = T_RAD(ndim,npoints,delx,lb,ub,g(j,0:ndim-1),coord,arry_i,arry_j)
      else if (sys .eq. 3) then
        T = T_POLR(ndim,npoints,delx,lb,ub,g(j,0:ndim-1),coord,arry_i,arry_j)
      else if (sys .eq. 6) then
        T = T_DNC(ndim,npoints,delx,lb,ub,om,coord,arry_i,arry_j)
      else
        write(*,*) "Sorry, that system type is not supported"
        stop
      end if
      H(i,j) = H(i,j) + T
    end do  
  end do

!  write(*,*) 
!  write(*,*) "Hamiltonian after T"
!  do i=0,Np-1
!    write(*,*) H(i,0:Np-1)
!  end do

  write(*,*)
end subroutine T_calc

!---------------------------------------------------------------------
! T_CART
!	- evalutes the kinetic matrix element in cartesian coordinates
!---------------------------------------------------------------------
! ndim		: int, number of dimensions
! npoints	: 1D int, number of points of each dimension
! delx		: 1D real*8, delta x of each dimension
! lb		: 1D real*8, lower bound of each dimension 
! ub		: 1D real*8, upper bound of each dimension
! g		: 1D real*8, gmatrix diagonal of each dimension
! coord		: 1D int, coordinate type of each dimension
! arry_i	: 1D int, i index array
! arry_j	: 1D int, j index array

real(kind=8) function T_CART(ndim,npoints,delx,lb,ub,g,&
                             coord,arry_i,arry_j)
  implicit none
  real(kind=8), dimension(0:), intent(in) :: delx,lb,ub,g
  integer, dimension(0:), intent(in) :: arry_i,arry_j,npoints,coord
  integer, intent(in) :: ndim
  real(kind=8) :: val,pi
  integer :: i,j,k,N
  pi = 3.14159265358979D0
  val = 0.0D0

  do k=0,ndim-1
    if (all(arry_i(0:k-1) .eq. arry_j(0:k-1)) .and. &
        all(arry_i(k+1:ndim-1) .eq. arry_j(k+1:ndim-1)) ) then
      i = arry_i(k)
      j = arry_j(k)

      !infinite cartesian
      if (coord(k) .eq. 1) then
        if (i .ne. j) then
          val = val + ((-1.0D0)**(i - j))*g(k)/((delx(k)*(i - j))**2.0D0)
        else
          val = val + (pi**2.0D0)*g(k)/(6.0D0*delx(k)**2.0D0)
        end if

      !boxed cartesian
      else if (coord(k) .eq. 5) then
        N = npoints(k) + 1
        if (i .ne. j) then
          val = val + 0.25D0*pi**2.0D0*g(k)*((-1.0D0)**(i-j))&
                      /((ub(k)-lb(k))**2.0D0)&
                      *(1.0D0/sin(pi*(i-j)/(2*N))**2.0D0 &
                        - 1.0D0/sin(pi*(i+j)/(2*N))**2.0D0)
        else 
          val = val + 0.25D0*pi**2.0D0*g(k)/((ub(k)-lb(k))**2.0D0) &
                      *((2.0D0*N**2.0D0+1.0D0)/3.0D0 &
                        - 1.0D0/sin(pi*i/N)**2.0D0) 
        end if
      end if
    end if   
  end do 
  T_CART = val

end function T_CART

!---------------------------------------------------------------------
! T_RAD
!	- calculates the kinetic energy matrix in 1D radial coord
!---------------------------------------------------------------------
! ndim		: int, number of dimensions
! npoints	: 1D int, number of points of each dimension
! delx		: 1D real*8, delta x of each dimension
! lb		: 1D real*8, lower bound of each dimension 
! ub		: 1D real*8, upper bound of each dimension
! g		: 1D real*8, gmatrix element for each dimension
! coord		: 1D int, coordinate type of each dimension
! arry_i	: 1D int, i index array
! arry_j	: 1D int, j index array

real(kind=8) function T_RAD(ndim,npoints,delx,lb,ub,g,&
                             coord,arry_i,arry_j)
  implicit none
  real(kind=8), dimension(0:), intent(in) :: delx,lb,ub,g
  integer, dimension(0:), intent(in) :: arry_i,arry_j,npoints,coord
  integer, intent(in) :: ndim
  real(kind=8) :: val,pi
  integer :: i,j,k
  pi = 3.14159265358979D0
  i = arry_i(0)
  j = arry_j(0)
  if (i .ne. j) then
    val = ((-1.0D0)**(i-j))*g(0)/(delx(0)**2.0D0)&
           *(1.0D0/(i-j)**2.0D0 - 1.0D0/(i+j)**2.0D0) 
  else
    val = ((-1.0D0)**(i-j))*g(0)/(delx(0)**2.0D0)&
           *(pi**2.0D0/6.0D0 - 0.25/i**2.0D0)
  end if
  T_RAD = val
end function T_RAD

!---------------------------------------------------------------------
! T_POLR
!       - calculates the kinetic energy matrix in 1D polar coord
!---------------------------------------------------------------------
! ndim          : int, number of dimensions
! npoints       : 1D int, number of points of each dimension
! delx          : 1D real*8, delta x of each dimension
! lb            : 1D real*8, lower bound of each dimension 
! ub            : 1D real*8, upper bound of each dimension
! g		: 1D real*8, gmatrix element for each dimension
! coord         : 1D int, coordinate type of each dimension
! arry_i        : 1D int, i index array
! arry_j        : 1D int, j index array

real(kind=8) function T_POLR(ndim,npoints,delx,lb,ub,g,&
                             coord,arry_i,arry_j)
  implicit none
  real(kind=8), dimension(0:), intent(in) :: delx,lb,ub,g
  integer, dimension(0:), intent(in) :: arry_i,arry_j,npoints,coord
  integer, intent(in) :: ndim
  real(kind=8) :: val,pi
  integer :: i,j,k,N
  pi = 3.14159265358979D0
  i = arry_i(0)
  j = arry_j(0)
  N = npoints(k) + 1
  if (i .ne. j) then
    val = 0.25D0*((-1.0D0)**(i-j))*g(0) &
          *(1.0D0/sin(pi*(i-j)/(2*N))**2.0D0 &
          - 1.0D0/sin(pi*(i+j)/(2*N))**2.0D0)
  else
    val = 0.25D0*g(0) &
          *((2.0D0*N**2.0D0+1.0D0)/3.0D0 &
          - 1.0D0/sin(pi*i/N)**2.0D0)
  end if
  T_POLR = val
end function T_POLR

!---------------------------------------------------------------------
! T_DNC
!	- evalutes the kinetic matrix element in dimensionless
!	  normal coordinates
!---------------------------------------------------------------------
! ndim		: int, number of dimensions
! npoints	: 1D int, number of points of each dimension
! delx		: 1D real*8, delta x of each dimension
! lb		: 1D real*8, lower bound of each dimension 
! ub		: 1D real*8, upper bound of each dimension
! om		: 1D real*8, omega for each dimension
! coord		: 1D int, coordinate type of each dimension
! arry_i	: 1D int, i index array
! arry_j	: 1D int, j index array

real(kind=8) function T_DNC(ndim,npoints,delx,lb,ub,om,&
                             coord,arry_i,arry_j)
  implicit none
  real(kind=8), dimension(0:), intent(in) :: delx,lb,ub,om
  integer, dimension(0:), intent(in) :: arry_i,arry_j,npoints,coord
  integer, intent(in) :: ndim
  real(kind=8) :: val,pi
  integer :: i,j,k,N
  pi = 3.14159265358979D0
  val = 0.0D0

  do k=0,ndim-1
    if (all(arry_i(0:k-1) .eq. arry_j(0:k-1)) .and. &
        all(arry_i(k+1:ndim-1) .eq. arry_j(k+1:ndim-1)) ) then
      i = arry_i(k)
      j = arry_j(k)

      if (i .ne. j) then
        val = val + om(k)*((-1.0D0)**(i - j))/((delx(k)*(i - j))**2.0D0)
      else
        val = val + om(k)*(pi**2.0D0)/(6.0D0*delx(k)**2.0D0)
      end if
    end if   
  end do 
  T_DNC = val

end function T_DNC

!---------------------------------------------------------------------
end module T
