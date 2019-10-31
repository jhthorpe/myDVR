module V
  use key
  implicit none

contains

!---------------------------------------------------------------------
! V_calc
!	- calculates and trims potential energy hamiltonian elements
!---------------------------------------------------------------------
! ndim		: int, number of dimensions
! npoints	: 1D int, number of points per dimension
! delx		: 1D real*8, stepsize for each dimension
! lb		: 1D real*8, lower bound for each dimension
! ub		: 1D real*8, upper bound for each dimension
! N		: int, size of hamiltonian
! pot		: int, potential type
! coord		: 1D int, coordinate type for each dimension
! Np		: int, size of trimmed hamiltonian
! id_vec	: 1D int, id vector of trimmed Hamiltonian
! H		: 2D real*8, trimmed Hamiltonian

subroutine V_calc(ndim,npoints,delx,lb,ub,N,pot,coord,Vc,Np,id_vec,H)
  implicit none
  real(kind=8), dimension(:,:), allocatable, intent(inout) :: H
  integer, dimension(:), allocatable, intent(inout) :: id_vec
  real(kind=8), dimension(0:), intent(in) :: delx,lb,ub
  integer, dimension(0:), intent(in) :: npoints,coord
  real(kind=8), intent(in) :: Vc
  integer, intent(inout) :: Np
  integer, intent(in) :: ndim,pot,N

  real(kind=8), dimension(:), allocatable :: R_temp
  integer, dimension(:), allocatable :: I_temp
  integer, dimension(0:ndim-1) :: arry_i,key
  real(kind=8), dimension(0:ndim-1) :: phi,L
  real(kind=8) :: V
  integer :: i,Vid

  write(*,*) "Calculating potential elements"
  !Harmonic Oscillator
  if (pot .eq. 1) then
    write(*,*) "Enter HO force constants..."
    do i=0,ndim-1
      write(*,*) "dimension ",i+1
      read(*,*) phi(i)
    end do
  !Particle in a Box
  else if (pot .eq. 2) then
!    write(*,*) "Enter box length..."
!    do i=0,ndim-1
!      write(*,*) "dimension ",i+1
!      read(*,*) L(i)
!    end do
  !Particle in a Stadium
  else if (pot .eq. 3) then
    if (ndim .ne. 2) then
      write(*,*) "You need 2 dimensions for PIS"
      stop
    end if
    write(*,*) "Enter X distance"
    read(*,*) L(0)
    write(*,*) "Enter Y distance"
    read(*,*) L(1) 
  !General Potential
  else if (pot .eq. 4) then
    Vid = 200
    write(*,*) "Reading potential from V.in"
    open(file='V.in',unit=Vid,status='old')
  else
    write(*,*) "Sorry, this potential not yet coded"
    STOP
  end if

  !Temporary storage of V_vec and id_vec
  allocate(R_temp(0:N-1))
  allocate(I_temp(0:N-1))
  Np = 0

  !loop scheme...
  call key_make(ndim,npoints,key)
  do i=0,N-1
    call key_eval(ndim,key,npoints,coord,i,arry_i)
    if (pot .eq. 1) then 
      V = V_HO(ndim,npoints,phi,delx,lb,arry_i)
    else if (pot .EQ. 2) then
      V = V_PIB(ndim,npoints,delx,lb,ub,arry_i)
    else if (pot .eq. 3) then
      V = V_PIS(ndim,npoints,L,delx,lb,arry_i)
    else if (pot .eq. 4) then
      V = V_GEN(ndim,npoints,Vid,delx,lb,arry_i,i)
    end if
    if (V .le. Vc) then
      R_temp(Np) = V
      I_temp(Np) = i       
      Np = Np + 1
    end if
  end do

  if (pot .eq. 4) then
    close(unit=Vid)
  end if

  !Store trimmed points and their ids
  write(*,*) "Trimmed Hamiltonian is dimension",Np

  allocate(id_vec(0:Np-1))
  id_vec(0:Np-1) = I_temp(0:Np-1)
  deallocate(I_temp)

  allocate(H(0:Np-1,0:Np-1))
  H = 0.0D0
  do i=0,Np-1
    H(i,i) = R_temp(i)
  end do
  deallocate(R_temp)
  write(*,*)

!  write(*,*) "Writing the hamiltonian"
!  do i=0,Np-1
!    write(*,*) H(i,0:Np-1)
!  end do

end subroutine V_calc

!---------------------------------------------------------------
! V_HO	
!	Calculates multidimensional HO potential
!	at a given point
!---------------------------------------------------------------
! ndim		: int, number of dimensions
! npoints	: 1D int, number of points per dimension
! phi		: 1D real*8, force constants
! delx		: 1D real*8, delta x for dimensions
! lb		: 1D real*8, lower bound for dimensions
! arry_i	: 1D int, gridpoint ids
real(kind=8) function V_HO(ndim,npoints,phi,delx,lb,arry_i)
  implicit none
  real(kind=8), dimension(0:), intent(in) :: delx,phi,lb
  integer, dimension(0:), intent(in) :: arry_i,npoints
  integer, intent(in) :: ndim
  real(kind=8) :: val,x
  integer :: k 
  val = 0.0D0
  do k=0,ndim-1
    !x = delx(k)*(arry_i(k)-npoints(k)/2)
    x = lb(k) + delx(k)*arry_i(k) 
    val = val + 0.5D0*phi(k)*x**2.0D0
  end do
  V_HO = val
end function V_HO

!---------------------------------------------------------------
! V_PIB	
!	Calculates multidimensional PIB potential
!	at a given point
!---------------------------------------------------------------
! ndim		: int, number of dimensions
! npoints	: 1D int, number of points per dimension
! delx		: 1D real*8, delta x for dimensions
! lb		: 1D real*8, lower bound for dimensions
! ub		: 1D real*8, upper bound for dimensions
! arry_i	: 1D int, gridpoint ids

real(kind=8) function V_PIB(ndim,npoints,delx,lb,ub,arry_i)
  implicit none
  real(kind=8), dimension(0:), intent(in) :: delx,lb,ub
  integer, dimension(0:), intent(in) :: arry_i,npoints
  integer, intent(in) :: ndim
  real(kind=8) :: val,x
  integer :: k 
  val = 0.0D0
  do k=0,ndim-1
    !x = delx(k)*(arry_i(k)-npoints(k)/2)
    !if ((x .GT. 0.5D0*L(k)) .OR. (x .LT. -0.5D0*L(k))) val = HUGE(val)
    x = lb(k) + delx(k)*arry_i(k)
    if ((x .lt. lb(k)) .OR. (x .gt. ub(k))) val = HUGE(val)
  end do
  V_PIB = val
end function V_PIB

!---------------------------------------------------------------
! V_PIS
!	Calculates multidimensional PIS potential
!	at a given point
! 	
!	The "student sections", or "field goals" (the
!	semicircles) protrude outsize of the X box length.
!---------------------------------------------------------------
! ndim		: int, number of dimensions
! npoints	: 1D int, number of points per dimension
! L		: 1D real*8, box lengths (0=x,1=y)
! delx		: 1D real*8, delta x for dimensions
! lb		: 1D real*8, lower bound for dimensions
! arry_i	: 1D int, gridpoint ids
real(kind=8) function V_PIS(ndim,npoints,L,delx,lb,arry_i)
  implicit none
  real(kind=8), dimension(0:), intent(in) :: delx,L,lb
  integer, dimension(0:), intent(in) :: arry_i,npoints
  integer, intent(in) :: ndim
  real(kind=8) :: val,x,y,r
  val = 0.0D0
  !y = delx(1)*(arry_i(1)-npoints(1)/2) 
  y = lb(1) + delx(1)*arry_i(1) 
  !if outside y-box
  if (abs(y) .gt. 0.5D0*L(1)) then
    val = HUGE(val)
  ! if inside y-box
  else
    !x = delx(0)*(arry_i(0)-npoints(0)/2)
    x = lb(0) + delx(0)*arry_i(0) 
    !if outside x-box
    if (abs(x) .gt. 0.5D0*L(0)) then
      r = sqrt((abs(x) - 0.5D0*L(0))**2.0D0 &
               + abs(y)**2.0D0)
      if (r .gt. 0.5D0*L(1)) val = HUGE(val)
    end if
  end if 

  V_PIS = val

end function V_PIS

!---------------------------------------------------------------
! V_GEN	
!	Reads in potential from V.in file	
!	Checks the geometries
!---------------------------------------------------------------
! ndim		: int, number of dimensions
! npoints	: 1D int, number of points per dimension
! Vid		: int, input file id 
! delx		: 1D real*8, delta x for dimensions
! lb		: 1D real*8, lower bound for dimensions
! arry_i	: 1D int, gridpoint ids
! num		: int, point number we are supposed to be on

real(kind=8) function V_GEN(ndim,npoints,Vid,delx,lb,arry_i,num)
  implicit none
  real(kind=8), dimension(0:), intent(in) :: delx,lb
  integer, dimension(0:), intent(in) :: arry_i,npoints
  integer, intent(in) :: ndim,Vid,num
  real(kind=8), dimension(0:ndim-1) :: xyz,xyz_real
  real(kind=8) :: val
  integer :: id,k
  read(Vid,*) id,xyz(0:ndim-1),val
  do k=0,ndim-1
!    xyz_real(k) = delx(k)*(arry_i(k)-npoints(k)/2)
    xyz_real(k) = lb(k) + delx(k)*arry_i(k) 
    if (ABS(xyz_real(k) - xyz(k)) .gt. 1.0D-15) then
      write(*,*) "Point",id,"in V.in had a bad xyz value"
      close(unit=Vid)
      stop
    end if
  end do
  if (id .ne. num) then
    write(*,*) "Point ",id," in V.in had a bad point number"
    close(unit=Vid)
    stop
  else 
    V_GEN = val
  end if
end function V_GEN

!---------------------------------------------------------------

end module V
