module V
  use key
  implicit none

contains

subroutine V_calc(ndim,npoints,delx,N,pot,Vc,Np,id_vec,H)
  implicit none
  real(kind=8), dimension(:,:), allocatable, intent(inout) :: H
  integer, dimension(:), allocatable, intent(inout) :: id_vec
  real(kind=8), dimension(0:), intent(in) :: delx
  integer, dimension(0:), intent(in) :: npoints
  real(kind=8), intent(in) :: Vc
  integer, intent(inout) :: Np
  integer, intent(in) :: ndim,pot,N

  real(kind=8), dimension(:), allocatable :: R_temp
  integer, dimension(:), allocatable :: I_temp
  integer, dimension(0:ndim-1) :: arry_i,key
  real(kind=8), dimension(0:ndim-1) :: phi,L
  real(kind=8) :: V
  integer :: i

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
    write(*,*) "Enter box length..."
    do i=0,ndim-1
      write(*,*) "dimension ",i+1
      read(*,*) L(i)
    end do
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
    call key_eval(ndim,key,npoints,i,arry_i)
    if (pot .eq. 1) then 
      V = V_HO(ndim,npoints,phi,delx,arry_i)
    else if (pot .EQ. 2) then
      V = V_PIB(ndim,npoints,L,delx,arry_i)
    else if (pot .eq. 3) then
      V = V_PIS(ndim,npoints,L,delx,arry_i)
    end if
    if (V .le. Vc) then
      R_temp(Np) = V
      I_temp(Np) = i       
      Np = Np + 1
    end if
  end do

  !Store trimmed points and their ids

  allocate(id_vec(0:Np-1))
  id_vec(0:Np-1) = I_temp(0:Np-1)
  deallocate(I_temp)

  allocate(H(0:Np-1,0:Np-1))
  H = 0.0D0
  do i=0,Np-1
    H(i,i) = R_temp(i)
  end do
  deallocate(R_temp)

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
! arry_i	: 1D int, gridpoint ids
real(kind=8) function V_HO(ndim,npoints,phi,delx,arry_i)
  implicit none
  real(kind=8), dimension(0:), intent(in) :: delx,phi
  integer, dimension(0:), intent(in) :: arry_i,npoints
  integer, intent(in) :: ndim
  real(kind=8) :: val,x
  integer :: k 
  val = 0.0D0
  do k=0,ndim-1
    x = delx(k)*(arry_i(k)-npoints(k)/2)
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
! L		: 1D real*8, box lengths
! delx		: 1D real*8, delta x for dimensions
! arry_i	: 1D int, gridpoint ids
real(kind=8) function V_PIB(ndim,npoints,L,delx,arry_i)
  implicit none
  real(kind=8), dimension(0:), intent(in) :: delx,L
  integer, dimension(0:), intent(in) :: arry_i,npoints
  integer, intent(in) :: ndim
  real(kind=8) :: val,x
  integer :: k 
  val = 0.0D0
  do k=0,ndim-1
    x = delx(k)*(arry_i(k)-npoints(k)/2)
    if ((x .GT. 0.5D0*L(k)) .OR. (x .LT. -0.5D0*L(k))) val = HUGE(val)
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
! arry_i	: 1D int, gridpoint ids
real(kind=8) function V_PIS(ndim,npoints,L,delx,arry_i)
  implicit none
  real(kind=8), dimension(0:), intent(in) :: delx,L
  integer, dimension(0:), intent(in) :: arry_i,npoints
  integer, intent(in) :: ndim
  real(kind=8) :: val,x,y,r
  val = 0.0D0
  y = delx(1)*(arry_i(1)-npoints(1)/2) 
  !if outside y-box
  if (abs(y) .gt. 0.5D0*L(1)) then
    val = HUGE(val)
  ! if inside y-box
  else
    x = delx(0)*(arry_i(0)-npoints(0)/2)
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

end module V
