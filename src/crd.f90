!---------------------------------------------------------------------
! crd	
!	- module containing subroutines for dealing with coordinate
!	  systems
!---------------------------------------------------------------------
module crd
  implicit none

contains

!---------------------------------------------------------------------
! crd_get
!	- determines the coordinate system
!---------------------------------------------------------------------
! ndim		: int, number of dimensions
! coord		: 1D int, coordinate type of each dim
! sys		: int, coordinate type of the system
! error		: int, error code

subroutine crd_get(ndim,coord,sys,error)
  implicit none
  integer, dimension(0:), intent(in) :: coord
  integer, intent(inout) :: sys,error
  integer, intent(in) :: ndim
  logical :: cart,rctl,dnc
  integer :: i
  error = 0

  write(*,*) "Determining coordinate system"
  !figure out coordinate system, sys
  ! 1 -> cartesisan
  ! 2 -> radial 1D
  ! 3 -> polar 1D
  ! 4 -> azimuthal 1D
  ! 5 -> rectilinear
  ! 6 -> dimensionless normal coordinates (DNC)

  !one dimensional coordiantes
  if (ndim .eq. 1) then
    if (coord(0) .eq. 1 .or. coord(0) .eq. 5) then
      sys = 1
    else if (coord(0) .eq. 2) then
      sys = 2
    else if (coord(0) .eq. 3) then
      sys = 3
    else if (coord(0) .eq. 4) then
      sys = 4
    else if (coord(0) .eq. 6) then
      sys = 6
    else
      sys = 0
    end if

  !multidimensional coordinates
  else if (ndim .gt. 1) then
    cart = .true.
    rctl = .true.
    dnc = .true.
    do i=0,ndim-1
      if (.not. cart .and. .not. rctl) then
        exit
      else if (coord(i) .ne. 1 .and. coord(i) .ne. 5) then
        cart = .false.
      else if (coord(i) .ne. 2 .and. coord(i) .ne. 3 &
               .and. coord(i) .ne. 4) then
        rctl = .false.
      else if (coord(i) .ne. 6) then
        dnc = .false.
      end if
    end do
    if (cart) then
      sys = 1
    else if (rctl) then
      sys = 5
    else if (dnc) then
      sys = 6
    else
      sys = 0
    end if
  end if

  !print
  if (sys .eq. 0) then
    write(*,*) "Coordinate system is unsupported"
    error = 1
    return
  else if (sys .eq. 1) then
    write(*,*) "Coordinate system is Cartesian"
  else if (sys .eq. 2) then
    write(*,*) "Coordinate system is 1D radial"
  else if (sys .eq. 3) then
    write(*,*) "Coordinate system is 1D polar"
  else if (sys .eq. 4) then
    write(*,*) "Coorindate system is 1D azimuthal"
  else if (sys .eq. 5) then
    write(*,*) "Coordinate system is rectilinear"
  else if (sys .eq. 6) then
    write(*,*) "Coordinate system is dimensionless normal coordinates"
  end if
  write(*,*)

end subroutine crd_get
!---------------------------------------------------------------------

end module crd
