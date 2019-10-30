module key
  implicit none

contains

!---------------------------------------------------------------------
! key_make
!	- generates the key to roll index arrays into one index
!---------------------------------------------------------------------
! ndim		: int, number of dimensions
! key		: 1D int, key
! npoints	: 1D int, number of points per dimension

subroutine key_make(ndim,npoints,key)
  implicit none
  integer, dimension(0:), intent(inout) :: key
  integer, dimension(0:), intent(in) :: npoints
  integer, intent(in) :: ndim 
  integer :: i
  key(0) = 1
  do i=1,ndim-1
    key(i) = key(i-1)*npoints(i-1) 
  end do
end subroutine key_make

!---------------------------------------------------------------------
! key_eval
!	- converts a rolled index into index array via key
!---------------------------------------------------------------------
! ndim		: int, number of dimensions
! key		: 1D int, key
! npoints	: 1D int, number of points per dimension
! coord		: 1D int, coordinate types per dimension
! idx		: int, rolled index
! arry		: 1D int, unrolled index array 

subroutine key_eval(ndim,key,npoints,coord,idx,arry)
  implicit none
  integer, dimension(0:), intent(inout) :: arry
  integer, dimension(0:), intent(in) :: key,npoints,coord
  integer, intent(in) :: ndim,idx
  integer :: i
  do i=0,ndim-1
    arry(i) = MOD(idx/key(i),npoints(i)) 
  end do
  arry = arry + 1
end subroutine key_eval

!---------------------------------------------------------------------
end module key
