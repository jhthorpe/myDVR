module key
  implicit none

contains

!Generate key
subroutine key_make(ndim,points,key)
  implicit none
  integer, dimension(0:), intent(inout) :: key
  integer, dimension(0:), intent(in) :: points
  integer, intent(in) :: ndim 
  integer :: i
  key(0) = 1
  do i=1,ndim-1
    key(i) = key(i-1)*points(i-1) 
  end do
end subroutine key_make

!convert key to index array 
subroutine key_eval(ndim,key,points,idx,arry)
  implicit none
  integer, dimension(0:), intent(inout) :: arry
  integer, dimension(0:), intent(in) :: key,points
  integer, intent(in) :: ndim,idx
  integer :: i
  do i=0,ndim-1
    arry(i) = MOD(idx/key(i),points(i)) 
  end do
end subroutine key_eval

end module key
