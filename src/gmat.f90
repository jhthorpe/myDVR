!---------------------------------------------------------------------
! gmat
!	module containing subrotuines for dealing with the gmatrix
!---------------------------------------------------------------------
module gmat

contains

!---------------------------------------------------------------------
! gmat_get
!	reads in the trimmed, diagonal gmatrix elements from g.in
!---------------------------------------------------------------------
! ndim		: int, number of dimensions
! Np		: int, number of trimmed gridpoints
! id_vec	: 1D int, ids of untrimmed points
! g		: 2D real*8, diagonal gmatrix elements for each point

subroutine gmat_get(ndim,Np,id_vec,g)
  implicit none
  real(kind=8), dimension(0:,0:), intent(inout) :: g
  integer, dimension(0:), intent(in) :: id_vec
  integer, intent(in) :: ndim,Np
  real(kind=8), dimension(0:ndim-1) :: dum2
  integer :: i,line,dum1

  write(*,*) "Reading trimmed g-matrix elements from g.in"
  open(file='g.in',unit=100,status='old')
  line = 0
  do i=0,Np-1
    !fast forward to next allowed point in g.in
    do while (line .lt. id_vec(i))
      read(100,*) 
      line = line + 1
    end do
    read(100,*) dum1,dum2,g(i,0:ndim-1)
    line = line + 1
  end do
  close(unit=100)
  write(*,*)

end subroutine gmat_get
!---------------------------------------------------------------------

end module gmat
