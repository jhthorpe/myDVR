module points
  use key 
  implicit none
  
contains
!---------------------------------------------------------------------
! points_full
!	- prints full list of possible gridpoints
!---------------------------------------------------------------------
! ndim		: int, number of dimensions
! N		: int, size of Hamiltonian
! npoints	: 1D int, number of points for each dim
! delx		: 1D real*8, delta x of each dim
! lb		: 1D real*8, lower bound for each dim
! ub		: 1D real*8, upper bound for each dim
! coord		: 1D int, coordinate type for each dim
subroutine points_full(ndim,N,npoints,delx,lb,ub,coord)
  implicit none
  real(kind=8), dimension(0:), intent(in) :: delx,lb,ub
  integer, dimension(0:), intent(in) :: npoints,coord
  integer, intent(in) :: ndim,N
  real(kind=8), dimension(0:ndim-1) :: xyz
  integer, dimension(0:ndim-1) :: arry,key
  integer :: i,j
  write(*,*) "Writing gridpoints to full_grid.txt"
  open(file='full_grid.txt',unit=100,status='replace')
  call key_make(ndim,npoints,key)
  do i=0,N-1
    call key_eval(ndim,key,npoints,coord,i,arry)
    do j=0,ndim-1
      xyz(j) = lb(j) + delx(j)*arry(j) 
    end do
    write(100,*) i,xyz(0:ndim-1)
  end do
  close(unit=100)
  write(*,*)
end subroutine points_full

!---------------------------------------------------------------------
! points_trim
!	- print list of trimmed gridpoints
!---------------------------------------------------------------------
! ndim		: int, number of dimensions
! Np		: int, size of trimmed Hamiltonian
! npoints	: 1D int, number of points for each dim
! delx		: 1D real*8, delta x of each dim
! lb		: 1D real*8, lower bound for each dim
! coord		: 1D int, coordinate type for each dim

subroutine points_trim(ndim,Np,npoints,delx,lb,coord,id_vec,H)
  implicit none
  real(kind=8), dimension(0:,0:), intent(in) :: H
  real(kind=8), dimension(0:), intent(in) :: delx,lb
  integer, dimension(0:), intent(in) :: npoints,id_vec,coord
  integer, intent(in) :: ndim,Np
  real(kind=8), dimension(0:ndim-1) :: xyz
  integer, dimension(0:ndim-1) :: arry,key
  integer :: i,j
  write(*,*) "Writing gridpoints to trim_grid.txt"
  write(*,*) "Writing potentials to trim_V.txt"
  open(file='trim_grid.txt',unit=100,status='replace')
  open(file='trim_V.txt',unit=101,status='replace')
  call key_make(ndim,npoints,key)
  do i=0,Np-1
    call key_eval(ndim,key,npoints,coord,id_vec(i),arry)
    do j=0,ndim-1
      !xyz(j) = delx(j)*(arry(j)-npoints(j)/2)
      xyz(j) = lb(j) + delx(j)*arry(j)
    end do
    write(100,*) i,xyz(0:ndim-1)
    write(101,*) i,H(i,i)
  end do
  close(unit=100)
  close(unit=101)
  write(*,*)
end subroutine points_trim
!---------------------------------------------------------------------

end module points
