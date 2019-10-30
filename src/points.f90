module points
  use key 
  implicit none
  
contains

subroutine points_full(ndim,N,npoints,delx,lb,ub,coord)
  implicit none
  real(kind=8), dimension(0:), intent(in) :: delx,lb,ub
  integer, dimension(0:), intent(in) :: npoints,coord
  integer, intent(in) :: ndim,N
  real(kind=8), dimension(0:ndim-1) :: xyz
  integer, dimension(0:ndim-1) :: arry,key
  integer :: i,j
  write(*,*)
  write(*,*) "Writing gridpoints to full_grid.txt"
  open(file='full_grid.txt',unit=100,status='replace')
  call key_make(ndim,npoints,key)
  do i=0,N-1
    call key_eval(ndim,key,npoints,coord,i,arry)
    do j=0,ndim-1
      !xyz(j) = delx(j)*(arry(j)-npoints(j)/2) 
      xyz(j) = lb(j) + delx(j)*arry(j) 
    end do
    write(100,*) i,xyz(0:ndim-1)
  end do
  close(unit=100)
end subroutine points_full

subroutine points_trim(ndim,Np,npoints,delx,coord,id_vec,H)
  implicit none
  real(kind=8), dimension(0:,0:), intent(in) :: H
  real(kind=8), dimension(0:), intent(in) :: delx
  integer, dimension(0:), intent(in) :: npoints,id_vec,coord
  integer, intent(in) :: ndim,Np
  real(kind=8), dimension(0:ndim-1) :: xyz
  integer, dimension(0:ndim-1) :: arry,key
  integer :: i,j
  write(*,*)
  write(*,*) "Writing gridpoints to trim_grid.txt"
  write(*,*) "Writing potentials to trim_V.txt"
  open(file='trim_grid.txt',unit=100,status='replace')
  open(file='trim_V.txt',unit=101,status='replace')
  call key_make(ndim,npoints,key)
  do i=0,Np-1
    call key_eval(ndim,key,npoints,coord,id_vec(i),arry)
    do j=0,ndim-1
      xyz(j) = delx(j)*(arry(j)-npoints(j)/2)
    end do
    write(100,*) i,xyz(0:ndim-1)
    write(101,*) i,H(i,i)
  end do
  close(unit=100)
  close(unit=101)
end subroutine points_trim

end module points
