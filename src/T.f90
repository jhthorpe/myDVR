module T
  USE key
  implicit none

contains

subroutine T_calc(ndim,npoints,delx,coord,Np,id_vec,H)
  implicit none
  real(kind=8), dimension(0:,0:), intent(inout) :: H
  real(kind=8), dimension(0:), intent(in) :: delx
  integer, dimension(0:), intent(in) :: npoints,id_vec,coord
  integer, intent(in) :: ndim,Np
  integer, dimension(0:ndim-1) :: arry_i,arry_j,key
  integer :: i,j
  write(*,*) "Calculating kinetic elements"
  !loop scheme...
  call key_make(ndim,npoints,key)
  do j=0,Np-1
    call key_eval(ndim,key,npoints,coord,id_vec(j),arry_j)
    do i=0,j
      call key_eval(ndim,key,npoints,coord,id_vec(i),arry_i)
      H(i,j) = H(i,j) + T_eval(ndim,delx,arry_i,arry_j)
    end do  
  end do
end subroutine T_calc

real(kind=8) function T_eval(ndim,delx,arry_i,arry_j)
  implicit none
  real(kind=8), dimension(0:), intent(in) :: delx
  integer, dimension(0:), intent(in) :: arry_i,arry_j
  integer, intent(in) :: ndim
  real(kind=8) :: val,pi
  integer :: i,j,k
  pi = 3.14159265358979D0
  val = 0.0D0
  do k=0,ndim-1
    if (all(arry_i(0:k-1) .eq. arry_j(0:k-1)) .and. &
        all(arry_i(k+1:ndim-1) .eq. arry_j(k+1:ndim-1)) ) then
      i = arry_i(k)
      j = arry_j(k)
      if (i .ne. j) then
        val = val + (1.0D0*(-1.0D0)**(i - j))/((delx(k)*(i - j))**2.0D0)
      else
        val = val + (pi**2.0D0)/(6.0D0*delx(k)**2.0D0)
      end if
    end if   
  end do 
  T_eval = val
end function 


end module T
