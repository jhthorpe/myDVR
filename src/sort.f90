module sort
  implicit none

contains

subroutine sort_eval(N,R_eval,I_eval,Psi)
  implicit none
  real(kind=8), dimension(0:,0:), intent(inout) :: Psi
  real(kind=8), dimension(0:), intent(inout) :: R_eval,I_eval
  integer, intent(in) :: N
  
  real(kind=8), dimension(0:N-1) :: vec1
  real(kind=8) :: val1
  integer :: i,j
  
  !Insertion Sort because I am lazy
  !Probably should be merge sort
  i = 1
  do while (i .lt. N) 
    j = i
    do while (j .ge. 1 .and. R_eval(j-1) .gt. R_eval(j)) 
      val1 = R_eval(j) 
      R_eval(j) = R_eval(j-1)
      R_eval(j-1) = val1

      val1 = I_eval(j) 
      I_eval(j) = I_eval(j-1)
      I_eval(j-1) = val1

      vec1 = Psi(0:N-1,j)
      Psi(0:N-1,j) = Psi(0:N-1,j-1)
      Psi(0:N-1,j-1) = vec1
      
      if (j .eq. 1) exit
      j = j - 1
    end do
    i = i + 1
  end do

end subroutine sort_eval

end module sort
