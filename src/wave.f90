module wave
  use key

contains

subroutine wave_print(ndim,Np,npoints,delx,id_vec,neig,eval,Psi)
  implicit none
  
  real(kind=8), dimension(0:,0:), intent(in) :: Psi
  real(kind=8), dimension(0:), intent(in) :: delx,eval
  integer, dimension(0:), intent(in) :: npoints,id_vec
  integer, intent(in) :: ndim,Np,neig
  real(kind=8), dimension(0:ndim-1) :: xyz
  integer, dimension(0:ndim-1) :: arry,key
  integer :: i,j,k,foff

  WRITE(*,*) "Wavefunction plots saved as 'wave_X.dat'"
  call execute_command_line('rm wave_[0-9]*.dat')

  foff = 100
  call wave_open(neig,foff)
  call key_make(ndim,npoints,key)

  do j=0,Np-1
    call key_eval(ndim,key,npoints,id_vec(j),arry)
    do i=0,neig-1
      do k=0,ndim-1
        xyz(k) = delx(k)*(arry(k)-npoints(k)/2)
      end do
      write(foff+i,*) j,xyz(0:ndim-1),Psi(j,i)
    end do
  end do

  call wave_close(neig,foff)

end subroutine wave_print

!---------------------------------------------------------------------
! wave_open
!       - opens wave_X.dat files 
!---------------------------------------------------------------------
! Np 		: int, number of points
! foff		: int, unit offset
subroutine wave_open(Np,foff)
  implicit none
  integer, intent(in) :: Np,foff
  character(len=1024) :: fname
  integer :: i,fid
  do i=0,Np-1
    fid = foff + i
    call wave_name(i,fname)
    open(file=TRIM(fname),unit=fid,status='replace')
  end do
end subroutine wave_open

!---------------------------------------------------------------------
! wave_close
!       - closes wave_X.dat files 
!---------------------------------------------------------------------
! Np		: int, number of points
! foff          : int, unit offset
subroutine wave_close(Np,foff)
  implicit none
  integer, intent(in) :: Np,foff
  integer :: i,fid
  do i=0,Np-1
    fid = foff + i
    close(unit=fid)
  end do
end subroutine wave_close

!---------------------------------------------------------------------
! wave_name 
!       - formats wave_X.dat file name
!---------------------------------------------------------------------
! id            : int, id
! fname         : int, filename
subroutine wave_name(id,fname)
  implicit none
  character(len=1024), intent(inout) :: fname
  integer, intent(in) :: id
  character(len=1024) :: str_fmt
  IF (id .LT. 10) THEN
    str_fmt = "(A5,I1,A4)"
  ELSE IF (id .GE. 10 .AND. id .LT. 100) THEN
    str_fmt = "(A5,I2,A4)"
  ELSE IF (id .GE. 100 .AND. id .LT. 1000) THEN
    str_fmt = "(A5,I3,A4)"
  ELSE IF (id .GE. 1000 .AND. id .LT. 10000) THEN
    str_fmt = "(A5,I4,A4)"
  ELSE
    WRITE(*,*) "More cases needed in wave_name"
    RETURN
  END IF
  WRITE(fname,str_fmt) "wave_",id,".dat"
end subroutine wave_name

end module wave
