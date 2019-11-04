module fprint
  use key
  use fname

contains

!---------------------------------------------------------------------
! fprint_espc
!	- generates eigenspctrum plotting files
!---------------------------------------------------------------------
! ndim		: int, number of dimensions
! lb		: 1D real*8, lower bounds
! ub		: 1D real*8, upper bounds
! neig		: int, number of eigenvalues
! eval		: 1D real*8, list of eigenvalues

subroutine fprint_espc(ndim,lb,ub,neig,eval)
  implicit none
  real(kind=8), dimension(0:), intent(in) :: lb,ub,eval
  integer, intent(in) :: ndim,neig
  integer :: i,j,foff
  write(*,*) "Eigenspectrum plots saved as 'espc_X.dat"
  write(*,*) "X is the x'th coordinate"
  call execute_command_line('rm espc_[0-9]*.dat')
  
  foff = 200
  call fprint_espc_open(ndim,foff)

  do i=0,neig-1
    do j=0,ndim-1
      write(foff+j,*) i,lb(j),eval(i)
      write(foff+j,*) i,ub(j),eval(i)
      write(foff+j,*)
    end do
  end do

  call fprint_espc_close(ndim,foff)
  write(*,*) 

end subroutine fprint_espc
!---------------------------------------------------------------------
! fprint_wave
!	- generates wavefunction plotting files
!---------------------------------------------------------------------
subroutine fprint_wave(ndim,Np,npoints,delx,lb,coord,id_vec,&
                      neig,eval,Psi)
  implicit none
  
  real(kind=8), dimension(0:,0:), intent(in) :: Psi
  real(kind=8), dimension(0:), intent(in) :: delx,eval,lb
  integer, dimension(0:), intent(in) :: npoints,id_vec,coord
  integer, intent(in) :: ndim,Np,neig
  real(kind=8), dimension(0:ndim-1) :: xyz
  integer, dimension(0:ndim-1) :: arry,key
  real(kind=8), dimension(0:neig-1) :: yold 
  integer :: i,j,k,foff

  write(*,*) "Wavefunction plots saved as 'wave_X.dat'"
  write(*,*) "X is the x'th wavefunction"
  call execute_command_line('rm wave_[0-9]*.dat')

  foff = 100
  call fprint_wave_open(neig,foff)
  call key_make(ndim,npoints,key)

  do j=0,Np-1
    call key_eval(ndim,key,npoints,coord,id_vec(j),arry)
    do i=0,neig-1
      do k=0,ndim-1
        xyz(k) = lb(k) + delx(k)*arry(k) 
      end do
      if (ndim .eq. 2) then
        if (yold(i) .ne. xyz(1) .and. j .gt. 1) write(foff+i,*) 
        yold(i) = xyz(1)
      end if
      write(foff+i,*) j,xyz(0:ndim-1),Psi(j,i)
    end do
  end do

  call fprint_wave_close(neig,foff)
  write(*,*)

end subroutine fprint_wave

!---------------------------------------------------------------------
! fprint_wave_open
!       - opens wave_X.dat files 
!---------------------------------------------------------------------
! Np 		: int, number of points
! foff		: int, unit offset
subroutine fprint_wave_open(Np,foff)
  implicit none
  integer, intent(in) :: Np,foff
  character(len=1024) :: fname
  integer :: i,fid
  do i=0,Np-1
    fid = foff + i
    call fname_wave(i,fname)
    open(file=TRIM(fname),unit=fid,status='replace')
  end do
end subroutine fprint_wave_open

!---------------------------------------------------------------------
! fprint_wave_close
!       - closes wave_X.dat files 
!---------------------------------------------------------------------
! Np		: int, number of points
! foff          : int, unit offset
subroutine fprint_wave_close(Np,foff)
  implicit none
  integer, intent(in) :: Np,foff
  integer :: i,fid
  do i=0,Np-1
    fid = foff + i
    close(unit=fid)
  end do
end subroutine fprint_wave_close

!---------------------------------------------------------------------
! fprint_espc_open
!       - opens espc_X.dat files 
!---------------------------------------------------------------------
! ndim 		: int, number of dimensions
! foff		: int, unit offset
subroutine fprint_espc_open(ndim,foff)
  implicit none
  integer, intent(in) :: ndim,foff
  character(len=1024) :: fname
  integer :: i,fid
  do i=0,ndim-1
    fid = foff + i
    call fname_espc(i,fname)
    open(file=TRIM(fname),unit=fid,status='replace')
  end do
end subroutine fprint_espc_open

!---------------------------------------------------------------------
! fprint_espc_close
!       - closes wave_X.dat files 
!---------------------------------------------------------------------
! ndim		: int, number of dimensions
! foff          : int, unit offset
subroutine fprint_espc_close(ndim,foff)
  implicit none
  integer, intent(in) :: ndim,foff
  integer :: i,fid
  do i=0,ndim-1
    fid = foff + i
    close(unit=fid)
  end do
end subroutine fprint_espc_close
!---------------------------------------------------------------------

end module fprint
