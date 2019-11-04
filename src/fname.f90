module fname

contains

!---------------------------------------------------------------------
! fname_wave
!       - formats wave_X.dat file name
!---------------------------------------------------------------------
! id            : int, id
! fname         : int, filename
subroutine fname_wave(id,fname)
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
    WRITE(*,*) "More cases needed in fname_wave"
    RETURN
  END IF
  WRITE(fname,str_fmt) "wave_",id,".dat"
end subroutine fname_wave

!---------------------------------------------------------------------
! fname_espc
!       - formats espc_X.dat file name
!---------------------------------------------------------------------
! id            : int, id
! fname         : int, filename
subroutine fname_espc(id,fname)
  implicit none
  character(len=1024), intent(inout) :: fname
  integer, intent(in) :: id
  character(len=1024) :: str_fmt
  IF (id+1 .LT. 10) THEN
    str_fmt = "(A5,I1,A4)"
  ELSE IF (id+1 .GE. 10 .AND. id+1 .LT. 100) THEN
    str_fmt = "(A5,I2,A4)"
  ELSE IF (id+1 .GE. 100 .AND. id+1 .LT. 1000) THEN
    str_fmt = "(A5,I3,A4)"
  ELSE IF (id+1 .GE. 1000 .AND. id+1 .LT. 10000) THEN
    str_fmt = "(A5,I4,A4)"
  ELSE
    WRITE(*,*) "More cases needed in fname_espc"
    RETURN
  END IF
  WRITE(fname,str_fmt) "espc_",id+1,".dat"
end subroutine fname_espc
!---------------------------------------------------------------------

end module fname
