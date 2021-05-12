MODULE Lecture1DModule

CONTAINS

FUNCTION GetNVal1D(inp, nxCol, nCol) RESULT(nVal)
  ! Determine number of values in input file connected to unit inp

  IMPLICIT NONE
  INTEGER :: nVal
  INTEGER, intent(in) :: inp, nxCol, nCol

  REAL(kind(0.d0)), dimension(1:nxCol-1) :: aTemp
  REAL(kind(0.d0)), dimension(1:nCol-nxCol-1) :: cTemp
  REAL(kind(0.d0)) :: x, f

  INTEGER :: io

  nVal = 0
  io = 0
  DO WHILE (io .EQ. 0)
           CALL Comment(inp)
           READ(inp, *, iostat=io) aTemp(1:nxCol-1), x,  cTemp(1:nCol-nxCol-1), f
           IF (io.EQ.0) nVal=nVal+1
  END DO

END FUNCTION GetNVal1D

SUBROUTINE ReadVal1D(inp, nxCol, nCol, nVal, xVal, fVal)

  IMPLICIT NONE
  INTEGER, intent(in) :: inp, nxCol, nCol, nVal
  ! Coordinates of all points: xVal(n)
  REAL(kind(0.d0)), dimension(:), intent(out) :: xVal
  ! Value of the function in the corresponding point: fVal(n)
  REAL(kind(0.d0)), dimension(:), intent(out) :: fVal

  REAL(kind(0.d0)), dimension(1:nxCol-1) :: aTemp
  REAL(kind(0.d0)), dimension(1:nCol-nxCol-1) :: cTemp

  INTEGER :: n

  DO n=1, nVal
           CALL Comment(inp)
           READ(inp, *) aTemp(1:nxCol-1), xVal(n),  cTemp(1:nCol-nxCol-1), fVal(n)
  END DO

END SUBROUTINE ReadVal1D

END MODULE Lecture1DModule
