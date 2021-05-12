MODULE Lecture2DModule

CONTAINS

FUNCTION GetNVal2D(inp, nxCol, nyCol, nCol) RESULT(nVal)
  ! Determine number of values in input file connected to unit inp

  IMPLICIT NONE
  INTEGER :: nVal
  INTEGER, intent(in) :: inp, nxCol, nyCol, nCol

  REAL(kind(0.d0)), dimension(1:nxCol-1) :: aTemp
  REAL(kind(0.d0)), dimension(1:nyCol-nxCol-1) :: bTemp
  REAL(kind(0.d0)), dimension(1:nCol-nyCol-1) :: cTemp
  REAL(kind(0.d0)) :: x, y, f

  INTEGER :: io

  nVal = 0
  io = 0
  DO WHILE (io .EQ. 0)
           CALL Comment(inp)
           READ(inp, *, iostat=io) aTemp(1:nxCol-1), x,  bTemp(1:nyCol-nxCol-1), y, cTemp(1:nCol-nyCol-1), f
           IF (io.EQ.0) nVal=nVal+1
  END DO

END FUNCTION GetNVal2D

SUBROUTINE ReadVal2D(inp, nxCol, nyCol, nCol, nVal, xVal, fVal)

  IMPLICIT NONE
  INTEGER, intent(in) :: inp, nxCol, nyCol, nCol, nVal
  ! Coordinates of all points: xVal(1:2,n)
  REAL(kind(0.d0)), dimension(:,:), intent(out) :: xVal
  ! Value of the function in the corresponding point: fVal(n)
  REAL(kind(0.d0)), dimension(:), intent(out) :: fVal

  REAL(kind(0.d0)), dimension(1:nxCol-1) :: aTemp
  REAL(kind(0.d0)), dimension(1:nyCol-nxCol-1) :: bTemp
  REAL(kind(0.d0)), dimension(1:nCol-nyCol-1) :: cTemp

  INTEGER :: n

  DO n=1, nVal
           CALL Comment(inp)
           READ(inp, *) aTemp(1:nxCol-1), xVal(1,n),  bTemp(1:nyCol-nxCol-1), xVal(2,n), cTemp(1:nCol-nyCol-1), fVal(n)
  END DO

END SUBROUTINE ReadVal2D

END MODULE Lecture2DModule
