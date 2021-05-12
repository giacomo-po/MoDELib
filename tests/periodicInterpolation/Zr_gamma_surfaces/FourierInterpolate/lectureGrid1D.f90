MODULE LectureGrid1DModule

CONTAINS

  SUBROUTINE ReadGridFile1D(inp, nxCol, nCol, nVal, T, f)

    IMPLICIT NONE

    ! Unité de lecture du fichier input
    INTEGER, intent(in) :: inp
    ! Numéro de colonne où lire la valeur de la position
    INTEGER, intent(in) :: nxCol    
    ! Numéro de colonne où lire la valeur de la fonction
    INTEGER, intent(in) :: nCol
    ! Nombre de valeurs de la fonction à interpoler
    INTEGER, intent(in) :: nVal
    ! Période
    REAL(kind(0.d0)), intent(in) :: T
    ! Tableau des valeurs de la fonction à interpoler
    REAL(kind(0.d0)), intent(out), dimension(0:nVal-1) :: f

    CHARACTER(len=50), dimension(1:nxCol-1) :: aTemp
    CHARACTER(len=50), dimension(1:nCol-nxCol-1) :: bTemp

    INTEGER :: i, io
    REAL(kind(0.d0)) :: dx, inv_dx, x, fx
    REAL(kind(0.d0)), parameter :: Zero=1.d-6
    LOGICAL, dimension(0:nVal-1) :: fRead


    dx = T/dble(nVal)
    inv_dx = 1.d0/dx

    fRead(:)=.FALSE.
    f(:)=0.d0

    io = 0
    DO WHILE (io .EQ. 0)
           CALL Comment(inp)
           READ(inp, *, iostat=io) aTemp(1:nxCol-1), x, bTemp(1:nCol-nxCol-1), fx           
           i = nInt(x*inv_dx)
           IF (Abs(dble(i)*dx-x) .GT. Zero*T) THEN
                   WRITE(0,'(a)') "Problème avec la ligne suivante"
                   WRITE(0,*) (Trim(aTemp(i))//" ",i=1,nxCol-1), x, (Trim(bTemp(i))//" ",i=1,nCol-nxCol-1), fx

                   WRITE(0, '(a,g20.6)') "Pas d'indice calculable pour x = ", x
                   STOP
           END IF
           i = Modulo(i, nVal)
           f(i) = fx
           fRead(i) = .TRUE.
    END DO

    DO i=0, nVal-1
       IF (fRead(i).EQV..False.) THEN
               WRITE(0,'(a,i0)') "Pas de valeur trouvée pour l'indice ", i
               STOP
       END IF 
    END DO

  END SUBROUTINE ReadGridFile1D

  SUBROUTINE ReadGridFile1D_linear(inp, nxCol, nCol, nVal, T, f1, f)

    IMPLICIT NONE

    ! Unité de lecture du fichier input
    INTEGER, intent(in) :: inp
    ! Numéro de colonne où lire la valeur de la position
    INTEGER, intent(in) :: nxCol    
    ! Numéro de colonne où lire la valeur de la fonction
    INTEGER, intent(in) :: nCol
    ! Nombre de valeurs de la fonction à interpoler
    INTEGER, intent(in) :: nVal
    ! Période
    REAL(kind(0.d0)), intent(in) :: T
    ! Variation linéaire à retrancher
    REAL(kind(0.d0)), intent(in) :: f1
    ! Tableau des valeurs de la fonction à interpoler
    REAL(kind(0.d0)), intent(out), dimension(0:nVal-1) :: f
    

    CHARACTER(len=50), dimension(1:nxCol-1) :: aTemp
    CHARACTER(len=50), dimension(1:nCol-nxCol-1) :: bTemp

    INTEGER :: i, io
    REAL(kind(0.d0)) :: dx, inv_dx, x, fx
    REAL(kind(0.d0)), parameter :: Zero=1.d-6
    LOGICAL, dimension(0:nVal-1) :: fRead


    dx = T/dble(nVal)
    inv_dx = 1.d0/dx

    fRead(:)=.FALSE.
    f(:)=0.d0

    io = 0
    DO WHILE (io .EQ. 0)
           CALL Comment(inp)
           READ(inp, *, iostat=io) aTemp(1:nxCol-1), x, bTemp(1:nCol-nxCol-1), fx           
           i = nInt(x*inv_dx)
           IF (Abs(dble(i)*dx-x) .GT. Zero*T) THEN
                   WRITE(0,'(a)') "Problème avec la ligne suivante"
                   WRITE(0,*) (Trim(aTemp(i))//" ",i=1,nxCol-1), x, (Trim(bTemp(i))//" ",i=1,nCol-nxCol-1), fx

                   WRITE(0, '(a,g20.6)') "Pas d'indice calculable pour x = ", x
                   STOP
           END IF
           i = Modulo(i, nVal)
           f(i) = fx - x*f1
           fRead(i) = .TRUE.
    END DO

    DO i=0, nVal-1
       IF (fRead(i).EQV..False.) THEN
               WRITE(0,'(a,i0)') "Pas de valeur trouvée pour l'indice ", i
               STOP
       END IF 
    END DO

  END SUBROUTINE ReadGridFile1D_linear


END MODULE LectureGrid1DModule
