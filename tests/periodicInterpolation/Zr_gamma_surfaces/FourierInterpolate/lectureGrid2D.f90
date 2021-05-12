MODULE LectureGrid2DModule

CONTAINS

  SUBROUTINE ReadGridFile2D(inp, nxCol, nCol, nVal, at, inv_at, f)

    IMPLICIT NONE

    ! Unité de lecture du fichier input
    INTEGER, intent(in) :: inp
    ! Numéro de colonne où lire la valeur de la position
    INTEGER, intent(in) :: nxCol
    ! Numéro de colonne où lire la valeur de la fonction
    INTEGER, intent(in) :: nCol
    ! Nombre de valeurs dans chaque direction de la fonction à interpoler
    INTEGER, dimension(2), intent(in) :: nVal
    ! Vecteurs de périodicité at(1:2,i) et matrice inverse
    REAL(kind(0.d0)), dimension(1:2,1:2), intent(in) :: at, inv_at
    ! Tableau des valeurs de la fonction à interpoler
    REAL(kind(0.d0)), intent(out), dimension(0:nVal(1)-1,0:nVal(2)-1) :: f

    CHARACTER(len=50), dimension(1:nxCol-1) :: aTemp
    CHARACTER(len=50), dimension(1:nCol-nxCol-2) :: bTemp

    INTEGER :: io, n, m, i
    REAL(kind(0.d0)) :: x, y, fx
    REAL(kind(0.d0)), dimension(1:2) :: ds
    REAL(kind(0.d0)), parameter :: Zero=1.d-4
    LOGICAL, dimension(0:nVal(1)-1,0:nVal(2)-1) :: fRead

    ! Matrice inverse des vecteurs de périodicité
    REAL(kind(0.d0)), dimension(2,2) :: long_inv_at

    long_inv_at(:,1) = dble(nVal(1))*inv_at(:,1)
    long_inv_at(:,2) = dble(nVal(2))*inv_at(:,2)

    fRead(:,:)=.FALSE.
    f(:,:)=0.d0

    io = 0
    DO WHILE (io .EQ. 0)
           CALL Comment(inp)
           READ(inp, *, iostat=io) aTemp(1:nxCol-1), x, y, bTemp(1:nCol-nxCol-2), fx
           ds(:) = MatMul( long_inv_at(:,:), (/ x, y /) )
           n=nInt(ds(1)) ; m=nInt(ds(2))
           IF ( (Abs(dble(n)*at(1,1)+dble(m)*at(1,2)-dble(nVal(1))*x) .GT. Zero*dble(nVal(1))) &
                .OR. (Abs(dble(n)*at(2,1)+dble(m)*at(2,2)-dble(nVal(2))*y) .GT. Zero*dble(nVal(2))) ) THEN
                   WRITE(0,'(a)') "Problème avec la ligne suivante"
                   WRITE(0,*) (Trim(aTemp(i))//" ",i=1,nxCol-1), x, y, (Trim(bTemp(i))//" ",i=1,nCol-nxCol-2), fx
                   WRITE(0, '(2(a,g20.6))') "Pas d'indice calculable pour x = ", x, &
                        "  -  y = ", y
                   WRITE(0,'(a,2g20.6)') "ds(1:2) = ", ds(:)
                   WRITE(0,'(a,2g20.6)') "at(:,:)*ds(:) = ", MatMul(at(:,:),ds(:))
                   WRITE(0,'(2(a,i0))') "n = ", n, "  -  m = ", m
                   WRITE(0,'(a,g20.6)') " n*at(1,1)+m*at(1,2) = ", dble(n)*at(1,1)+dble(m)*at(1,2)
                   WRITE(0,'(a,g20.6)') " n*at(2,1)+m*at(2,2) = ", dble(n)*at(2,1)+dble(m)*at(2,2)
                   STOP
           END IF
           n = Modulo(n, nVal(1))
           m = Modulo(m, nVal(2))
           f(n,m) = fx
           fRead(n,m) = .TRUE.
    END DO

    DO n=0, nVal(1)-1  ;  DO m=0, nVal(2)-1
       IF (fRead(n,m).EQV..False.) THEN
               WRITE(0,'(2(a,i0))') "Pas de valeur trouvée pour les indices n = ", n, &
                " et m = ", m
               STOP
       END IF 
    END DO ; END DO

    !! ==== DEBUG  ============================
    !DO n=0, nVal(1)-1
        !DO m=0, nVal(2)-1
            !r(1:2) = dble(n)/dble(nVal(1))*at(1:2,1) + dble(m)/dble(nVal(2))*at(1:2,2)
            !WRITE(6,'(2i5,3g26.12)') n, m, r(:), f(n,m)
        !END DO
        !WRITE(6,*)
    !END DO
    !! ==== End of DEBUG  =====================

  END SUBROUTINE ReadGridFile2D

END MODULE LectureGrid2DModule
