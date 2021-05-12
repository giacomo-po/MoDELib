MODULE MinModule


CONTAINS

  SUBROUTINE Minimize(f, f0, x0, dxMin, dxMax, nIterMax)

     IMPLICIT NONE
     ! Fonction à minimiser
     EXTERNAL :: f
     ! Valeur du minimum
     REAL(kind(0.d0)), intent(out) :: f0
     ! Variable: en entrée, point de départ
     !           en sortie, position du minimum
     REAL(kind(0.d0)), intent(inout) :: x0
     ! Tolérance sur x
     REAL(kind(0.d0)), intent(in) :: dxMin
     ! Variation maximale de x sur une itération
     REAL(kind(0.d0)), intent(in) :: dxMax
     ! Maximal number of iterations
     INTEGER, intent(in) :: nIterMax

     REAL(kind(0.d0)) :: dx, x1, f1, df0, d2f0
     INTEGER :: n

     IF (dxMin.LE.0.d0)  STOP '< Minimize > : dxMin < 0'
     IF (dxMax.LE.dxMin) STOP '< Minimize > : dxMax < dxMin'

     ! Valeur de la fonction au point x0
     CALL f(x0, f0, df0, d2f0)
     !WRITE(6,'(i6,1x,4g20.12)') 0, x0, f0, df0, d2f0    ! DEBUG 

     DO n=1, nIterMax

             ! Variation de x
             IF ( Abs(df0) .GT. dxMax*Abs(d2f0) ) THEN
                     dx = -Sign(dxmax,df0)
             ELSE
                     dx = -df0/d2f0
             END IF

             ! Nouvelle position
             x1 = x0 + dx
             CALL f(x1, f1, df0, d2f0)
             
             ! Vérification qu'on se déplace bien vers un minimum
             IF (f1.GT.f0) THEN ! vers un maximum
                     x0 = x0 - dx
                     CALL f(x0, f0, df0, d2f0)
             ELSE               ! vers un minimum
                     x0 = x1 ; f0 = f1
                     ! Test de sortie
                     IF (Abs(dx).LE.dxMin) Return
             END IF

             !WRITE(6,'(i6,1x,4g20.12)') n, x0, f0, df0, d2f0    ! DEBUG 
     END DO

     IF (n.GE.nIterMax) THEN
             WRITE(0,'(a)') "< Minimize > : nombre maximal d'itérations atteint"
             WRITE(0,'(i6,1x,4g20.12)') n, x0, f0, df0, d2f0 
     END IF

     


  END SUBROUTINE Minimize


END MODULE MinModule
