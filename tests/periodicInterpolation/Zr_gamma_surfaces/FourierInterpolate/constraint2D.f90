MODULE Constraint2DModule

        SAVE

        ! Vecteurs de périodicité de la fonction à interpoler: at(1:2,i)
        !   et matrice inverse
        REAL(kind(0.d0)), dimension(1:2,1:2) :: at, inv_at
        ! Coef de Fourier cos et sin : cf(0:nFC-1) et sf(0:nFC-1)
        INTEGER, dimension(2) :: nFC
        REAL(kind(0.d0)), dimension(:,:), allocatable :: An, Bn, Cn, Dn

        ! Origine et direction de la droite définissant la fonction contrainte
        REAL(kind(0.d0)), dimension(1:2), private :: x0, nDir0

CONTAINS

        SUBROUTINE Init_fConstrained(x, nDir)

           IMPLICIT NONE
           ! Origine et direction de la droite définissant la fonction contrainte
           REAL(kind(0.d0)), dimension(:), intent(in) :: x
           REAL(kind(0.d0)), dimension(:), intent(in) :: nDir

           x0(:) = x(:)
           nDir0(:) = nDir(:)

        END SUBROUTINE Init_fConstrained

        SUBROUTINE fConstrained(t, f, df, d2f)
        
           USE Fourier2DModule
           IMPLICIT NONE
           REAL(kind(0.d0)), intent(in) :: t
           REAL(kind(0.d0)), intent(out) :: f, df, d2f

           REAL(kind(0.d0)), dimension(1:2) :: x, dfdx
           REAL(kind(0.d0)), dimension(1:2,1:2) :: d2fdx2

           ! Point où calculer la fonction
           x(:) = x0(:) + t*nDir0(:)
        
           ! Fonction et dérivées dans l'espace à 2D
           CALL TF2D(x, An, Bn, Cn, Dn, nFC, inv_at, f, dfdx, d2fdx2)

           ! Dérivées première et seconde de la fonction contrainte
           df = Sum( dfdx(:)*nDir0(:) )
           d2f = Dot_Product(nDir0(:), MatMul(d2fdx2(:,:),nDir0(:)) )
           !d2f = nDir0(1)*d2fdx2(1,1)*nDir0(1) &
               !+ nDir0(1)*d2fdx2(1,2)*nDir0(2) &
               !+ nDir0(2)*d2fdx2(2,1)*nDir0(1) &
               !+ nDir0(2)*d2fdx2(2,2)*nDir0(2)
        
        END SUBROUTINE fConstrained

END MODULE Constraint2DModule
