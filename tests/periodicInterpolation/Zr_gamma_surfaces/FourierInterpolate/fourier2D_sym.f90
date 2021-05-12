MODULE Fourier2DSymModule

       USE Symmetry2DModule
       IMPLICIT NONE
      
       CONTAINS

       SUBROUTINE nFC2D_sym(nFCsym, at, inv_at, sym2D, nSym, nFC)
         ! Calcule le nombre de coefficients de la série de Fourier symmétrisée 
         ! réécrite sous sa forme régulière

         IMPLICIT NONE
         ! Nombre de coefficients de Fourier de la série symétrisée
         INTEGER, dimension(2), intent(in) :: nFCsym
         ! Vecteurs de périodicité at(1:2,i) et matrice inverse
         REAL(kind(0.d0)), dimension(1:2,1:2), intent(in) :: at, inv_at
         ! Opérations de symétrie
         TYPE(sym2Dtype), dimension(1:max_nSym), intent(in) :: sym2D
         ! Nombre d'opérations de symétrie
         INTEGER, intent(in) :: nSym
         ! Nombre de coefficients de Fourier de la série régulière
         INTEGER, dimension(2), intent(out) :: nFC

         INTEGER :: k, l, ns
         INTEGER, dimension(1:2) :: mn
         REAL(kind(0.d0)), dimension(1:2,1:2) :: RotReduced

         nFC(:) = nFCsym

         DO ns=1, nSym
            RotReduced(:,:) = MatMul(inv_at(:,:), MatMul( sym2D(ns)%rot(:,:) , at(:,:) ) )
            DO k=-nFCsym(1)+1, nFCsym(1)-1
               DO l=-nFCsym(2)+1, nFCsym(2)-1
                  mn(:) = k*NInt( RotReduced(1,:) ) + l*NInt( RotReduced(2,:) )
                  IF (mn(1).GT.nFC(1)-1) nFC(1)=mn(1)+1
                  IF (mn(2).GT.nFC(2)-1) nFC(2)=mn(2)+1
                  !WRITE(6,'(5(a,i0))') 'ns=', ns, ' k=', k, ' l=', l, ' m=', mn(1), ' n=', mn(2)
               END DO
            END DO
         END DO

       END SUBROUTINE nFC2D_sym

       SUBROUTINE MinimizeCostFunction2D_sym(fVal, xVal, nVal, FnSym, nFCsym, inv_at, sym2D, nSym)
         ! Ajustement des coefficients de la série de Fourier par minimisation
         ! de la fonction coût sans inversion du système linéaire
         ! Les coefficients complexes Fn sont données en sortie pour la forme
         ! symétrisée de la série de Fourier

         IMPLICIT NONE
         ! Nombre de valeurs de f
         INTEGER, intent(in) :: nVal
         ! fVal(x) valeur de la fonction en x (valeurs discrètes)
         REAL(kind(0.d0)), dimension(1:2,1:nVal), intent(in) :: xVal
         REAL(kind(0.d0)), dimension(1:nVal), intent(in) :: fVal
         ! Nombre de coefficients de Fourier pour la fonction symmétrisée
         INTEGER, dimension(2), intent(in) :: nFCsym
         ! Coefficients complexes de la série de Fourier symmétrisée
         COMPLEX(kind(0.d0)), dimension( -nFCsym(1)+1:nFCsym(1)-1, -nFCsym(2)+1:nFCsym(2)-1 ), intent(inout) :: FnSym
         ! Vecteurs de périodicité at(1:2,i) et matrice inverse
         REAL(kind(0.d0)), dimension(1:2,1:2), intent(in) :: inv_at
         ! Opérations de symétrie
         TYPE(sym2Dtype), dimension(1:max_nSym), intent(in) :: sym2D
         ! Nombre d'opérations de symétrie
         INTEGER, intent(in) :: nSym

         INTEGER :: n, nMax
         REAL(kind(0.d0)) :: fCost0, fCost1, dt, inv_mass, power, ftol
         COMPLEX(kind(0.d0)), dimension( -nFCsym(1)+1:nFCsym(1)-1, -nFCsym(2)+1:nFCsym(2)-1 ) :: FnSym0, FnSym1, FnSym2, dFnSym
         COMPLEX(kind(0.d0)), dimension( -nFCsym(1)+1:nFCsym(1)-1, -nFCsym(2)+1:nFCsym(2)-1 ) :: dfCost

         ! Initialisation des paramètres
         inv_mass = 1.d0/Dble(200*nFCsym(1)*nFCsym(2))
         dt = 1.d0
         nMax = 1000
         ftol = 1.d-4

         ! Valeurs initiales
         FnSym0(:,:) = FnSym(:,:)
         CALL CostFunction2D_sym(fVal, xVal, nVal, FnSym0, nFCsym, inv_at, sym2D, nSym, fCost0, dfCost)
         WRITE(6,'(a,g13.6)') 'n = 0 :  fcost = ', fCost0

         FnSym1(:,:) = FnSym0(:,:) - dt**2*inv_mass*dfCost(:,:)
         dFnSym = Cmplx( 0.d0, 0.d0 )
         CALL CostFunction2D_sym(fVal, xVal, nVal, FnSym1, nFCsym, inv_at, sym2D, nSym, fCost1, dfCost)
         WRITE(6,'(a,g13.6)') 'n = 1 :  fcost = ', fCost1

         n=1
         DO 
                 n = n + 1
                 FnSym2(:,:) = FnSym1(:,:) - dt**2*inv_mass*dfCost(:,:) + dt*dFnSym(:,:)
                 dFnSym = 0.5d0*( FnSym2 - FnSym0 )
                 FnSym0 = FnSym1 ; fCost0 = fCost1
                 FnSym1 = FnSym2
                 power = Dble( Sum( dfCost(:,:)*Conjg(dFnSym(:,:)) ))
                 
                 CALL CostFunction2D_sym(fVal, xVal, nVal, FnSym1, nFCsym, inv_at, sym2D, nSym, fCost1, dfCost)

                 IF (power.GT.0.d0) THEN
                         dFnSym(:,:) = Cmplx(0.d0,0.d0)
                         WRITE(6,'(a,i0,a,g13.6,a)') 'n = ', n, ' :  fcost = ', fCost1,  ' power < 0'
                 ELSE
                         WRITE(6,'(a,i0,a,g13.6)') 'n = ', n, ' :  fcost = ', fCost1
                 END IF

                 IF (Abs(fCost1-fCost0).LE.ftol*Abs(fCost0)) EXIT
                 IF (n.GT.nMax) EXIT
         END DO

         FnSym = FnSym2

       END SUBROUTINE MinimizeCostFunction2D_sym

       SUBROUTINE CostFunction2D_sym(fVal, xVal, nVal, FnSym, nFCsym, inv_at, sym2D, nSym, fCost, dfCost)
         ! Fonction coût à minimiser pour ajuster la série de Fourier
         ! Les coefficients complexes Fn données en entrée correspondent 
         ! à la forme symétrisée de la série de Fourier
         ! Résultats: fcost: fonction coût
         !           dfCost: dérivée par rapport aux coefficients Fn

         IMPLICIT NONE
         ! Nombre de valeurs de f
         INTEGER, intent(in) :: nVal
         ! fVal(x) valeur de la fonction en x (valeurs discrètes)
         REAL(kind(0.d0)), dimension(1:2,1:nVal), intent(in) :: xVal
         REAL(kind(0.d0)), dimension(1:nVal), intent(in) :: fVal
         ! Nombre de coefficients de Fourier pour la fonction symmétrisée
         INTEGER, dimension(2), intent(in) :: nFCsym
         ! Coefficients complexes de la série de Fourier symmétrisée
         COMPLEX(kind(0.d0)), dimension( -nFCsym(1)+1:nFCsym(1)-1, -nFCsym(2)+1:nFCsym(2)-1 ), intent(in) :: FnSym
         ! Vecteurs de périodicité at(1:2,i) et matrice inverse
         REAL(kind(0.d0)), dimension(1:2,1:2), intent(in) :: inv_at
         ! Opérations de symétrie
         TYPE(sym2Dtype), dimension(1:max_nSym), intent(in) :: sym2D
         ! Nombre d'opérations de symétrie
         INTEGER, intent(in) :: nSym
         ! Valeur de la fonction cout
         REAL(kind(0.d0)), intent(out) :: fCost
         ! Dérivée de la fonction coût 
         COMPLEX(kind(0.d0)), dimension( -nFCsym(1)+1:nFCsym(1)-1, -nFCsym(2)+1:nFCsym(2)-1 ), intent(out) :: dfCost

         REAL(kind(0.d0)), parameter :: pi=3.14159265358979d0
         INTEGER :: n, k, l, kl, ns
         REAL(kind(0.d0)), dimension(1:2) :: u
         REAL(kind(0.d0)) :: fx
         COMPLEX(kind(0.d0)) :: e1, e2
         COMPLEX(kind(0.d0)), dimension( 1:nVal, (2*nFCsym(1)-1)*(2*nFCsym(2)-1) ) :: ekl
         COMPLEX(kind(0.d0)), dimension( (2*nFCsym(1)-1)*(2*nFCsym(2)-1), (2*nFCsym(1)-1)*(2*nFCsym(2)-1) ) :: Amat
         COMPLEX(kind(0.d0)), dimension( (2*nFCsym(1)-1)*(2*nFCsym(2)-1) ) :: Bvec, FnVec, dfCostVec

         ! Calcul de la fonction coût
         fCost = 0.d0
         DO n=1, nVal
           CALL TF2DSymCmplx(xVal(:,n), FnSym, nFCsym, inv_at, sym2D, nSym, fx)
           fCost = fCost + ( fVal(n) - fx )**2
         END DO

         ! Calcul de la dérivée
         dfCost(:,:) = Cmplx( 0.d0, 0.d0 )

         ! Construit le système linéaire à résoudre
         ekl(:,:) = cmplx( 0.d0, 0.d0 )
         DO n=1, nVal
            DO ns=1, nSym
               u(:) = 2.d0*pi*MatMul(inv_at(:,:), & 
                       MatMul( sym2D(ns)%rot(:,:), xVal(:,n) ) + sym2D(ns)%u(:) )
               kl=0
               DO k=-nFCsym(1)+1, nFCsym(1)-1
                  e1 = Cmplx( cos(dble(k)*U(1)), sin(dble(k)*U(1)) )
                  DO l=-nFCsym(2)+1, nFCsym(2)-1
                      e2 = Cmplx( cos(dble(l)*U(2)), sin(dble(l)*U(2)) )
                      kl=kl+1
                      ekl(n,kl) = ekl(n,kl) + e1*e2 
                  END DO
               END DO
            END DO
         END DO
         ekl(:,:) = 1.d0/dble(nSym)*ekl(:,:)

         ! Les coefficients complexes de la série de Fourier sont solutions de
         ! Amat(:,:)*X(:) = Bvec(:)
         Amat(:,:) = MatMul( Transpose( Conjg(ekl(:,:)) ), ekl(:,:) )
         Bvec(:) = MatMul( Transpose( Conjg(ekl(:,:)) ), Cmplx( fVal(:), 0.d0) )

         ! Écriture sous forme d'un vecteur du tableau contenant les coef de Fourier
         kl=0
         DO k=-nFCsym(1)+1, nFCsym(1)-1
            DO l=-nFCsym(2)+1, nFCsym(2)-1
                kl=kl+1
                FnVec(kl) = FnSym(k,l)
            END DO
         END DO

         ! Vecteur dérivée
         dfCostVec(:) = MatMul( Amat(:,:), FnVec(:) ) - BVec(:)

         ! Réécriture sous forme d'un tableau du vecteur dérivée
         kl=0
         DO k=-nFCsym(1)+1, nFCsym(1)-1
            DO l=-nFCsym(2)+1, nFCsym(2)-1
                kl=kl+1
                dfCost(k,l) = dfCostVec(kl)
            END DO
         END DO


       END SUBROUTINE CostFunction2D_sym

       SUBROUTINE FitFourier2DSymCmplx(fVal, xVal, nVal, FnSym, nFCsym, inv_at, sym2D, nSym)
         ! Ajustement moindre carré d'une série de Fourier à partir de nVal
         ! points qui ne sont pas forcément répartis sur une grille régulière
         ! en imposant les opérations de symétrie sym2D
         ! Les coefficients complexes Fn sont données en sortie pour la forme
         ! symétrisée de la série de Fourier

         IMPLICIT NONE
         ! Nombre de valeurs de f
         INTEGER, intent(in) :: nVal
         ! fVal(x) valeur de la fonction en x (valeurs discrètes)
         REAL(kind(0.d0)), dimension(1:2,1:nVal), intent(in) :: xVal
         REAL(kind(0.d0)), dimension(1:nVal), intent(in) :: fVal
         ! Nombre de coefficients de Fourier pour la fonction symmétrisée
         INTEGER, dimension(2), intent(in) :: nFCsym
         ! Coefficients complexes de la série de Fourier symmétrisée
         COMPLEX(kind(0.d0)), dimension( -nFCsym(1)+1:nFCsym(1)-1, -nFCsym(2)+1:nFCsym(2)-1 ), intent(out) :: FnSym
         ! Vecteurs de périodicité at(1:2,i) et matrice inverse
         REAL(kind(0.d0)), dimension(1:2,1:2), intent(in) :: inv_at
         ! Opérations de symétrie
         TYPE(sym2Dtype), dimension(1:max_nSym), intent(in) :: sym2D
         ! Nombre d'opérations de symétrie
         INTEGER, intent(in) :: nSym

         REAL(kind(0.d0)), parameter :: pi=3.14159265358979d0
         INTEGER :: n, k, l, kl, nDim, info, ns
         REAL(kind(0.d0)), dimension(1:2) :: u
         REAL(kind(0.d0)) :: f0
         COMPLEX(kind(0.d0)) :: e1, e2
         COMPLEX(kind(0.d0)), dimension( 1:nVal, (2*nFCsym(1)-1)*(2*nFCsym(2)-1) ) :: ekl
         COMPLEX(kind(0.d0)), dimension( (2*nFCsym(1)-1)*(2*nFCsym(2)-1), (2*nFCsym(1)-1)*(2*nFCsym(2)-1) ) :: Amat
         COMPLEX(kind(0.d0)), dimension( (2*nFCsym(1)-1)*(2*nFCsym(2)-1) ) :: Bvec
         INTEGER, dimension( (2*nFCsym(1)-1)*(2*nFCsym(2)-1) ) :: iPiv

         ! Calcul de la valeur moyenne
         f0 = Sum( fVal(1:nVal) )/dble(nVal)

         ! Nombre de coefficients complexes à trouver
         nDim = (2*nFCsym(1)-1)*(2*nFCsym(2)-1)

         ! Construit le système linéaire à résoudre
         ekl(:,:) = cmplx( 0.d0, 0.d0 )
         DO n=1, nVal
            DO ns=1, nSym
               u(:) = 2.d0*pi*MatMul(inv_at(:,:), & 
                       MatMul( sym2D(ns)%rot(:,:), xVal(:,n) ) + sym2D(ns)%u(:) )
               kl=0
               DO k=-nFCsym(1)+1, nFCsym(1)-1
                  e1 = Cmplx( cos(dble(k)*U(1)), sin(dble(k)*U(1)) )
                  DO l=-nFCsym(2)+1, nFCsym(2)-1
                      e2 = Cmplx( cos(dble(l)*U(2)), sin(dble(l)*U(2)) )
                      kl=kl+1
                      ekl(n,kl) = ekl(n,kl) + e1*e2 
                  END DO
               END DO
            END DO
         END DO
         ekl(:,:) = 1.d0/dble(nSym)*ekl(:,:)

         ! Les coefficients complexes de la série de Fourier sont solutions de
         ! Amat(:,:)*X(:) = Bvec(:)
         Amat(:,:) = MatMul( Transpose( Conjg(ekl(:,:)) ), ekl(:,:) )
         Bvec(:) = MatMul( Transpose( Conjg(ekl(:,:)) ), Cmplx( fVal(:)-f0, 0.d0) )
         CALL zgesv(nDim, 1, Amat, nDim, iPiv, Bvec, nDim, info)
         IF (info.NE.0) THEN
                 WRITE(0,'(a,i0)') 'zsysv: Info = ', info
                 STOP '< FitFourier >'
         END IF

         ! Bvec(:) contient les coefficients complexes de la série de Fourier
         ! symmétrisée
         !   => réécriture sous forme d'un tableau
         kl=0
         DO k=-nFCsym(1)+1, nFCsym(1)-1
            DO l=-nFCsym(2)+1, nFCsym(2)-1
                kl=kl+1
                FnSym(k,l) = Bvec(kl)
            END DO
         END DO

         ! On rajoute la valeur moyenne
         FnSym(0,0) = FnSym(0,0) + Cmplx(f0,0.d0)

       END SUBROUTINE FitFourier2DSymCmplx

       SUBROUTINE coefFourier2D_unsymmetrize(FnSym, nFCsym, Fn, nFC, at, inv_at, sym2D, nSym)
         ! Construit les coefficients de la série de Fourier régulière 
         ! à partir des coefficients de la sérier symétrisée      

         IMPLICIT NONE
         ! Nombre de coefficients de Fourier pour la fonction symétrisée
         INTEGER, dimension(2), intent(in) :: nFCsym
         ! Coefficients complexes de la série de Fourier symétrisée
         COMPLEX(kind(0.d0)), dimension( -nFCsym(1)+1:nFCsym(1)-1, -nFCsym(2)+1:nFCsym(2)-1 ), intent(in) :: FnSym
         ! Nombre de coefficients de Fourier pour la série standarde
         INTEGER, dimension(2), intent(in) :: nFC
         ! Coefficients complexes de la série de Fourier standarde
         COMPLEX(kind(0.d0)), dimension( -nFC(1)+1:nFC(1)-1, -nFC(2)+1:nFC(2)-1 ), intent(out) :: Fn
         ! Vecteurs de périodicité at(1:2,i) et matrice inverse
         REAL(kind(0.d0)), dimension(1:2,1:2), intent(in) :: at, inv_at
         ! Opérations de symétrie
         TYPE(sym2Dtype), dimension(1:max_nSym), intent(in) :: sym2D
         ! Nombre d'opérations de symétrie
         INTEGER, intent(in) :: nSym

         REAL(kind(0.d0)), parameter :: pi=3.14159265358979d0
         INTEGER :: m, n, k, l, ns
         INTEGER, dimension(1:2) :: mn
         COMPLEX(kind(0.d0)) :: e1, e2
         REAL(kind(0.d0)), dimension(1:2,1:2) :: RotReduced
         REAL(kind(0.d0)), dimension(1:2) :: u

         ! Réécriture de la série de Fourier symmétrisée en série de Fourier
         ! standard
         Fn(:,:) = Cmplx( 0.d0, 0.d0 )
         DO ns=1, nSym
            u(:) = 2.d0*pi*MatMul(inv_at(:,:), sym2D(ns)%u(:) )
            RotReduced(:,:) = MatMul(inv_at(:,:), MatMul( sym2D(ns)%rot(:,:) , at(:,:) ) )
            DO k=-nFCsym(1)+1, nFCsym(1)-1
               e1 = Cmplx( cos(dble(k)*U(1)), sin(dble(k)*U(1)) )
               DO l=-nFCsym(2)+1, nFCsym(2)-1
                  e2 = Cmplx( cos(dble(l)*U(2)), sin(dble(l)*U(2)) )
                  mn(:) = k*NInt( RotReduced(1,:) ) + l*NInt( RotReduced(2,:) )
                  m = mn(1)
                  n = mn(2)
                  Fn(m,n) = Fn(m,n) + FnSym(k,l)*e1*e2
               END DO
            END DO
         END DO
         Fn(:,:) = Fn(:,:)/dble(nSym)

       END SUBROUTINE coefFourier2D_unsymmetrize

       SUBROUTINE FitFourier2D_sym(fVal, xVal, nVal, nFCsym, An, Bn, Cn, Dn, nFC, at, inv_at, sym2D, nSym)
         ! Ajustement moindre carré d'une série de Fourier à partir de nVal
         ! points qui ne sont pas forcément répartis sur une grille régulière
         ! en imposant les opérations de symétrie sym2D

         USE Fourier2DModule
         IMPLICIT NONE
         ! Nombre de valeurs de f
         INTEGER, intent(in) :: nVal
         ! fVal(x) valeur de la fonction en x (valeurs discrètes)
         REAL(kind(0.d0)), dimension(1:2,1:nVal), intent(in) :: xVal
         REAL(kind(0.d0)), dimension(1:nVal), intent(in) :: fVal
         ! Nombre de coefficients de Fourier pour la fonction symmétrisée
         INTEGER, dimension(2), intent(in) :: nFCsym
         ! Nombre de coefficients de Fourier pour la série standarde
         INTEGER, dimension(2), intent(in) :: nFC
         ! Coefficients de la série de Fourier
         REAL(kind(0.d0)), dimension(0:nFC(1)-1,0:nFC(2)-1), intent(out) :: An, Bn, Cn, Dn
         ! Vecteurs de périodicité at(1:2,i) et matrice inverse
         REAL(kind(0.d0)), dimension(1:2,1:2), intent(in) :: at, inv_at
         ! Opérations de symétrie
         TYPE(sym2Dtype), dimension(1:max_nSym), intent(in) :: sym2D
         ! Nombre d'opérations de symétrie
         INTEGER, intent(in) :: nSym

         COMPLEX(kind(0.d0)), dimension( -nFCsym(1)+1:nFCsym(1)-1, -nFCsym(2)+1:nFCsym(2)-1 ) :: FnSym
         COMPLEX(kind(0.d0)), dimension( -nFC(1)+1:nFC(1)-1, -nFC(2)+1:nFC(2)-1 ) :: Fn

         CALL FitFourier2DSymCmplx(fVal, xVal, nVal, FnSym, nFCsym, inv_at, sym2D, nSym)
         CALL coefFourier2D_unsymmetrize(FnSym, nFCsym, Fn, nFC, at, inv_at, sym2D, nSym)
         CALL coefFourier2D_cmplx2real(Fn, nFC, An, Bn, Cn, Dn)

        END SUBROUTINE FitFourier2D_sym

        SUBROUTINE TF2DSymCmplx(x, FnSym, nFCsym, inv_at, sym2D, nSym, f, df, d2f)
         ! Calcul de la fonction f en un point x au moyen des séries de
         ! Fourier, ainsi que ses dérivées premières et secondes
         ! en utilisant les coefficients complexes 
         ! et en utilisant la version symétrisée de la série

         IMPLICIT NONE

         ! x: point auquel on calcule f
         REAL(kind(0.d0)), dimension(2), intent(in) :: x
         ! Nombre de coefficients de Fourier
         INTEGER, dimension(2), intent(in) :: nFCsym
         ! Coefficients de Fourier complexes
         COMPLEX(kind(0.d0)), dimension( -nFCsym(1)+1:nFCsym(1)-1, -nFCsym(2)+1:nFCsym(2)-1 ), intent(in) :: FnSym
         ! Matrice inverse des vecteurs de périodicité at(1:2,i)
         REAL(kind(0.d0)), dimension(1:2,1:2), intent(in) :: inv_at
         ! Opérations de symétrie
         TYPE(sym2Dtype), dimension(1:max_nSym), intent(in) :: sym2D
         ! Nombre d'opérations de symétrie
         INTEGER, intent(in) :: nSym
         ! Valeur de la fonction en x 
         REAL(kind(0.d0)), intent(out) :: f
         ! Derivée première
         REAL(kind(0.d0)), dimension(1:2), intent(out), optional :: df
         ! Derivée seconde
         REAL(kind(0.d0)), dimension(1:2,1:2), intent(out), optional :: d2f

         ! Fréquence de la fonction
         REAL(kind(0.d0)), dimension(1:2) :: U
         REAL(kind(0.d0)), parameter :: pi=3.14159265358979d0
         COMPLEX(kind(0.d0)) :: e1, e2
         REAL(kind(0.d0)) :: f_kl
         REAL(kind(0.d0)), dimension(1:2) :: dfdu, df_kl
         REAL(kind(0.d0)), dimension(1:2,1:2) :: g, d2fdu2, d2f_kl
         INTEGER :: ns, k, l
         LOGICAL :: calcul_df, calcul_d2f
         
         calcul_df =  Present(df) 
         calcul_d2f = Present(d2f)


         f = 0.d0
         IF (calcul_df) df(:) = 0.d0
         IF (calcul_d2f) d2f(:,:) = 0.d0

         DO k = -nFCsym(1)+1, nFCsym(1)-1
            DO l = -nFCsym(2)+1, nFCsym(2)-1
               f_kl = 0.d0
               df_kl(:) = 0.d0
               d2f_kl(:,:) = 0.d0   
               symmetry_loop: DO ns=1, nSym
                  U(:) = 2.d0*pi*MatMul(inv_at(:,:), &
                             MatMul( sym2D(ns)%rot(:,:), x(:) ) + sym2D(ns)%u(:) )
                  e1 = Cmplx( cos(dble(k)*U(1)), sin(dble(k)*U(1)) )
                  e2 = Cmplx( cos(dble(l)*U(2)), sin(dble(l)*U(2)) )
                  f_kl = f_kl + Dble( FnSym(k,l)*e1*e2 )
                  IF (calcul_df) THEN
                           dfdU(1) = dfdU(1) + Dble(k)*aImag( FnSym(k,l)*e1*e2 )
                           dfdU(2) = dfdU(2) + Dble(l)*aImag( FnSym(k,l)*e1*e2 )
                           g(:,:) = 2.d0*pi*Transpose(MatMul( inv_at(:,:), sym2D(ns)%rot(:,:) ))
                           df_kl(:) = df_kl(:) + MatMul(g(:,:) , dfdU(:) )
                  END IF
                  IF (calcul_d2f) THEN
                           d2fdU2(1,1) = d2fdU2(1,1) - Dble(k*k)*Dble( FnSym(k,l)*e1*e2 )
                           d2fdU2(1,2) = d2fdU2(1,2) + Dble(k*l)*Dble( FnSym(k,l)*e1*e2 )
                           d2fdU2(2,1) = d2fdU2(1,2)
                           d2fdU2(2,2) = d2fdU2(2,2) + Dble(l*l)*Dble( FnSym(k,l)*e1*e2 )
                           g(:,:) = 2.d0*pi*Transpose(MatMul( inv_at(:,:), sym2D(ns)%rot(:,:) ))
                           d2f_kl(:,:) = d2f_kl(:,:) + MatMul(  MatMul( g(:,:), g(:,:) ), d2fDU2(:,:) )
                  END IF
               END DO symmetry_loop
               f = f + f_kl/dble(nSym)
               IF (calcul_df) df(:) = df(:) + df_kl(:)/dble(nSym)
               IF (calcul_d2f) d2f(:,:) = d2f(:,:) + d2f_kl(:,:)/dble(nSym)
            END DO
         END DO

        END SUBROUTINE TF2DSymCmplx

        !!$ SUBROUTINE TF2DSymCmplx(x, FnSym, nFCsym, inv_at, sym2D, nSym, f, df, d2f)
        !!$  ! Calcul de la fonction f en un point x au moyen des séries de
        !!$  ! Fourier, ainsi que ses dérivées premières et secondes
        !!$  ! en utilisant les coefficients complexes 
        !!$  ! et en utilisant la version symétrisée de la série

        !!$  IMPLICIT NONE

        !!$  ! x: point auquel on calcule f
        !!$  REAL(kind(0.d0)), dimension(2), intent(in) :: x
        !!$  ! Nombre de coefficients de Fourier
        !!$  INTEGER, dimension(2), intent(in) :: nFCsym
        !!$  ! Coefficients de Fourier complexes
        !!$  COMPLEX(kind(0.d0)), dimension( -nFCsym(1)+1:nFCsym(1)-1, -nFCsym(2)+1:nFCsym(2)-1 ), intent(in) :: FnSym
        !!$  ! Matrice inverse des vecteurs de périodicité at(1:2,i)
        !!$  REAL(kind(0.d0)), dimension(1:2,1:2), intent(in) :: inv_at
        !!$  ! Opérations de symétrie
        !!$  TYPE(sym2Dtype), dimension(1:max_nSym), intent(in) :: sym2D
        !!$  ! Nombre d'opérations de symétrie
        !!$  INTEGER, intent(in) :: nSym
        !!$  ! Valeur de la fonction en x 
        !!$  REAL(kind(0.d0)), intent(out) :: f
        !!$  ! Derivée première
        !!$  REAL(kind(0.d0)), dimension(1:2), intent(out), optional :: df
        !!$  ! Derivée seconde
        !!$  REAL(kind(0.d0)), dimension(1:2,1:2), intent(out), optional :: d2f

        !!$  ! Fréquence de la fonction
        !!$  REAL(kind(0.d0)), dimension(1:2) :: U
        !!$  REAL(kind(0.d0)), parameter :: pi=3.14159265358979d0
        !!$  COMPLEX(kind(0.d0)) :: e1, e2
        !!$  REAL(kind(0.d0)), dimension(1:2) :: dfdu
        !!$  REAL(kind(0.d0)), dimension(1:2,1:2) :: g, d2fdu2
        !!$  INTEGER :: ns, k, l
        !!$  LOGICAL :: calcul_df, calcul_d2f
        !!$  
        !!$  calcul_df =  Present(df) 
        !!$  calcul_d2f = Present(d2f)


        !!$  f = 0.d0
        !!$  IF (calcul_df) df(:) = 0.d0
        !!$  IF (calcul_d2f) d2f(:,:) = 0.d0

        !!$  DO ns=1, nSym
        !!$     U(:) = 2.d0*pi*MatMul(inv_at(:,:), &
        !!$                MatMul( sym2D(ns)%rot(:,:), x(:) ) + sym2D(ns)%u(:) )
        !!$     dfdU(:) = 0.d0
        !!$     d2fdU2(:,:) = 0.d0   
        !!$     DO k = -nFCsym(1)+1, nFCsym(1)-1
        !!$        e1 = Cmplx( cos(dble(k)*U(1)), sin(dble(k)*U(1)) )
        !!$        DO l = -nFCsym(2)+1, nFCsym(2)-1
        !!$            e2 = Cmplx( cos(dble(l)*U(2)), sin(dble(l)*U(2)) )
        !!$            f = f + Dble( FnSym(k,l)*e1*e2 )
        !!$            IF (calcul_df) THEN
        !!$                    dfdU(1) = dfdU(1) + Dble(k)*aImag( FnSym(k,l)*e1*e2 )
        !!$                    dfdU(2) = dfdU(2) + Dble(l)*aImag( FnSym(k,l)*e1*e2 )
        !!$            END IF
        !!$            IF (calcul_d2f) THEN
        !!$                    d2fdU2(1,1) = d2fdU2(1,1) - Dble(k*k)*Dble( FnSym(k,l)*e1*e2 )
        !!$                    d2fdU2(1,2) = d2fdU2(1,2) + Dble(k*l)*Dble( FnSym(k,l)*e1*e2 )
        !!$                    d2fdU2(2,2) = d2fdU2(2,2) + Dble(l*l)*Dble( FnSym(k,l)*e1*e2 )
        !!$            END IF
        !!$        END DO
        !!$     END DO
        !!$     IF (calcul_df.OR.calcul_d2f) &
        !!$             g(:,:) = 2.d0*pi/dble(nSym)*Transpose(MatMul( inv_at(:,:), sym2D(ns)%rot(:,:) ))
        !!$     IF (calcul_df) &
        !!$             df(:) = MatMul(g(:,:) , dfdU(:) )
        !!$     IF (calcul_d2f) THEN
        !!$             d2fdU2(2,1) = d2fdU2(1,2)
        !!$             d2f(:,:) = MatMul(  MatMul( g(:,:), g(:,:) ), d2fDU2(:,:) )
        !!$     END IF
        !!$  END DO
        !!$  f = f/dble(nSym)

        !!$ END SUBROUTINE TF2DSymCmplx

      END MODULE Fourier2DSymModule
