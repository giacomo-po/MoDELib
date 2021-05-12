      MODULE Fourier2DModule

       IMPLICIT NONE
      
       CONTAINS

       SUBROUTINE FitFourier2D(fVal, xVal, nVal, nFC, An, Bn, Cn, Dn, inv_at)
         ! Ajustement moindre carré d'une série de Fourier à partir de nVal
         ! points qui ne sont pas forcément répartis sur une grille régulière

         ! Nombre de valeurs de f
         INTEGER, intent(in) :: nVal
         ! fVal(x) valeur de la fonction en x (valeurs discrètes)
         REAL(kind(0.d0)), dimension(1:2,1:nVal), intent(in) :: xVal
         REAL(kind(0.d0)), dimension(1:nVal), intent(in) :: fVal
         ! Nombre de coefficients de Fourier
         INTEGER, dimension(2), intent(in) :: nFC
         ! Coefficients de la série de Fourier
         REAL(kind(0.d0)), DIMENSION(0:nFC(1)-1,0:nFC(2)-1), intent(out) :: An, Bn, Cn, Dn
         ! Matrice inverse des vecteurs de périodicité at(1:2,i)
         REAL(kind(0.d0)), dimension(1:2,1:2), intent(in) :: inv_at

         COMPLEX(kind(0.d0)), dimension( -nFC(1)+1:nFC(1)-1, -nFC(2)+1:nFC(2)-1 ) :: Fn

         CALL FitFourier2Dcmplx(fVal, xVal, nVal, Fn, nFC, inv_at)
         CALL coefFourier2D_cmplx2real(Fn, nFC, An, Bn, Cn, Dn)

        END SUBROUTINE FitFourier2D

       SUBROUTINE FitFourier2Dcmplx(fVal, xVal, nVal, Fn, nFC, inv_at)
         ! Ajustement moindre carré d'une série de Fourier à partir de nVal
         ! points qui ne sont pas forcément répartis sur une grille régulière

         IMPLICIT NONE
         ! Nombre de valeurs de f
         INTEGER, intent(in) :: nVal
         ! fVal(x) valeur de la fonction en x (valeurs discrètes)
         REAL(kind(0.d0)), dimension(1:2,1:nVal), intent(in) :: xVal
         REAL(kind(0.d0)), dimension(1:nVal), intent(in) :: fVal
         ! Nombre de coefficients de Fourier
         INTEGER, dimension(2), intent(in) :: nFC
         ! Coefficients de la série de Fourier complexes
         COMPLEX(kind(0.d0)), dimension( -nFC(1)+1:nFC(1)-1, -nFC(2)+1:nFC(2)-1 ), intent(out) :: Fn
         ! Matrice inverse des vecteurs de périodicité at(1:2,i)
         REAL(kind(0.d0)), dimension(1:2,1:2), intent(in) :: inv_at

         REAL(kind(0.d0)), parameter :: pi=3.14159265358979d0
         INTEGER :: n, k, l, kl, nDim, info
         REAL(kind(0.d0)), dimension(1:2) :: u
         REAL(kind(0.d0)) :: f0
         COMPLEX(kind(0.d0)) :: e1, e2
         COMPLEX(kind(0.d0)), dimension( 1:nVal, (2*nFC(1)-1)*(2*nFC(2)-1) ) :: ekl
         COMPLEX(kind(0.d0)), dimension( (2*nFC(1)-1)*(2*nFC(2)-1), (2*nFC(1)-1)*(2*nFC(2)-1) ) :: Amat
         COMPLEX(kind(0.d0)), dimension( (2*nFC(1)-1)*(2*nFC(2)-1) ) :: Bvec
         INTEGER, dimension( (2*nFC(1)-1)*(2*nFC(2)-1) ) :: iPiv

         ! Calcul de la valeur moyenne
         f0 = Sum( fVal(1:nVal) )/dble(nVal)

         ! Nombre de coefficients complexes à trouver
         nDim = (2*nFC(1)-1)*(2*nFC(2)-1)

         ! Construit le système linéaire à résoudre
         DO n=1, nVal
            u(:) = 2.d0*pi*MatMul(inv_at(:,:),xVal(:,n))

            kl=0
            DO k=-nFC(1)+1, nFC(1)-1
               e1 = Cmplx( cos(dble(k)*U(1)), sin(dble(k)*U(1)) )
               DO l=-nFC(2)+1, nFC(2)-1
                   e2 = Cmplx( cos(dble(l)*U(2)), sin(dble(l)*U(2)) )
                   kl=kl+1
                   ekl(n,kl) = e1*e2 
               END DO
            END DO

         END DO

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
         !   => réécriture sous forme d'un tableau
         kl=0
         DO k=-nFC(1)+1, nFC(1)-1
            DO l=-nFC(2)+1, nFC(2)-1
                kl=kl+1
                Fn(k,l) = Bvec(kl)
            END DO
         END DO

         ! On rajoute la valeur moyenne
         Fn(0,0) = Fn(0,0) + Cmplx(f0, 0.d0)

        END SUBROUTINE FitFourier2Dcmplx

        SUBROUTINE coefFourier2D_cmplx2real(Fn, nFC, An, Bn, Cn, Dn)
         ! Calcule les coefficients réels de la série de Fourier en cos() et
         ! sin() à partir des coefficients complexes
         ! points qui ne sont pas forcément répartis sur une grille régulière

         IMPLICIT NONE
         ! Nombre de coefficients de Fourier
         INTEGER, dimension(2), intent(in) :: nFC
         ! Coefficients complexes de la série de Fourier
         COMPLEX(kind(0.d0)), dimension( -nFC(1)+1:nFC(1)-1, -nFC(2)+1:nFC(2)-1 ), intent(in) :: Fn
         ! Coefficients réels de la série de Fourier
         REAL(kind(0.d0)), DIMENSION(0:nFC(1)-1,0:nFC(2)-1), intent(out) :: An, Bn, Cn, Dn

         INTEGER :: k, l

         ! Coefficients réels des séries de Fourier en cos() et sin()
         DO k=0, nFC(1)-1
            DO l=0, nFC(2)-1
               An(k,l) =   Dble( Fn(k,l) + Fn(-k,l) + Fn(k,-l) + Fn(-k,-l) )
               Bn(k,l) = -aImag( Fn(k,l) + Fn(-k,l) - Fn(k,-l) - Fn(-k,-l) )
               Cn(k,l) = -aImag( Fn(k,l) - Fn(-k,l) + Fn(k,-l) - Fn(-k,-l) )
               Dn(k,l) = - Dble( Fn(k,l) - Fn(-k,l) - Fn(k,-l) + Fn(-k,-l) )
            END DO
         END DO

         An(0,:) = 0.5d0*An(0,:) 
         Bn(0,:) = 0.5d0*Bn(0,:) 
         Cn(0,:) = 0.5d0*Cn(0,:) 
         Dn(0,:) = 0.5d0*Dn(0,:) 
         An(:,0) = 0.5d0*An(:,0) 
         Bn(:,0) = 0.5d0*Dn(:,0) 
         Cn(:,0) = 0.5d0*Cn(:,0) 
         Dn(:,0) = 0.5d0*Dn(:,0) 

       END SUBROUTINE coefFourier2D_cmplx2real 

        SUBROUTINE coefFourier2D(f, nVal, nFC, An, Bn, Cn, Dn)
         ! Calcul des coefficients d'une série de Fourier à partir d'une grille
         ! régulière

         IMPLICIT NONE
         ! Nombre de valeurs de f
         INTEGER, dimension(2), intent(in) :: nVal
         ! f(x) valeur de la fonction en x (valeurs discrètes)
         REAL(kind(0.d0)), dimension(0:nVal(1)-1,0:nVal(2)-1), intent(in) :: f
         ! Nombre de coefficients de Fourier
         INTEGER, dimension(2), intent(in) :: nFC
         ! Coefficients de la série de Fourier
         REAL(kind(0.d0)), DIMENSION(0:nFC(1)-1,0:nFC(2)-1), intent(out) :: An, Bn, Cn, Dn

         REAL(kind(0.d0)), parameter :: pi=3.14159265358979d0
         REAL(kind(0.d0)) :: omega1, omega2, f0
         INTEGER :: n, m, k, l

         omega1=2.d0*pi/dble(nVal(1))
         omega2=2.d0*pi/dble(nVal(2))

         f0=f(0,0)
         DO k=0, nFC(1)-1
            DO l=0, nFC(2)-1
               An(k,l)=0.d0 ; Bn(k,l)=0.d0 ; Cn(k,l)=0.d0 ; Dn(k,l)=0.d0
               DO n=0, nVal(1)-1
                  DO m=0, nVal(2)-1
                     An(k,l) = An(k,l) &
                        + (f(n,m)-f0)*Cos(dble(k*n)*omega1)*Cos(dble(l*m)*omega2)
                     Bn(k,l) = Bn(k,l) &
                        + (f(n,m)-f0)*Cos(dble(k*n)*omega1)*Sin(dble(l*m)*omega2)
                     Cn(k,l) = Cn(k,l) &
                        + (f(n,m)-f0)*Sin(dble(k*n)*omega1)*Cos(dble(l*m)*omega2)
                     Dn(k,l) = Dn(k,l) &
                        + (f(n,m)-f0)*Sin(dble(k*n)*omega1)*Sin(dble(l*m)*omega2)
                  END DO
               END DO
            END DO
         END DO
         An(:,:) = 4.d0*An(:,:)/dble(nVal(1)*nVal(2))
         Bn(:,:) = 4.d0*Bn(:,:)/dble(nVal(1)*nVal(2))
         Cn(:,:) = 4.d0*Cn(:,:)/dble(nVal(1)*nVal(2))
         Dn(:,:) = 4.d0*Dn(:,:)/dble(nVal(1)*nVal(2))
         An(0,:) = 0.5d0*An(0,:) 
         Bn(0,:) = 0.5d0*Bn(0,:) 
         Cn(0,:) = 0.5d0*Cn(0,:) 
         Dn(0,:) = 0.5d0*Dn(0,:) 
         An(:,0) = 0.5d0*An(:,0) 
         Bn(:,0) = 0.5d0*Dn(:,0) 
         Cn(:,0) = 0.5d0*Cn(:,0) 
         Dn(:,0) = 0.5d0*Dn(:,0) 
         An(0,0) = An(0,0)+f0

        END SUBROUTINE coefFourier2D

        SUBROUTINE TF2D(x, An, Bn, Cn, Dn, nFC, inv_at, f, df, d2f)
         ! Calcul de la fonction f en un point x au moyen des séries de
         ! Fourier, ainsi que ses dérivées premières et secondes

         IMPLICIT NONE

         ! x: point auquel on calcule f
         REAL(kind(0.d0)), dimension(2), intent(in) :: x
         ! Nombre de coefficients de Fourier
         INTEGER, dimension(2), intent(in) :: nFC
         ! Coefficients de Fourier
         REAL(kind(0.d0)), DIMENSION(0:nFC(1)-1,0:nFC(2)-1), intent(in) :: An, Bn, Cn, Dn
         ! Matrice inverse des vecteurs de périodicité at(1:2,i)
         REAL(kind(0.d0)), dimension(1:2,1:2), intent(in) :: inv_at
         ! Valeur de la fonction en x 
         REAL(kind(0.d0)), intent(out) :: f
         ! Derivée première
         REAL(kind(0.d0)), dimension(1:2), intent(out), optional :: df
         ! Derivée seconde
         REAL(kind(0.d0)), dimension(1:2,1:2), intent(out), optional :: d2f

         ! Fréquence de la fonction
         REAL(kind(0.d0)), dimension(1:2) :: U
         REAL(kind(0.d0)), parameter :: pi=3.14159265358979d0
         REAL(kind(0.d0)) :: c1, c2, s1, s2
         REAL(kind(0.d0)), dimension(1:2) :: dfdu
         REAL(kind(0.d0)), dimension(1:2,1:2) :: d2fdu2
         INTEGER :: k, l
         LOGICAL :: calcul_df, calcul_d2f
         
         calcul_df =  Present(df) 
         calcul_d2f = Present(d2f)

         U(:) = 2.d0*pi*MatMul(inv_at(:,:),x(:))

         f = 0.d0
         dfdU(:) = 0.d0
         d2fdU2(:,:) = 0.d0
         DO k = 0, nFC(1)-1
            c1 = cos(dble(k)*U(1))   ; s1 = sin(dble(k)*U(1))
            DO l = 0, nFC(2)-1
                c2 = cos(dble(l)*U(2))   ; s2 = sin(dble(l)*U(2))
                f = f + An(k,l)*c1*c2 + Bn(k,l)*c1*s2 &
                      + Cn(k,l)*s1*c2 + Dn(k,l)*s1*s2
                IF (calcul_df) THEN
                        dfdU(1) = dfdU(1) + dble(k)*( - An(k,l)*s1*c2 - Bn(k,l)*s1*s2 &
                              + Cn(k,l)*c1*c2 + Dn(k,l)*c1*s2 )
                        dfdU(2) = dfdU(2) + dble(l)*( - An(k,l)*c1*s2 + Bn(k,l)*c1*c2 &
                              - Cn(k,l)*s1*s2 + Dn(k,l)*s1*c2 )
                END IF
                IF (calcul_d2f) THEN
                        d2fdU2(1,1) = d2fdU2(1,1) + dble(k*k)*( - An(k,l)*c1*c2 - Bn(k,l)*c1*s2 &
                              - Cn(k,l)*s1*c2 - Dn(k,l)*s1*s2 )
                        d2fdU2(1,2) = d2fdU2(1,2) + dble(k*l)*( An(k,l)*s1*s2 - Bn(k,l)*s1*c2 &
                              - Cn(k,l)*c1*s2 + Dn(k,l)*c1*c2 )
                        d2fdU2(2,2) = d2fdU2(2,2) + dble(l*l)*( - An(k,l)*c1*c2 - Bn(k,l)*c1*s2 &
                              - Cn(k,l)*s1*c2 - Dn(k,l)*s1*s2 )
                END IF
            END DO
         END DO

         IF (calcul_df) THEN
                 df(:) = 2.d0*pi*MatMul( Transpose(inv_at(:,:)), dfdU(:) )
         END IF
         IF (calcul_d2f) THEN
                 d2fdU2(2,1) = d2fdU2(1,2)
                 d2f(:,:) = (2.d0*pi)**2*MatMul( Transpose( MatMul( inv_at(:,:),inv_at(:,:) ) ), d2fDU2(:,:) )
         END IF
        
        END SUBROUTINE TF2D

        SUBROUTINE TF2Dcmplx(x, Fn, nFC, inv_at, f, df, d2f)
         ! Calcul de la fonction f en un point x au moyen des séries de
         ! Fourier, ainsi que ses dérivées premières et secondes
         ! en utilisant les coefficients complexes

         IMPLICIT NONE

         ! x: point auquel on calcule f
         REAL(kind(0.d0)), dimension(2), intent(in) :: x
         ! Nombre de coefficients de Fourier
         INTEGER, dimension(2), intent(in) :: nFC
         ! Coefficients de Fourier complexes
         COMPLEX(kind(0.d0)), dimension( -nFC(1)+1:nFC(1)-1, -nFC(2)+1:nFC(2)-1 ), intent(in) :: Fn
         ! Matrice inverse des vecteurs de périodicité at(1:2,i)
         REAL(kind(0.d0)), dimension(1:2,1:2), intent(in) :: inv_at
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
         REAL(kind(0.d0)), dimension(1:2) :: dfdu
         REAL(kind(0.d0)), dimension(1:2,1:2) :: d2fdu2
         INTEGER :: k, l
         LOGICAL :: calcul_df, calcul_d2f
         
         calcul_df =  Present(df) 
         calcul_d2f = Present(d2f)

         U(:) = 2.d0*pi*MatMul(inv_at(:,:),x(:))

         f = 0.d0
         dfdU(:) = 0.d0
         d2fdU2(:,:) = 0.d0
         DO k = -nFC(1)+1, nFC(1)-1
            e1 = Cmplx( cos(dble(k)*U(1)), sin(dble(k)*U(1)) )
            DO l = -nFC(2)+1, nFC(2)-1
                e2 = Cmplx( cos(dble(l)*U(2)), sin(dble(l)*U(2)) )
                f = f + Dble( Fn(k,l)*e1*e2 )
                IF (calcul_df) THEN
                        !dfdU(1) = dfdU(1) + Cmplx(0.d0, dble(k))*Fn(k,l)*e1*e2
                        !dfdU(2) = dfdU(2) + Cmplx(0.d0, dble(l))*Fn(k,l)*e1*e2
                        dfdU(1) = dfdU(1) + Dble(k)*aImag( Fn(k,l)*e1*e2 )
                        dfdU(2) = dfdU(2) + Dble(l)*aImag( Fn(k,l)*e1*e2 )
                END IF
                IF (calcul_d2f) THEN
                        d2fdU2(1,1) = d2fdU2(1,1) - Dble(k*k)*Dble( Fn(k,l)*e1*e2 )
                        d2fdU2(1,2) = d2fdU2(1,2) + Dble(k*l)*Dble( Fn(k,l)*e1*e2 )
                        d2fdU2(2,2) = d2fdU2(2,2) + Dble(l*l)*Dble( Fn(k,l)*e1*e2 )
                END IF
            END DO
         END DO

         IF (calcul_df) THEN
                 df(:) = 2.d0*pi*MatMul( Transpose(inv_at(:,:)), dfdU(:) )
         END IF
         IF (calcul_d2f) THEN
                 d2fdU2(2,1) = d2fdU2(1,2)
                 d2f(:,:) = (2.d0*pi)**2*MatMul( Transpose( MatMul( inv_at(:,:),inv_at(:,:) ) ), d2fDU2(:,:) )
         END IF
        
        END SUBROUTINE TF2Dcmplx

        SUBROUTINE Write_coefFourier2D(out, An, Bn, Cn, Dn, nFC, at)
         ! Ecriture des coefficients de la série de Fourier dans le fichier
         ! connecté à l'unité out

         IMPLICIT NONE

         ! Unité du fichier de sortie
         INTEGER, intent(in) :: out
         ! Nombre de coefficients de Fourier
         INTEGER, dimension(2), intent(in) :: nFC
         ! Coefficients de Fourier
         REAL(kind(0.d0)), DIMENSION(0:nFC(1)-1,0:nFC(2)-1), intent(in) :: An, Bn, Cn, Dn
         ! Matrice inverse des vecteurs de périodicité at(1:2,i)
         REAL(kind(0.d0)), dimension(1:2,1:2), intent(in) :: at

         INTEGER :: k, l

         WRITE(out,'(a)')        '# Nombre de coefficients de Fourier à calculer dans chaque direction: '
         WRITE(out,'(2(i0,1x))') nFC(1:2)
         WRITE(out,'(a)')       '# Vecteurs de périodicité de la fonction a interpoler: '
         WRITE(out,'(2g20.12)')  at(1:2,1)
         WRITE(out,'(2g20.12)')  at(1:2,2)
         WRITE(out,*)
         WRITE(out,'(a)') "#k, l, An(k,l), Bn(k,l), Cn(k,l), Dn(k,l)"
         DO k=0, nFC(1)-1 
            DO l=0, nFC(2)-1
               WRITE(out,'(2i5,4g24.12)') k, l, An(k,l), Bn(k,l), Cn(k,l), Dn(k,l)
            END DO 
            WRITE(out,*)
         END DO

        END SUBROUTINE Write_coefFourier2D

        SUBROUTINE Read_coefFourier2D(inp, An, Bn, Cn, Dn, nFC, at)
         ! Lecture des coefficients de la série de Fourier dans le fichier
         ! connecté à l'unité inp

         IMPLICIT NONE

         ! Unité du fichier d'entrée
         INTEGER, intent(in) :: inp
         ! Nombre de coefficients de Fourier
         INTEGER, dimension(2), intent(out) :: nFC
         ! Coefficients de Fourier
         REAL(kind(0.d0)), DIMENSION(:,:), allocatable, intent(out) :: An, Bn, Cn, Dn
         ! Matrice inverse des vecteurs de périodicité at(1:2,i)
         REAL(kind(0.d0)), dimension(1:2,1:2), intent(out) :: at

         INTEGER :: k, l, k1, l1

         ! Nombre de coefficients de Fourier à calculer dans chaque direction
         CALL Comment(inp)
         READ(inp,*) nFC(1:2)

         ! Allocation des tableaux
         IF (Allocated(An)) Deallocate(An) ; Allocate(An(0:nFC(1)-1,0:nFC(2)-1)) 
         IF (Allocated(Bn)) Deallocate(Bn) ; Allocate(Bn(0:nFC(1)-1,0:nFC(2)-1)) 
         IF (Allocated(Cn)) Deallocate(Cn) ; Allocate(Cn(0:nFC(1)-1,0:nFC(2)-1)) 
         IF (Allocated(Dn)) Deallocate(Dn) ; Allocate(Dn(0:nFC(1)-1,0:nFC(2)-1)) 

         ! Vecteurs de périodicité de la fonction a interpoler
         CALL Comment(inp); READ(inp,*)  at(1:2,1)
         CALL Comment(inp); READ(inp,*)  at(1:2,2)

         DO k=0, nFC(1)-1 
            DO l=0, nFC(2)-1
               CALL Comment(inp)
               READ(inp,*) k1, l1, An(k,l), Bn(k,l), Cn(k,l), Dn(k,l)
               IF ( (k.NE.k1).OR.(l.NE.l1) ) THEN
                       WRITE(0,'(2(a,i0))') 'read k = ', k1, ' instead of ', k
                       WRITE(0,'(2(a,i0))') 'read l = ', l1, ' instead of ', l
                       STOP '< Read_coefFourier2D >'
               END IF
            END DO 
         END DO

        END SUBROUTINE Read_coefFourier2D

      END MODULE Fourier2DModule
