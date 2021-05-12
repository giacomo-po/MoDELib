      MODULE Fourier1DModule

       IMPLICIT NONE
      
       CONTAINS

       SUBROUTINE fitFourier1D(fVal, xVal, nVal, nFC, An, Bn, T)
         ! Ajustement moindre carré d'une série de Fourier à partir de nVal
         ! points qui ne sont pas forcément répartis sur une grille régulière

         ! Nombre de valeurs de f
         INTEGER, intent(in) :: nVal
         ! fVal(x) valeur de la fonction en x (valeurs discrètes)
         REAL(kind(0.d0)), dimension(1:nVal), intent(in) :: fVal, xVal
         ! Coefficients de la série de Fourier
         REAL(kind(0.d0)), intent(out), DIMENSION(0:nFC-1) :: An, Bn
         ! Nombre de coefficients de Fourier
         INTEGER, intent(in) :: nFC
         ! Période de la fonction
         REAL(kind(0.d0)), intent(in) :: T

         REAL(kind(0.d0)), parameter :: pi=3.14159265358979d0
         REAL(kind(0.d0)) :: f0, fmean, u
         INTEGER :: n, nDim, k, kl, info
         COMPLEX(kind(0.d0)), dimension( 1:nVal, 2*nFC-1 ) :: akl
         COMPLEX(kind(0.d0)), dimension( 2*nFC-1, 2*nFC-1 ) :: Amat
         COMPLEX(kind(0.d0)), dimension( 2*nFC-1 ) :: Bvec
         COMPLEX(kind(0.d0)), dimension( -nFC+1:nFC-1 ) :: Fn
         INTEGER, dimension( 2*nFC-1 ) :: iPiv
         
         ! Calcul de la valeur moyenne
         f0 = fVal(1)
         fmean = f0 + Sum( fVal(1:nVal)-f0 )/dble(nVal) 

         ! Nombre de coefficients complexes à trouver
         nDim = 2*nFC-1

        ! Construit le système linéaire à résoudre
         DO n=1, nVal
            u = 2.d0*pi*xVal(n)/T

            kl=0
            DO k=-nFC+1, nFC-1
               kl=kl+1
               akl(n,kl) = Cmplx( cos(dble(k)*u), sin(dble(k)*u) )
            END DO

         END DO

         ! Les coefficients complexes de la série de Fourier sont solutions de
         ! Amat(:,:)*X(:) = Bvec(:)
         Amat(:,:) = MatMul( Transpose( Conjg(akl(:,:)) ), akl(:,:) )
         Bvec(:) = MatMul( Transpose( Conjg(akl(:,:)) ), Cmplx( fVal(:)-fmean, 0.d0) )
         CALL zgesv(nDim, 1, Amat, nDim, iPiv, Bvec, nDim, info)
         IF (info.NE.0) THEN
                 WRITE(0,'(a,i0)') 'zsysv: Info = ', info
                 STOP '< FitFourier >'
         END IF

         ! Bvec(:) contient les coefficients complexes de la série de Fourier
         !   => réécriture sous forme d'un tableau
         kl=0
         DO k=-nFC+1, nFC-1
                kl=kl+1
                Fn(k) = Bvec(kl)
         END DO         

        ! Coefficients réels des séries de Fourier en cos() et sin()
         DO k=0, nFC-1
               An(k) =   Dble( Fn(k) + Fn(-k) )
               Bn(k) = -aImag( Fn(k) - Fn(-k) )
         END DO

         An(0) = 0.5d0*An(0) 
         Bn(0) = 0.d0

         ! On rajoute la valeur moyenne
         An(0) = An(0)+fmean


       END SUBROUTINE fitFourier1D

        SUBROUTINE coefFourier1D(f, nVal, nFC, An, Bn)
         ! Calcul des coefficients d'une série de Fourier

         ! f(x) valeur de la fonction en x (valeurs discrètes)
         REAL(kind(0.d0)), dimension(0:nVal-1), intent(in) :: f
         ! Nombre de valeurs de f
         INTEGER, intent(in) :: nVal
         ! Coefficients de la série de Fourier
         REAL(kind(0.d0)), intent(out), DIMENSION(0:nFC-1) :: An, Bn
         ! Nombre de coefficients de Fourier
         INTEGER, intent(in) :: nFC

         REAL(kind(0.d0)), parameter :: pi=3.14159265358979d0
         REAL(kind(0.d0)) :: omega, f0, fmean
         INTEGER :: x, k

         omega=2.d0*pi/dble(nVal)

         f0 = f(0)
         fmean = f0 + Sum( f(0:nVal-1)-f0 )/dble(nVal) 
         An(0) = fmean
         Bn(0) = 0.d0

         DO k=1, nFC-1
            An(k)=0.d0 ; Bn(k)=0.d0
            DO x=0, nVal-1
               An(k) = An(k) + (f(x)-fmean)*Cos(dble(k*x)*omega)
               Bn(k) = Bn(k) + (f(x)-fmean)*Sin(dble(k*x)*omega)
            END DO
         END DO
         An(1:) = An(1:)/dble(nVal)
         Bn(1:) = Bn(1:)/dble(nVal)


        END SUBROUTINE coefFourier1D
       
        SUBROUTINE TF1D(x, An, Bn, nFC, T, f, df, d2f)
         ! Calcul de la fonction f en un point x et de ses dérivées premières
         ! et secondes au moyen des séries de Fourier

         ! x: point auquel on calcule f
         REAL(kind(0.d0)), intent(in) :: x
         ! Coefficients de Fourier
         REAL(kind(0.d0)), intent(in), DIMENSION(0:nFC-1) :: An, Bn
         ! Nombre de coefficients de Fourier
         INTEGER, intent(in) :: nFC
         ! Période de la fonction f
         REAL(kind(0.d0)), intent(in) :: T
         ! Valeur de la fonction en x 
         REAL(kind(0.d0)), intent(out) :: f, df, d2f

         ! Fréquence de la fonction
         REAL(kind(0.d0)) :: U, c, s
         REAL(kind(0.d0)), parameter :: pi=3.14159265358979d0
         INTEGER :: i
         
         U = 2.d0*pi/T*x

         f = An(0)  ;  df = 0.d0  ;  d2f = 0.d0
         DO i = 1, nFC-1
                c = cos(dble(i)*U)  ; s = sin(dble(i)*U)
                f = f + An(i)*c + Bn(i)*s
                df = df + dble(i)*( - An(i)*s + Bn(i)*c )
                d2f = d2f - dble(i*i)*( An(i)*c + Bn(i)*s )
         END DO
         df = df*2.d0*pi/T
         d2f = d2f*(2.d0*pi/T)**2
        
        END SUBROUTINE TF1D

      END MODULE Fourier1DModule
