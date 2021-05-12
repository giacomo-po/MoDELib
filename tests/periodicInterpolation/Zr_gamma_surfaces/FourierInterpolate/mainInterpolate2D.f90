PROGRAM Interpolation2D

        USE LectureGrid2DModule
        USE Math
        USE Fourier2DModule

        IMPLICIT NONE

        CHARACTER(len = 100) :: inpFile, outFile, coefFile
        INTEGER :: n, n0, m, m0, nCol, nxCol
        REAL(kind(0.d0)), dimension(1:2) :: x, x0, x1, x2
        REAL(kind(0.d0)) :: fx, f0, c, sd

        ! Vecteurs de périodicité de la fonction à interpoler: at(1:2,i)
        !   et matrice inverse
        REAL(kind(0.d0)), dimension(1:2,1:2) :: at, inv_at

        ! Fonction à interpoler : f(0:nVal(1)-1,0:nVal(2)-1)
        INTEGER, dimension(2) :: nVal
        REAL(kind(0.d0)), dimension(:,:), allocatable :: f

        ! Coef de Fourier cos et sin : cf(0:nFC-1) et sf(0:nFC-1)
        INTEGER, dimension(2) :: nFC
        REAL(kind(0.d0)), dimension(:,:), allocatable :: An, Bn, Cn, Dn
        
        WRITE(6,'(a)') 'Nombre de valeurs dans chaque direction de la fonction a interpoler'
        READ(5,*) nVal(1), nVal(2)
        Allocate(f(0:nVal(1)-1,0:nVal(2)-1))

        WRITE(6,'(a)') 'Nombre de coefficients de Fourier à calculer dans chaque direction'
        READ(5,*) nFC(1), nFC(2)
        Allocate(An(0:nFC(1)-1,0:nFC(2)-1)) 
        Allocate(Bn(0:nFC(1)-1,0:nFC(2)-1)) 
        Allocate(Cn(0:nFC(1)-1,0:nFC(2)-1)) 
        Allocate(Dn(0:nFC(1)-1,0:nFC(2)-1)) 
        
        WRITE(6,'(a)') '1er vecteur de périodicité de la fonction a interpoler'
        READ(5,*) at(1:2,1)
        WRITE(6,'(a)') '2nd vecteur de périodicité de la fonction a interpoler'
        READ(5,*) at(1:2,2)
        CALL Mat2Inv(at, inv_at)

        Write(6,'(a)') 'Nom du fichier de données où lire la fonction'
        READ(5,*) inpFile
        WRITE(6,'(a)') 'Numéro de colonne où lire la position'
        READ(5,*) nxCol
        WRITE(6,'(a)') 'Numéro de colonne où lire la fonction'
        READ(5,*) nCol
        OPEN(file=inpFile, unit=50, status='old', action='read')
        CALL ReadGridFile2D(50, nxCol, nCol, nVal, at, inv_at, f)
        CLOSE(50)

        Write(6,'(a)') 'Constante à retrancher de la fonction'
        READ(5,*) f0
        f(:,:) = f(:,:) - f0

        WRITE(6,'(a)') 'Facteur multiplicatif'
        READ(5,*) c
        f(:,:) = c*f(:,:)

        ! Calcul des coefficients de Fourier
        CALL coefFourier2D(f, nVal, nFC, An, Bn, Cn, Dn)

        ! Fichier de sortie pour les coef de Fourier
        WRITE(6,'(a)') "Fichier de sortie pour les coef de Fourier"
        READ(5,*) coefFile
        OPEN(file=coefFile, unit=61, action='write')
        CALL Write_coefFourier2D(61, An, Bn, Cn, Dn, nFC, at)
        CLOSE(61)
       
        ! Fichier de sortie pour l'interpolation
        WRITE(6,'(a)') 'Interpolation: origine de la grille'
        READ(5,*) x0
        WRITE(6,'(a)') 'Interpolation: extrémité du premier côté de la grille'
        READ(5,*) x1
        WRITE(6,'(a)') 'Interpolation: extrémité du second côté de la grille'
        READ(5,*) x2
        WRITE(6,'(a)') "Interpolation: nombre de valeurs dans chaque direction"
        READ(5,*) n0, m0
        WRITE(6,'(a)') "Fichier de sortie pour l'interpolation"
        READ(5,*) outFile
        OPEN(file=outFile, unit=62, action='write')
        WRITE(62,'(a,2(i0,1x))') '# Nombre de valeurs dans chaque direction de la fonction a interpoler: ', nVal(1:2)
        WRITE(62,'(a,2(i0,1x))') '# Nombre de coefficients de Fourier à calculer dans chaque direction: ', nFC(1:2)
        WRITE(62,'(a,2g20.12)')  '# 1er vecteur de périodicité de la fonction a interpoler: ', at(1:2,1)
        WRITE(62,'(a,2g20.12)')  '# 2nd vecteur de périodicité de la fonction a interpoler: ', at(1:2,2)
        Write(62,'(2a)')         '# Nom du fichier de données où lire la fonction: ', Trim(inpFile)
        WRITE(62,'(a,i0)')       '# Numéro de colonne où lire la position: ', nxCol
        WRITE(62,'(a,i0)')       '# Numéro de colonne où lire la fonction: ', nCol
        Write(62,'(a,g20.12)')   '# Constante à retrancher de la fonction: ', f0
        WRITE(62,'(a,g20.12)')   '# Facteur multiplicatif: ', c
        WRITE(62,'(a,2g20.12)')  '# Interpolation: origine de la grille :', x0(:)
        WRITE(62,'(a,2g20.12)')  '# Interpolation: extrémité du premier côté de la grille :', x1(:)
        WRITE(62,'(a,2g20.12)')  '# Interpolation: extrémité du second côté de la grille :', x2(:)
        WRITE(62,'(a,2(i0,1x))') "# Interpolation: nombre de valeurs dans chaque direction : ", n0, m0
        WRITE(62,*)
        WRITE(62,'(a)') '#x(1), x(2), f(x)'
        DO n=0, n0
           DO m=0, m0
              !x(:) = dble(n)/dble(n0)*at(:,1) + dble(m)/dble(m0)*at(:,2)
              x(:) = x0(:) + dble(n)/dble(n0)*(x1(:)-x0(:)) + dble(m)/dble(m0)*(x2(:)-x0(:))
              CALL TF2D(x, An, Bn, Cn, Dn, nFC, inv_at, fx)
              WRITE(62,'(3g24.12)') x(:), fx
           END DO
           WRITE(62,*)
        END DO
        CLOSE(62)

         ! Calcul de l'écart type
         sd = 0.d0
         DO n=0, nVal(1)-1
            DO m=0, nVal(2)-1
               x(:) = dble(n)/dble(nVal(1))*at(:,1) + dble(m)/dble(nVal(2))*at(:,2)
               CALL TF2d(x(:), An, Bn, Cn, Dn, nFC, inv_at, fx)
               sd = sd + ( fx -f(n,m) )**2
            END DO
         END DO
         sd = sqrt( sd/dble(nVal(1)*nVal(2)) )
         WRITE(6,*)
         WRITE(6,'(a,g20.12)') 'Ecart type: ', sd

        Deallocate(f) 
        Deallocate(An, Bn, Cn, Dn)
        

END PROGRAM Interpolation2D
