PROGRAM Interpolation1D

        USE LectureGrid1DModule
        USE Fourier1DModule

        IMPLICIT NONE

        CHARACTER(len = 100) :: inpFile, outFile, coefFile
        INTEGER :: n, n0, nCol, nxCol
        REAL(kind(0.d0)) :: x0, x1, x, fx, dfdx, d2fdx2, f0, c, sd

        ! Période de la fonction à interpoler
        REAL(kind(0.d0)) :: T

        ! Fonction à interpoler : f(0:nVal-1)
        INTEGER :: nVal
        REAL(kind(0.d0)), dimension(:), allocatable :: f

        ! Coef de Fourier cos et sin : cf(0:nFC-1) et sf(0:nFC-1)
        INTEGER :: nFC
        REAL(kind(0.d0)), dimension(:), allocatable :: cf, sf
        
        WRITE(6,'(a)') 'Nombre de valeurs de la fonction a interpoler'
        READ(5,*) nVal
        Allocate(f(0:nVal-1))

        WRITE(6,'(a)') 'Nombre de coefficients de Fourier a calculer'
        READ(5,*) nFC
        Allocate(cf(0:nFC-1)) ; Allocate(sf(0:nFC-1))
        
        WRITE(6,'(a)') 'Periode de la fonction a interpoler'
        READ(5,*) T

        Write(6,'(a)') 'Nom du fichier de données où lire la fonction'
        READ(5,*) inpFile
        WRITE(6,'(a)') 'Numéro de colonne où lire la position'
        READ(5,*) nxCol
        WRITE(6,'(a)') 'Numéro de colonne où lire la fonction'
        READ(5,*) nCol
        OPEN(file=inpFile, unit=50, status='old', action='read')
        CALL ReadGridFile1D(50, nxCol, nCol, nVal, T, f)
        CLOSE(50)

        Write(6,'(a)') 'Constante à retrancher de la fonction'
        READ(5,*) f0
        f(:) = f(:) - f0

        WRITE(6,'(a)') 'Facteur multiplicatif'
        READ(5,*) c
        f(:) = c*f(:)

        ! Calcul des coefficients de Fourier
        CALL coefFourier1D(f, nVal, nFC, cf, sf)

        ! Fichier de sortie pour les coef de Fourier
        WRITE(6,'(a)') "Fichier de sortie pour les coef de Fourier"
        READ(5,*) coefFile
        OPEN(file=coefFile, unit=61, action='write')
        WRITE(61,'(a,i0)')     '# Nombre de valeurs de la fonction a interpoler: ', nVal
        WRITE(61,'(a,i0)')     '# Nombre de coefficients de Fourier a calculer: ', nFC
        WRITE(61,'(a,g20.12)') '# Periode de la fonction a interpoler: ', T
        Write(61,'(2a)')       '# Nom du fichier de données où lire la fonction: ', inpFile
        WRITE(61,'(a,i0)')     '# Numéro de colonne où lire la position: ', nxCol
        WRITE(61,'(a,i0)')     '# Numéro de colonne où lire la fonction: ', nCol
        WRITE(61,'(a,g20.12)') '# Constante à retrancher de la fonction: ', f0
        WRITE(61,'(a,g20.12)') '# Facteur multiplicatif: ', c
        WRITE(61,'(a)')        "#"
        WRITE(61,'(a)') "#n, c(n), s(n)"
        DO n=0, nFC-1
           WRITE(61,'(i0,2g24.12)') n, cf(n), sf(n)
        END DO
        CLOSE(61)
       
        ! Fichier de sortie pour l'interpolation
        WRITE(6,'(a)') "Interpolation: origine du segment"
        READ(5,*) x0
        WRITE(6,'(a)') "Interpolation: extrémité du segment"
        READ(5,*) x1
        WRITE(6,'(a)') "Interpolation: nombre de valeurs"
        READ(5,*) n0
        WRITE(6,'(a)') "Fichier de sortie pour l'interpolation"
        READ(5,*) outFile

        OPEN(file=outFile, unit=62, action='write')
        WRITE(62,'(a,i0)')     '# Nombre de valeurs de la fonction a interpoler: ', nVal
        WRITE(62,'(a,i0)')     '# Nombre de coefficients de Fourier a calculer: ', nFC
        WRITE(62,'(a,g20.12)') '# Periode de la fonction a interpoler: ', T
        Write(62,'(2a)')       '# Nom du fichier de données où lire la fonction: ', Trim(inpFile)
        WRITE(62,'(a,i0)')     '# Numéro de colonne où lire la position: ', nxCol
        WRITE(62,'(a,i0)')     '# Numéro de colonne où lire la fonction: ', nCol
        WRITE(62,'(a,g20.12)') '# Constante à retrancher de la fonction: ', f0
        WRITE(62,'(a,g20.12)') '# Facteur multiplicatif: ', c
        WRITE(62,'(a,g20.12)') "# Interpolation: origine du segment: ", x0
        WRITE(62,'(a,g20.12)') "# Interpolation: extrémité du segment: ", x1
        WRITE(62,'(a,i0)')     "# Interpolation: nombre de valeurs: ", n0
        WRITE(62,'(a)')        "#"
        WRITE(62,'(a)') '# 1: x'
        WRITE(62,'(a)') '# 2: f(x)'
        WRITE(62,'(a)') '# 3: df/dx'
        WRITE(62,'(a)') '# 4: d^2 f / d x^2'
        DO n=0, n0
           !x=dble(n)/dble(n0)*T
           x = x0 + dble(n)/dble(n0)*(x1-x0)
           CALL TF1D(x, cf, sf, nFC, T, fx, dfdx, d2fdx2)
           WRITE(62,'(4g24.12)') x, fx, dfdx, d2fdx2
        END DO
        CLOSE(62)

         ! Calcul de l'écart type
         sd = 0.d0
         DO n=0, nVal-1
            x = dble(n)/dble(nVal)*T
            CALL TF1D(x, cf, sf, nFC, T, fx, dfdx, d2fdx2)
            sd = sd + ( fx -f(n) )**2
         END DO
         sd = sqrt( sd/dble(nVal) )
         WRITE(6,*)
         WRITE(6,'(a,g20.12)') 'Ecart type: ', sd

        Deallocate(f) ; Deallocate(cf) ; Deallocate(sf)

END PROGRAM Interpolation1D
