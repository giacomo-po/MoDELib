PROGRAM Fit1D

        USE Lecture1DModule
        USE Fourier1DModule

        IMPLICIT NONE

        CHARACTER(len = 100) :: inpFile, outFile, coefFile
        INTEGER :: n, n0, nCol, nxCol
        REAL(kind(0.d0)) :: x0, x1, x, fx, dfdx, d2fdx2, f0, c, sd

        ! Période de la fonction à interpoler
        REAL(kind(0.d0)) :: T

        ! Fonction à interpoler : fVal(:) en xVal(:)
        INTEGER :: nVal
        REAL(kind(0.d0)), dimension(:), allocatable :: fVal, xVal

        ! Coef de Fourier cos et sin : cf(0:nFC-1) et sf(0:nFC-1)
        INTEGER :: nFC
        REAL(kind(0.d0)), dimension(:), allocatable :: cf, sf
        
        WRITE(6,'(a)') 'Number of Fourier coefficients'
        READ(5,*) nFC
        Allocate(cf(0:nFC-1)) ; Allocate(sf(0:nFC-1))
        
        WRITE(6,'(a)') 'Periodicity of the funtion'
        READ(5,*) T

        WRITE(6,'(a)') 'File name where to read the function'
        READ(5,*) inpFile
        WRITE(6,'(a)') 'Column index where to read x position'
        READ(5,*) nxCol

        WRITE(6,'(a)') 'Column index where to read function value'
        READ(5,*) nCol
        OPEN(file=inpFile, unit=50, status='old', action='read')
        nVal = GetNVal1D(50, nxCol, nCol)
        REWIND(50)
        IF (nVal.LE.0) STOP "nVal should be greate than 0" 
        ALLOCATE(xVal(1:nVal)) ; ALLOCATE(fVal(1:nVal))
        xVal(:)=0.d0             ; fVal(:)=0.d0
        CALL ReadVal1D(50, nxCol, nCol, nVal, xVal, fVal)
        CLOSE(50)

        WRITE(6,'(a)') 'Constant to subtract from the function (0. is a good choice)'
        READ(5,*) f0
        fVal(:) = fVal(:) - f0

        WRITE(6,'(a)') 'Factor used to myltiply the function by (1. is a good choice)'
        READ(5,*) c
        fVal(:) = c*fVal(:)

        ! Calcul des coefficients de Fourier
        CALL fitFourier1D(fVal, xVal, nVal, nFC, cf, sf, T)

        ! Fichier de sortie pour les coef de Fourier
        WRITE(6,'(a)') 'Name of the files where to write Fourier coefficients'
        READ(5,*) coefFile
        OPEN(file=coefFile, unit=61, action='write')
        WRITE(61,'(a,i0)')    '# Number of values of the function to be interpolated: ', nVal
        WRITE(61,'(a,i0)')     '# Number of Fourier coefficients: ', nFC
        WRITE(61,'(a,g20.12)') '# Periodicity of the funtion: ', T
        Write(61,'(2a)')       '# File name where to read the function: ', Trim(inpFile)
        WRITE(61,'(a,i0)')     '# Column index where to read x position: ', nxCol
        WRITE(61,'(a,i0)')     '# Column index where to read function value: ', nCol
        WRITE(61,'(a,g20.12)') '# Constant to subtract from the function (0. is a good choice): ', f0
        WRITE(61,'(a,g20.12)') '# Factor used to myltiply the function by (1. is a good choice): ', c
        WRITE(61,'(a)')        "#"
        WRITE(61,'(a)') "#n, c(n), s(n)"
        DO n=0, nFC-1
           WRITE(61,'(i0,2g24.12)') n, cf(n), sf(n)
        END DO
        CLOSE(61)
       
        ! Fichier de sortie pour l'interpolation
        WRITE(6,'(a)') 'Interpolation: origin of the segment'
        READ(5,*) x0
        WRITE(6,'(a)') 'Interpolation: end of the segment'
        READ(5,*) x1
        WRITE(6,'(a)') 'Interpolation: number of points'
        READ(5,*) n0
        WRITE(6,'(a)') 'Name of the file where to write the Fourier interpolation'
        READ(5,*) outFile

        OPEN(file=outFile, unit=62, action='write')
        WRITE(62,'(a,i0)')    '# Number of values of the function to be interpolated: ', nVal
        WRITE(62,'(a,i0)')     '# Number of Fourier coefficients: ', nFC
        WRITE(62,'(a,g20.12)') '# Periodicity of the funtion: ', T
        Write(62,'(2a)')       '# File name where to read the function: ', Trim(inpFile)
        WRITE(62,'(a,i0)')     '# Column index where to read x position: ', nxCol
        WRITE(62,'(a,i0)')     '# Column index where to read function value: ', nCol
        WRITE(62,'(a,g20.12)') '# Constant to subtract from the function (0. is a good choice): ', f0
        WRITE(62,'(a,g20.12)') '# Factor used to myltiply the function by (1. is a good choice): ', c
        WRITE(62,'(a,g20.12)') "# Interpolation: origin of the segment: ", x0
        WRITE(62,'(a,g20.12)') "# Interpolation: end of the segment: ", x1
        WRITE(62,'(a,i0)')     "# Interpolation: number of points: ", n0
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
         DO n=1, nVal
            CALL TF1D(xVal(n), cf, sf, nFC, T, fx, dfdx, d2fdx2)
            sd = sd + ( fx -fVal(n) )**2
         END DO
         sd = sqrt( sd/dble(nVal) )
         WRITE(6,*)
         WRITE(6,'(a,g20.12)') 'Ecart type: ', sd

        Deallocate(fVal) ; Deallocate(cf) ; Deallocate(sf)

END PROGRAM Fit1D
