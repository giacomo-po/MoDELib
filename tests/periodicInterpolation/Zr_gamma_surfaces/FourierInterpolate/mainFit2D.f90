PROGRAM Fit2D

        USE Lecture2DModule
        USE Math
        USE Fourier2DModule

        IMPLICIT NONE

        CHARACTER(len = 100) :: inpFile, outFile, coefFile
        INTEGER :: n, n0, m, m0, nCol, nxCol, nyCol
        REAL(kind(0.d0)), dimension(1:2) :: x, x0, x1, x2
        REAL(kind(0.d0)) :: fx, f0, c, sd

        ! Vecteurs de périodicité de la fonction à interpoler: at(1:2,i)
        !   et matrice inverse
        REAL(kind(0.d0)), dimension(1:2,1:2) :: at, inv_at

        ! Fonction à interpoler : f() en x(,:)
        INTEGER :: nVal
        REAL(kind(0.d0)), dimension(:,:), allocatable :: xVal
        REAL(kind(0.d0)), dimension(:), allocatable :: fVal

        ! Coef de Fourier cos et sin : cf(0:nFC-1) et sf(0:nFC-1)
        INTEGER, dimension(2) :: nFC
        REAL(kind(0.d0)), dimension(:,:), allocatable :: An, Bn, Cn, Dn
        
        WRITE(6,'(a)') 'Number of Fourier coefficients in each direction'
        READ(5,*) nFC(1), nFC(2)
        Allocate(An(0:nFC(1)-1,0:nFC(2)-1)) 
        Allocate(Bn(0:nFC(1)-1,0:nFC(2)-1)) 
        Allocate(Cn(0:nFC(1)-1,0:nFC(2)-1)) 
        Allocate(Dn(0:nFC(1)-1,0:nFC(2)-1)) 
        
        WRITE(6,'(a)') '1st periodicity vector of the funtion'
        READ(5,*) at(1:2,1)
        WRITE(6,'(a)') '2nd periodicity vector of the function'
        READ(5,*) at(1:2,2)
        CALL Mat2Inv(at, inv_at)

        WRITE(6,'(a)') 'File name where to read the function'
        READ(5,*) inpFile
        WRITE(6,'(a)') 'Column index where to read x position (y position should be in the next column)'
        READ(5,*) nxCol
        nyCol=nxCol+1
        WRITE(6,'(a)') 'Column index where to read function value'
        READ(5,*) nCol
        OPEN(file=inpFile, unit=50, status='old', action='read')
        nVal = GetNVal2D(50, nxCol, nyCol, nCol)
        REWIND(50)
        IF (nVal.LE.0) STOP "nVal should be greate than 0" 
        ALLOCATE(xVal(1:2,1:nVal)) ; ALLOCATE(fVal(1:nVal))
        xVal(:,:)=0.d0             ; fVal(:)=0.d0
        CALL ReadVal2D(50, nxCol, nyCol, nCol, nVal, xVal, fVal)
        CLOSE(50)

        WRITE(6,'(a)') 'Constant to subtract from the function (0. is a good choice)'
        READ(5,*) f0
        fVal(:) = fVal(:) - f0

        WRITE(6,'(a)') 'Factor used to myltiply the function by (1. is a good choice)'
        READ(5,*) c
        fVal(:) = c*fVal(:)

        ! Calcul des coefficients de Fourier
         CALL FitFourier2D(fVal, xVal, nVal, nFC, An, Bn, Cn, Dn, inv_at) 

        ! Fichier de sortie pour les coef de Fourier
        WRITE(6,'(a)') 'Name of the files where to write Fourier coefficients'
        READ(5,*) coefFile
        OPEN(file=coefFile, unit=61, action='write')
        CALL Write_coefFourier2D(61, An, Bn, Cn, Dn, nFC, at)
        CLOSE(61)
       
        ! Fichier de sortie pour l'interpolation
        WRITE(6,'(a)') 'Interpolation: origin of the grid'
        READ(5,*) x0
        WRITE(6,'(a)') 'Interpolation: first corner of the grid'
        READ(5,*) x1
        WRITE(6,'(a)') 'Interpolation: second corner of the grid'
        READ(5,*) x2
        WRITE(6,'(a)') 'Interpolation: number of points on the grid in each direction'
        READ(5,*) n0, m0
        WRITE(6,'(a)') 'Name of the file where to write the Fourier interpolation'
        READ(5,*) outFile
        OPEN(file=outFile, unit=62, action='write')
        WRITE(62,'(a,1(i0,1x))') '# Number of values of the function to be interpolated: ', nVal
        WRITE(62,'(a,2(i0,1x))') '# Number of Fourier coefficients in each direction: ', nFC(1:2)
        WRITE(62,'(a,2g20.12)')  '# 1st periodicity vector of the funtion: ', at(1:2,1)
        WRITE(62,'(a,2g20.12)')  '# 2nd periodicity vector of the function: ', at(1:2,2)
        Write(62,'(2a)')         '# File name where to read the function: ', Trim(inpFile)
        WRITE(62,'(a,i0)')       '# Column index where to read x position (y position should be in the next column): ', nxCol
        WRITE(62,'(a,i0)')       '# Column index where to read function value: ', nCol
        Write(62,'(a,g20.12)')   '# Constant to subtract from the function (0. is a good choice): ', f0
        WRITE(62,'(a,g20.12)')   '# Factor used to myltiply the function by (1. is a good choice): ', c
        WRITE(62,'(a,2g20.12)')  '# Interpolation: origin of the grid:', x0(:)
        WRITE(62,'(a,2g20.12)')  '# Interpolation: first corner of the grid:', x1(:)
        WRITE(62,'(a,2g20.12)')  '# Interpolation: second corner of the grid:', x2(:)
        WRITE(62,'(a,2(i0,1x))') "# Interpolation: number of points on the grid in each direction: ", n0, m0
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
         DO n=1, nVal
            CALL TF2d(xVal(:,n), An, Bn, Cn, Dn, nFC, inv_at, fx)
            sd = sd + ( fx -fVal(n) )**2
         END DO
         sd = sqrt( sd/dble(nVal) )
         WRITE(6,*)
         WRITE(6,'(a,g20.12)') 'Standard deviation: ', sd

        Deallocate(fVal, xVal) 
        Deallocate(An, Bn, Cn, Dn)
        

END PROGRAM Fit2D
