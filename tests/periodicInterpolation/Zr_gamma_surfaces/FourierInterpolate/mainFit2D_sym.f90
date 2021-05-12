PROGRAM Fit2D

        ! Ajustement d'une série de Fourier 2D 
        !   avec des symétires imposées

        USE Lecture2DModule
        USE Math
        USE Fourier2DModule
        USE Fourier2DSymModule
        USE Symmetry2DModule

        IMPLICIT NONE

        CHARACTER(len = 100) :: inpFile, outFile, coefFile, symFile
        INTEGER :: n, n0, m, m0, nCol, nxCol, nyCol
        REAL(kind(0.d0)), dimension(1:2) :: x, x0, x1, x2
        REAL(kind(0.d0)) :: fx, f0, c, sd

        ! Symmetry operations
        TYPE(sym2Dtype), dimension(1:max_nSym) :: sym2D0, sym2D
        ! Total number of symmetry operations
        INTEGER :: nSym0, nSym

        ! Vecteurs de périodicité de la fonction à interpoler: at(1:2,i)
        !   et matrice inverse
        REAL(kind(0.d0)), dimension(1:2,1:2) :: at, inv_at

        ! Fonction à interpoler : f() en x(,:)
        INTEGER :: nVal
        REAL(kind(0.d0)), dimension(:,:), allocatable :: xVal
        REAL(kind(0.d0)), dimension(:), allocatable :: fVal

        ! Coef de Fourier cos et sin pour la série standard
        INTEGER, dimension(2) :: nFC
        REAL(kind(0.d0)), dimension(:,:), allocatable :: An, Bn, Cn, Dn
        ! Coef de Fourier complexes pour la série standard
        COMPLEX(kind(0.d0)), dimension(:,:), allocatable :: Fn
        ! Coef de Fourier complexes pour la série symmétrisée
        INTEGER, dimension(2) :: nFCsym
        COMPLEX(kind(0.d0)), dimension(:,:), allocatable :: FnSym
        
        WRITE(6,'(a)') 'Number of Fourier coefficients in each direction'
        READ(5,*) nFCsym(1), nFCsym(2)
        
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

        ! ==== Symmetry operations =======================
        WRITE(6,'(a)') 'Nom du fichier contenant les opérations de symétrie'
        READ(5,*) symFile
        OPEN(unit=51, file=symFile, status='old', action='read')
        CALL ReadSymmetry2D(sym2D0, nSym0, 51)
        CLOSE(51)
        CALL CheckBasisSymmetry2D(at, inv_at, sym2D0, nSym0, 6)
        WRITE(6,*)
        WRITE(6,'(a)') '========================='
        WRITE(6,'(a)') 'Input symmetry operations'
        CALL PrintSpaceGroup2D(sym2D0, nSym0, 6)
        ! Generate space group
        WRITE(6,*)
        WRITE(6,'(a)') '========================='
        WRITE(6,'(a)') 'Space group'
        CALL BuildSpaceGroup2D(sym2D0, nSym0, sym2D, nSym, at, inv_at)
        CALL PrintSpaceGroup2D(sym2D, nSym, 6)
        CALL CheckBasisSymmetry2D(at, inv_at, sym2D, nSym, 6)
        ! ================================================

        ! Calcule le nombre de coefficients de la série de Fourier symmétrisée 
        ! réécrite sous sa forme régulière
        CALL nFC2D_sym(nFCsym, at, inv_at, sym2D, nSym, nFC)
        write(6,'(a,2(i0,1x))') 'nFCsym(1:2) = ', nFCsym(1:2)
        write(6,'(a,2(i0,1x))') 'nFC(1:2)    = ', nFC(1:2)

        Allocate(An(0:nFC(1)-1,0:nFC(2)-1)) 
        Allocate(Bn(0:nFC(1)-1,0:nFC(2)-1)) 
        Allocate(Cn(0:nFC(1)-1,0:nFC(2)-1)) 
        Allocate(Dn(0:nFC(1)-1,0:nFC(2)-1)) 
        Allocate(Fn( -nFC(1)+1:nFC(1)-1, -nFC(2)+1:nFC(2)-1 ))
        Allocate(FnSym( -nFCsym(1)+1:nFCsym(1)-1, -nFCsym(2)+1:nFCsym(2)-1 ))

        ! Calcul des coefficients de Fourier
        FnSym(:,:) = Cmplx( 0.d0, 0.d0 )
        !    - par inversion du système linéaire
        CALL FitFourier2DSymCmplx(fVal, xVal, nVal, FnSym, nFCsym, inv_at, sym2D, nSym)
        !    - par minimisation en utilisant l'algorithme Fire
        CALL MinimizeCostFunction2D_sym(fVal, xVal, nVal, FnSym, nFCsym, inv_at, sym2D, nSym)

        CALL coefFourier2D_unsymmetrize(FnSym, nFCsym, Fn, nFC, at, inv_at, sym2D, nSym)
        CALL coefFourier2D_cmplx2real(Fn, nFC, An, Bn, Cn, Dn)

        ! Fichier de sortie pour les coef de Fourier
        WRITE(6,'(a)') 'Name of the files where to write Fourier coefficients'
        READ(5,*) coefFile
        OPEN(file=coefFile, unit=61, action='write')
        nFC(1)=10       !! DEBUG
        nFC(2)=10       !! DEBUG
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
              x(:) = x0(:)
              IF (n0.NE.0) x(:) = x(:) + dble(n)/dble(n0)*(x1(:)-x0(:))
              IF (m0.NE.0) x(:) = x(:) + dble(m)/dble(m0)*(x2(:)-x0(:))
              !x(:) = x0(:) + dble(n)/dble(n0)*(x1(:)-x0(:)) + dble(m)/dble(m0)*(x2(:)-x0(:))
              !CALL TF2D(x, An, Bn, Cn, Dn, nFC, inv_at, fx)
              !CALL TF2Dcmplx(x, Fn, nFC, inv_at, fx)
              CALL TF2DSymCmplx(x, FnSym, nFCsym, inv_at, sym2D, nSym, fx)
              WRITE(62,'(3g24.12)') x(:), fx
           END DO
           IF (m0.NE.0) WRITE(62,*)
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
