PROGRAM MinPath2D

        USE Constraint2DModule
        USE Fourier2DModule
        USE Math
        USE MinModule

        IMPLICIT NONE

        CHARACTER(len=100) :: coefFile, outFile1, outFile2

        REAL(kind(0.d0)), dimension(1:2) :: x0, x1, x, nDir
        REAL(kind(0.d0)) :: z, t, E, dE, d2E
        INTEGER :: nVal, n, outUnit

        ! Lecture des coefficients de Fourier
        WRITE(6,'(a)') 'Fichier où lire les coefficients de Fourier'
        READ(5,*) coefFile
        OPEN(file=coefFile,unit=51,action='read',status='old')
        CALL Read_coefFourier2D(51, An, Bn, Cn, Dn, nFC, at)
        CLOSE(51)

        ! Inversion des vecteurs de périodicité
        CALL Mat2Inv(at, inv_at)

        ! Extrémités entre lesquelles calculer le chemin
        WRITE(6,'(a)') 'Point initial'
        READ(5,*) x0(:)
        WRITE(6,'(a)') 'Point final'
        READ(5,*) x1(:)

        ! Nombre de points intermédiaires
        WRITE(6,'(a)') 'Nombre de valeurs'
        READ(5,*) nVal

        ! Chemin non relaxé
        WRITE(6,'(a)') 'Fichier de sortie pour le chemin non relaxé'
        READ(5,*) outFile1
        IF (Trim(outFile1).EQ.'-') THEN
                outUnit=6
                WRITE(outUnit,*)
        ELSE
                outUnit=61
                OPEN(file=outFile1,unit=outUnit,action='write')
        END IF
        WRITE(outUnit,'(a)')         '# Chemin non relaxé entre les points'
        WRITE(outUnit,'(a,2g20.12)') '#   x0(:) = ', x0(1:2)
        WRITE(outUnit,'(a,2g20.12)') '#   x1(:) = ', x1(1:2)
        WRITE(outUnit,*)
        WRITE(outUnit,'(a)') '#  1: coordonnée de réaction, z'
        WRITE(outUnit,'(a)') '#  2: vecteur position, x(1)'
        WRITE(outUnit,'(a)') '#  3:                   x(2)'
        WRITE(outUnit,'(a)') '#  4: énergie, E'
        WRITE(outUnit,'(a)') '#  5: dérivée première projetée, dE'
        WRITE(outUnit,'(a)') '#  6:         seconde          , d2E'

        ! Direction de la droite
        nDir(:) = x1(:)-x0(:)
        nDir(:) = nDir(:)/sqrt( Sum( nDir(:)**2 ) )

        DO n=0, nVal
          z = dble(n)/dble(nVal)
          x(:) = x0(:) + z*(x1(:)-x0(:))
          CALL Init_fConstrained(x, nDir)
          CALL fConstrained(0.d0, E, dE, d2E)
          WRITE(outUnit,'(6g20.12)') z, x(1:2), E, dE, d2E
        END DO
        IF (outUnit.NE.6) CLOSE(outUnit)

        ! Chemin relaxé
        WRITE(6,'(a)') 'Fichier de sortie pour le chemin relaxé'
        READ(5,*) outFile2
        IF (Trim(outFile2).EQ.'-') THEN
                outUnit=6
                WRITE(outUnit,*)
        ELSE
                outUnit=62
                OPEN(file=outFile2,unit=outUnit,action='write')
        END IF
        WRITE(outUnit,'(a)')         '# Chemin non relaxé entre les points'
        WRITE(outUnit,'(a,2g20.12)') '#   x0(:) = ', x0(1:2)
        WRITE(outUnit,'(a,2g20.12)') '#   x1(:) = ', x1(1:2)
        WRITE(outUnit,*)
        WRITE(outUnit,'(a)') '#  1: coordonnée de réaction, z'
        WRITE(outUnit,'(a)') '#  2: vecteur position, x(1)'
        WRITE(outUnit,'(a)') '#  3:                   x(2)'
        WRITE(outUnit,'(a)') '#  4: énergie, E'
        WRITE(outUnit,'(a)') '#  5: dérivée première projetée, dE'
        WRITE(outUnit,'(a)') '#  6:         seconde          , d2E'

        ! Direction perpendiculaire à la droite
        nDir(1) =  x1(2)-x0(2)
        nDir(2) = -x1(1)+x0(1)
        nDir(:) = nDir(:)/sqrt( Sum( nDir(:)**2 ) )

        DO n=0, nVal
          z = dble(n)/dble(nVal)
          x(:) = x0(:) + z*(x1(:)-x0(:))
          CALL Init_fConstrained(x, nDir)
          t = 0.d0
          CALL Minimize(fConstrained, E, t, 1.d-6, 1.d-2, 1000)
          CALL fConstrained(t, E, dE, d2E)
          WRITE(outUnit,'(6g20.12)') z, x(1:2)+t*nDir(1:2), E, dE, d2E
        END DO
        IF (outUnit.NE.6) CLOSE(outUnit)



END PROGRAM MinPath2D
