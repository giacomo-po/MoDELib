MODULE Symmetry2DModule

  SAVE

  ! Maximal number of symmetry operations
  INTEGER, parameter :: max_nSym=100


  TYPE :: sym2Dtype
          ! Rotation matrix 
          REAL(kind(0.d0)), dimension(1:2,1:2) :: rot
          ! Translation vector
          REAL(kind(0.d0)), dimension(1:2) :: u
  END TYPE

  REAL(kind(0.d0)), parameter, private :: zero=1.d-4
  REAL(kind(0.d0)), parameter, private :: zero2=zero*zero

  PUBLIC:: ReadSymmetry2D, ApplySymmetry2D, CheckBasisSymmetry2D
  PRIVATE :: InitSymmetry2D

CONTAINS

  SUBROUTINE ReadSymmetry2D(sym2Dout, nSym, inp)
  
    IMPLICIT NONE
    ! Symmetry operations
    TYPE(sym2Dtype), dimension(1:max_nSym), intent(out) :: sym2Dout
    ! Number of symmetry operations
    INTEGER, intent(out) :: nSym
    ! Input unit
    INTEGER, intent(in) :: inp
  
    ! Local symmetry type to allow definition of an application point when reading
    ! symmetry operations
    TYPE :: sym2DReadType
            ! Rotation matrix 
            REAL(kind(0.d0)), dimension(1:2,1:2) :: rot
            ! Application point
            REAL(kind(0.d0)), dimension(1:2) :: x0
            ! Translation vector
            REAL(kind(0.d0)), dimension(1:2) :: u
    END TYPE
    TYPE(sym2DReadType), dimension(1:max_nSym) :: sym2D
    INTEGER :: ns
  
    NAMELIST /symmetry/ nSym, sym2D
  
    ! Initialize symmetry operations before reading
    nSym = 0
    DO ns=1, max_nSym
       sym2D(ns)%rot(:,:) = 0.d0
       sym2D(ns)%x0(:)    = 0.d0
       sym2D(ns)%u(:)     = 0.d0
    END DO
  
    ! Read namelist
    READ(inp,nml=symmetry)
    IF (nSym.GT.max_nSym) STOP '< ReadSymmetry >: increase max_nSym'
  
    ! Eliminate application point for each symmetry operation by converting
    ! it in a translation vector
    DO ns=1, nSym
       sym2Dout(ns)%rot(:,:) = sym2D(ns)%rot(:,:)
       sym2Dout(ns)%u(:)     = sym2D(ns)%u(:) + sym2D(ns)%x0(:) &
               - MatMul( sym2D(ns)%rot(:,:), sym2D(ns)%x0(:) )
    END DO
    DO ns=nSym+1, max_nSym
       sym2Dout(ns)%rot(:,:) = 0.d0
       sym2Dout(ns)%u(:)     = 0.d0
    END DO
  
  END SUBROUTINE ReadSymmetry2D
  
  SUBROUTINE InitSymmetry2D(sym2D, nSym)
  
    IMPLICIT NONE
    ! Symmetry operation
    TYPE(sym2Dtype), dimension(1:max_nSym), intent(out) :: sym2D
    ! Number of symmetry operations
    INTEGER, intent(out) :: nSym
  
    INTEGER :: ns
  
    nSym=0
    DO ns=1, max_nSym
       sym2D(ns)%rot(:,:)=0.d0
       sym2D(ns)%u(:)=0.d0
    END DO
  
  END SUBROUTINE InitSymmetry2D
  
  SUBROUTINE ApplySymmetry2D(nVal, xVal, fVal, nNew, xNew, fNew, sym2D, nSym)
  
    IMPLICIT NONE
    ! Number of values
    INTEGER, intent(in) :: nVal
    INTEGER, intent(out) :: nNew
    ! Coordinates of all points: xVal(1:2,n)
    REAL(kind(0.d0)), dimension(:,:), intent(in) :: xVal
    REAL(kind(0.d0)), dimension(:,:), allocatable, intent(out) :: xNew
    ! Value of the function in the corresponding point: fVal(n)
    REAL(kind(0.d0)), dimension(:), intent(in) :: fVal
    REAL(kind(0.d0)), dimension(:), allocatable, intent(out) :: fNew
    ! Symmetry operation
    TYPE(sym2Dtype), dimension(1:max_nSym), intent(in) :: sym2D
    ! Number of symmetry operations
    INTEGER, intent(in) :: nSym
  
    INTEGER :: n, ns, m
  
    IF (nSym.LE.0) STOP '< ApplySymmetry >: nSym lower than 1'
  
    nNew = nVal*nSym
    ALLOCATE(xNew(1:2,1:nNew))
    ALLOCATE(fNew(1:nNew))
  
    m = 0
    DO n=1, nVal
       DO ns=1, nSym
          !m = (n-1)*nSym + ns
          m = m + 1
          xNew(1:2,m) = MatMul( sym2D(ns)%rot(1:2,1:2), xVal(1:2,n) ) &
                  + sym2D(ns)%u(1:2)
          fNew(m) = fVal(n)
       END DO
    END DO
  
  END SUBROUTINE ApplySymmetry2D
  
  SUBROUTINE CheckBasisSymmetry2D(at, inv_at, sym2D, nSym, out)
          ! Check that periodicity vectors defining basis
          !   are compatible with symmetry operations
  
          IMPLICIT NONE
          ! Periodicity vectors: at(1:2,i)
          !   and inverse matrix
          REAL(kind(0.d0)), dimension(1:2,1:2), intent(in) :: at, inv_at
          ! Symmetry operation
          TYPE(sym2Dtype), dimension(1:max_nSym), intent(in) :: sym2D
          ! Number of symmetry operations
          INTEGER, intent(in) :: nSym
          ! Output unit
          INTEGER, intent(in) :: out
  
          INTEGER :: i, ns
          REAL(kind(0.d0)), dimension(1:2) :: atNew, s, ds
          LOGICAL :: test
  
          !WRITE(out,'(a,i0)') 'nSym = ', nSym           !! DEBUG
          !WRITE(out,'(a,i0)') 'max_nSym = ', max_nSym           !! DEBUG
          !WRITE(out,'(a,2g20.12)') 'at(1:2,1) = ', at(1:2,1)    !! DEBUG
          !WRITE(out,'(a,2g20.12)') 'at(1:2,2) = ', at(1:2,2)    !! DEBUG
          !WRITE(out,'(a,2g20.12)') 'inv_at(1:2,1) = ', inv_at(1:2,1)    !! DEBUG
          !WRITE(out,'(a,2g20.12)') 'inv_at(1:2,2) = ', inv_at(1:2,2)    !! DEBUG
          !CALL PrintSpaceGroup2D(sym2D, nSym, out)      !! DEBUG
          test=.true.
          DO i=1, 2
             DO ns=1, nSym
  
                ! Apply symmetry operations
                atNew(1:2) = MatMul(sym2D(ns)%rot(1:2,1:2), at(1:2,i) )
  
                ! Reduced units
                s(1:2) = MatMul( inv_at(1:2,1:2), atNew(1:2) )
                ds(1:2) = s(1:2) - aNInt(s(1:2))
  
                IF ( (Abs(ds(1)).GT.zero).OR.(Abs(ds(2)).GT.zero) ) THEN
                        WRITE(out,'(2(a,i0))') 'Periodicity vector ', i, &
                                ' is not compatible with symetry operation ', ns
                        test=.FALSE.
                END IF
  
             END DO
          END DO
  
          IF (.NOT.test) stop
  
  END SUBROUTINE CheckBasisSymmetry2D

  FUNCTION EquivalentSymmetry2D(sym1, sym2, at, inv_at) RESULT(test)
     ! Function returns .true. or .false. if the symmetry operations are
     ! equivalent or not

     IMPLICIT NONE
     TYPE(sym2Dtype), intent(in) :: sym1, sym2
     LOGICAL :: test
     ! Periodicity vectors: at(1:2,i)
     !   and inverse matrix
     REAL(kind(0.d0)), dimension(1:2,1:2), intent(in) :: at, inv_at

     REAL(kind(0.d0)) :: rotNorm2, uNorm2
     REAL(kind(0.d0)), dimension(1:2) :: du, dur

     ! Compare rotation matrices
     rotNorm2 = Sum( ( sym1%rot(:,:) - sym2%rot(:,:) )**2 )
     IF (rotNorm2.GT.zero2) THEN
             test=.FALSE.
             RETURN
     END IF

     ! Compare translation vectors, taking into account periodicity
     du(:) = sym1%u(:) - sym2%u(:)
     dur(:) = MatMul( inv_at(:,:), du(:) )
     dur(:) = aNInt( dur(:) )
     du(:) = du(:) - MatMul( at(:,:), dur(:) )
     uNorm2 = Sum( du(:)**2 )
     IF (uNorm2.GT.zero2) THEN
             test=.FALSE.
             RETURN
     END IF

     test=.TRUE.
     RETURN

  END FUNCTION EquivalentSymmetry2D

  FUNCTION InsideSpaceGroup2D(sym0, sym2D, nSym, at, inv_at) RESULT(test)
     ! Returns .true. if symmetry operation sym0 belongs to sym2D(1:nSym)

     IMPLICIT NONE
     TYPE(sym2Dtype) :: sym0
     TYPE(sym2DType), dimension(1:max_nSym), intent(in) :: sym2D
     INTEGER, intent(in) :: nSym
     ! Periodicity vectors: at(1:2,i)
     !   and inverse matrix
     REAL(kind(0.d0)), dimension(1:2,1:2), intent(in) :: at, inv_at
     LOGICAL :: test

     INTEGER :: ns

     DO ns=1, nSym
        IF ( EquivalentSymmetry2D(sym0, sym2D(ns), at, inv_at) ) THEN
                test=.TRUE.
                RETURN
        END IF
     END DO
     test=.FALSE.
     RETURN

  END FUNCTION InsideSpaceGroup2D

  FUNCTION ComposeSymmetry2D( sym1, sym2, at, inv_at) RESULT(sym3)
     ! Generate symmetry operation sym3 by composing sym1 and sym2

     IMPLICIT NONE
     ! Input symmetry operations
     TYPE(sym2DType), intent(in) :: sym1, sym2
     ! Periodicity vectors: at(1:2,i)
     !   and inverse matrix
     REAL(kind(0.d0)), dimension(1:2,1:2), intent(in) :: at, inv_at
     ! Output symmetry operation
     TYPE(sym2DType) :: sym3

     REAL(kind(0.d0)), dimension(1:2) :: u, ur

     sym3%rot(:,:) = MatMul( sym1%rot(:,:), sym2%rot(:,:) )

     u(:) = sym1%u(:) + MatMul( sym1%rot(:,:), sym2%u(:) )
     ur(:) = MatMul( inv_at(:,:), u(:) )
     ur(:) = aNInt( ur(:) )
     u(:) = u(:) - MatMul( at(:,:), ur(:) )
     sym3%u(:) = u(:)

  END FUNCTION ComposeSymmetry2D

  SUBROUTINE BuildSpaceGroup2D(sym2D0, nSym0, sym2D, nSym, at, inv_at)
     ! Builds space group sym2D (nSym symmetry operations) containing the input
     ! symmetry operations sym2D0 (nSym0 input symmetry operations)

     IMPLICIT NONE
     ! Input symmetry operations
     TYPE(sym2DType), dimension(1:max_nSym), intent(in) :: sym2D0
     INTEGER, intent(in) :: nSym0
     ! Output symmetry operations (space group)
     TYPE(sym2DType), dimension(1:max_nSym), intent(out) :: sym2D
     INTEGER, intent(out) :: nSym
     ! Periodicity vectors: at(1:2,i)
     !   and inverse matrix
     REAL(kind(0.d0)), dimension(1:2,1:2), intent(in) :: at, inv_at

     INTEGER :: ns0, ns1, ns2
     TYPE(sym2DType) :: sym3

     ! Check that symmetry operations are compatible with periodicity vectors
     CALL CheckBasisSymmetry2D(at, inv_at, sym2D0, nSym0, 6)

     ! First symmetry operations in point group = identity
     nSym = 1
     sym2D(1)%rot(:,:) = Reshape( (/ 1.d0, 0.d0, 0.d0, 1.d0 /), (/2,2/) )
     sym2D(1)%u(:) = 0.d0

     ! Add one instance of each input symmetry operation
     DO ns0=1, nSym0
        IF (.NOT.InsideSpaceGroup2D(sym2D0(ns0), sym2D, nSym, at, inv_at)) THEN
                nSym = nSym + 1
                sym2D(nSym) = sym2D0(ns0)
                !WRITE(6,'(a,i0,a)') 'Structure ', ns0, ' added'          ! DEBUG
        !ELSE                                                             ! DEBUG
                !WRITE(6,'(a,i0,a)') 'Structure ', ns0, ' already inside' ! DEBUG
        END IF
     END DO

     ! Generate new symmetry operations by composing them and add them to the
     ! space group if they are not already included
     ns1 = 2
     loop_sym1: DO

        ns2 = 2
        loop_sym2: DO
         
           sym3 =  ComposeSymmetry2D( sym2D(ns1), sym2D(ns2), at, inv_at)  
           !WRITE(6,*)                                                               ! DEBUG
           !WRITE(6,'(2(a,i0))') "Result of symmetry operations ", ns1, " and ", ns2 ! DEBUG
           !CALL PrintSymmetry2D(sym3, 6)                                            ! DEBUG
           IF (.NOT.InsideSpaceGroup2D(sym3, sym2D, nSym, at, inv_at)) THEN
                   nSym = nSym + 1
                   sym2D(nSym) = sym3
                   !WRITE(6,'(a,i0)') 'Structure added as number ', nSym ! DEBUG
           !ELSE                                                         ! DEBUG
                   !WRITE(6,'(a)') 'Structure already inside'            ! DEBUG
           END IF

           ! Increment symmetry operation 
           ns2 = ns2 + 1
           IF (ns2.GT.nSym) EXIT
        END DO loop_sym2

        ! Increment symmetry operation 
        ns1 = ns1 + 1
        IF (ns1.GT.nSym) EXIT

     END DO loop_sym1

  END SUBROUTINE BuildSpaceGroup2D

  SUBROUTINE PrintSpaceGroup2D(sym2D, nSym, out)

          IMPLICIT NONE

          ! Symmetry operation
          TYPE(sym2Dtype), dimension(1:max_nSym), intent(in) :: sym2D
          ! Number of symmetry operations
          INTEGER, intent(in) :: nSym
          ! Output unit
          INTEGER, intent(in) :: out

          INTEGER :: ns

          WRITE(out,*)
          WRITE(out,'(a,i0)') 'Number of symmetry operations, nSym = ', nSym
          DO ns=1, nSym
             WRITE(out,*)
             WRITE(out,'(a,i0)') '  operation ', ns
             CALL PrintSymmetry2D(sym2D(ns), out)
          END DO
          WRITE(out,*)
  
  END SUBROUTINE PrintSpaceGroup2D

  SUBROUTINE PrintSymmetry2D(sym2D, out)

          IMPLICIT NONE

          ! Symmetry operation
          TYPE(sym2Dtype), intent(in) :: sym2D
          ! Output unit
          INTEGER, intent(in) :: out

          WRITE(out,'(a,2(g13.6,1x))') '  rot = ', sym2D%rot(1,1:2)
          WRITE(out,'(a,2(g13.6,1x))') '        ', sym2D%rot(2,1:2)
          WRITE(out,'(a,2(g13.6,1x))') '    u = ', sym2D%u(1:2)
  
  END SUBROUTINE PrintSymmetry2D

END MODULE Symmetry2DModule
