PROGRAM SymmetrizeEnergy

  ! Apply symmetry x -> -x and z -> -z to file containing stacking fault
  ! energies for prismatic plane
  IMPLICIT NONE

  ! c/a ratio of the hcp structure
  REAL(kind(0.d0)), parameter :: coa=1.601d0

  ! Number of discretization steps in each direction
  INTEGER, parameter :: nx=10, nz=10

  REAL(kind(0.d0)), dimension(1:2,0:nx,0:nz) :: uFault
  REAL(kind(0.d0)), dimension(0:nx,0:nz) :: energy
  REAL(kind(0.d0)), dimension(1:6,0:nx,0:nz) :: sigma

  INTEGER :: out, inp        ! Output and input units
  INTEGER, external :: iArgc
  CHARACTER(len=100) :: input_file
  REAL(kind(0.d0)), dimension(1:11) :: tampon
  INTEGER :: ix, iz, io

  out=6 ! Output is written on screen

  ! Read name of the input file and open it
  IF (iArgc().LE.0) THEN
          WRITE(6,'(a)') 'Name of the input file'
          READ(5,*) input_file
  ELSE
          CALL getArg(1,input_file)
  END IF
  IF (input_file.EQ.'-') THEN
          inp=5 ! Input file read from keyboard
  ELSE
          inp=50
          OPEN(file=input_file,unit=inp,action='read',status='old')
  END IF

  ! Read input energy file
  CALL Comment(inp)
  READ(inp,*,iostat=io) tampon(1:11)

  DO WHILE (io.EQ.0)
      ix=nInt(tampon(1)*dble(nx))
      iz=nInt(tampon(2)*dble(nz))
      uFault(:,ix,iz)=tampon(3:4)     ! Fault vector
      Energy(ix,iz)=tampon(5)         ! Energy
      sigma(:,ix,iz)=tampon(6:11)     ! Stress tensor

      ! Apply symmetry x -> -x
      uFault(:,nx-ix,iz) = (/dble(nx-ix)/dble(nx), dble(iz)/dble(nz)*coa /)
      Energy(nx-ix,iz) = Energy(ix,iz)
      sigma(1:3,nx-ix,iz) =  sigma(1:3,ix,iz)
      sigma(4,nx-ix,iz)   =  sigma(4,ix,iz)
      sigma(5,nx-ix,iz)   = -sigma(5,ix,iz)
      sigma(6,nx-ix,iz)   = -sigma(6,ix,iz)

      ! Apply symmetry z -> -z
      uFault(:,ix,nz-iz) = (/dble(ix)/dble(nx), dble(nz-iz)/dble(nz)*coa /)
      Energy(ix,nz-iz) = Energy(ix,iz)
      sigma(1:3,ix,nz-iz) =  sigma(1:3,ix,iz)
      sigma(4,ix,nz-iz)   = -sigma(4,ix,iz)
      sigma(5,ix,nz-iz)   = -sigma(5,ix,iz)
      sigma(6,ix,nz-iz)   =  sigma(6,ix,iz)

      ! Apply symmetry x -> -x and z -> -z
      uFault(:,nx-ix,nz-iz) = (/dble(nx-ix)/dble(nx), dble(nz-iz)/dble(nz)*coa /)
      Energy(nx-ix,nz-iz) = Energy(ix,iz)
      sigma(1:3,nx-ix,nz-iz) =  sigma(1:3,ix,iz)
      sigma(4,nx-ix,nz-iz)   = -sigma(4,ix,iz)
      sigma(5,nx-ix,nz-iz)   =  sigma(5,ix,iz)
      sigma(6,nx-ix,nz-iz)   = -sigma(6,ix,iz)

      CALL Comment(inp)
      READ(inp,*,iostat=io) tampon(1:11)

  END DO

  IF (inp.NE.5) CLOSE(inp)


  ! Write new energy file
  WRITE(out,'(a)') '# 1:  u1'
  WRITE(out,'(a)') '# 2:  u2'
  WRITE(out,'(a)') '# 3:  ux'
  WRITE(out,'(a)') '# 4:  uy'
  WRITE(out,'(a)') '# 5: Efre'
  WRITE(out,'(a)') '# 6: stress s11'
  WRITE(out,'(a)') '# 7:        s22'
  WRITE(out,'(a)') '# 8:        s33'
  WRITE(out,'(a)') '# 9:        s23'
  WRITE(out,'(a)') '# 10:       s13'
  WRITE(out,'(a)') '# 11:       s12'

  DO ix=0, nx
      DO iz=0, nz
         WRITE(out,'(2(f0.2,1x),2(f0.4,1x),f0.4,1x,6(f0.2,1x))') &
                dble(ix)/dble(nx), dble(iz)/dble(nz), uFault(1:2,ix,iz), &
                Energy(ix,iz), sigma(1:6,ix,iz)
      END DO
      WRITE(out,*)

  END DO

CONTAINS

  SUBROUTINE Comment(unit)

    ! INPUT: unit to read from
    
    ! comment cares about comments and blank lines
    !   and avoids reading of them,
    !   comment lines start with a #

    implicit none

    character(LEN=50) :: phrase
    logical ::  com 
    integer :: unit, i

    com = .true.
    do while ( com .eqv. .true.)
       read(unit,50,END=49) phrase
       if (phrase .NE. ' ') then
          i=index(phrase, '#')
          if (i .ne. 0) then
             if (i .ne. 1) then
                phrase= phrase(:i-1)
                if (phrase .ne. ' ') com= .false.
             endif
          else
             com = .false.
          endif
       endif
    end do

  49 backspace unit
     
    return

  50 format (50a)

  END SUBROUTINE Comment
END PROGRAM SymmetrizeEnergy
