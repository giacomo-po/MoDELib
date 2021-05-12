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
!!$     write(6,*) phrase
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
