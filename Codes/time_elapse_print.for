      subroutine time_elapsed(s)
!----------------------------------------------------------------------
!     Subroutine to calculate the elapsed CPU Time

      use dfport
      implicit none

      integer, parameter :: output = 6
      real(8) :: s

      s = rtc()

      end subroutine time_elapsed



      subroutine time_print

      implicit none
      integer, parameter :: output = 6

      ! . local scalars.
      character ( len =  8 ) :: datstr
      character ( len = 10 ) :: timstr

      ! . Get the current date and time.
      call date_and_time ( datstr, timstr )

      ! . Write out the date and time.
      write ( output, "(/A)"  ) "                  Date = " // 
     &datstr(7:8) // "/" // 
     &datstr(5:6) // "/" // 
     &datstr(1:4)
      write ( output, "(A)"   ) "                  Time = " // 
     &timstr(1:2) // ":" // 
     &timstr(3:4) // ":" // 
     &timstr(5:10)

      write ( output, *)

      end subroutine time_print

