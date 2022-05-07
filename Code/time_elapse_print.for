      subroutine time_elapsed(s)
*################################################################################
*-------------------------------------------------------------------------------*
** TENO-SPH CODE
** Name of program: SPH
** Purpose: to simulate 2D Taylor-Green vortex flow by TENO-SPH
** Input file: input.for
** Output file: output_Binary.for or output_ASCII.for
** Programer: Zifei Meng and Pingping Wang
** Last revision: May 02 2022
** The introduction of TENO-SPH can be found in Reference [1]
** [1] Meng, Z. F., Zhang, A. M., Wang, P. P., Ming, F. R., & Khoo, B. C. (2022). 
** A targeted essentially non-oscillatory (TENO) SPH method and its applications 
** in hydrodynamics. Ocean Engineering, 243, 110100.
** If you use this code to write a paper, please cite the above Reference.  
** Copyright 2022 Harbin Engineering University All rights reserved.
*-------------------------------------------------------------------------------*
*################################################################################
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

