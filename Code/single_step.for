      subroutine single_step
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
!     Subroutine to calculate the pressure and update the state of governing equations 
      use global
      implicit none
      include 'param.inc'

      integer i
      real*8 det_r
      integer num_tnt
      real*8 mass_total, energy_total, volume_total
      real*8 vp_max
      real*8 pmax_water

      call periodic_BC

      do i=nstart_water,ntotal!+nvirt !calculate the pressure of water
          call water_Tait(p(i),rho(i),c(i))
      enddo
      
      !!!!call TENO-SPH scheme to update the governing equations
      call link_teno5

      do i=1,ntotal!+nvirt
          dvxdt(i) = 0.
          dvydt(i) = 0. 
          dvzdt(i) = 0. 
      enddo


      do i=nair,ntotal!+nvirt
          dvxdt(i) = indvxdt(i) + avdvxdt(i)+grx         
          dvydt(i) = indvydt(i) + avdvydt(i)+gry  
      enddo

      return
      end
