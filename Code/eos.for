      subroutine water_Tait(pii,rhoi,ci)
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
      !--- Tait EoS for nearly incompressible water ---
      use global
      implicit none
      include 'param.inc'

      real*8 rho0,pii,rhoi,ci,gamma   
      integer i_gamma 
      gamma=7.
      rho0= 1000.0
      i_gamma=int(gamma)
      pii = B * ( (rhoi/rho0)**i_gamma - 1.)  
      ci = cs0*((rhoi/rho0)**3)   
      return 
      end

      subroutine p_air(pii,rhoi,ci)
      use global
      implicit none
      include 'param.inc'

      real*8 pii,rhoi,ci,rhoo0,gamma,bair

      rhoo0=1.0
      gamma=1.4
      bair=cs1*cs1*1.0/1.4
      pii =bair * ( (rhoi/rhoo0)**gamma - 1.)
      ci=cs1*(rhoi/rhoo0)**0.2
      return
      end


      subroutine p_gas(pii,rhoi,ei,ci)
      !--- Ideal EoS for gas ---
      implicit none
      real*8 pii,rhoi,ei,gamma,ci
      gamma=3
      pii=(gamma-1)*rhoi*ei
      ci=sqrt((gamma-1) * ei)
      return
      end


