      subroutine grid_geom (i,xg,yg,ngridx,ngridy,
     &dgeomx,dgeomy,xxcell,yycell)
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
!     Subroutine to calculate the coordinates of the cell for linked list search method (see Liu & Liu 2003) 
!     Liu, G.R., Liu, M.B., 2003. Smoothed Particle Hydrodynamics: A Meshfree Particle Method. World Scientific
      use global
      implicit none
      include 'param.inc'

      real*8 xg,yg,zg
      integer xxcell,yycell,zzcell
      integer ngridx,ngridy,ngridz
      real*8 dgeomx,dgeomy,dgeomz,hh2
      integer i

      xxcell=1
      yycell=1

      hh2=scale_k*h_max
      xxcell=int((xg-x_mingeom)/hh2 + 1.e0)         
      yycell=int((yg-y_mingeom)/hh2 + 1.e0)
      
      return
      end