      subroutine init_grid (ngridx,ngridy,
     &dgeomx,dgeomy)
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
!     Subroutine to establish the initial grid for linked list search method (see Liu & Liu 2003) 
!     Liu, G.R., Liu, M.B., 2003. Smoothed Particle Hydrodynamics: A Meshfree Particle Method. World Scientific

      use global
      implicit none
      include 'param.inc'
      integer ngridx,ngridy,ngridz
      real*8 dgeomx,dgeomy,dgeomz,hh2
      integer i,j,k

      ngridx = 1
      ngridy = 1    

      dgeomx = x_maxgeom- x_mingeom
      dgeomy = y_maxgeom- y_mingeom

      hh2=scale_k*h_max
      ngridx=min(int(dgeomx/hh2) + 1,maxngx)      
      ngridy=min(int(dgeomy/hh2) + 1,maxngy)


      do i=1,ngridx
          do j=1,ngridy
              grid(i,j)=0
          enddo
      enddo

      return
      end
