      subroutine real_visc(i,j,tdwdx,tdwdy,dvx,dvy)
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
!     Subroutine to calculate the real viscosity (see Sun, 2015) 
!     Sun, P., Ming, F., & Zhang, A. (2015). 
!     Numerical simulation of interactions between free surface and rigid body using a robust SPH method. Ocean Engineering, 98, 32-49.
      use global
      implicit none
      include 'param.inc'

      real*8  alpha, beta, etq, piv,muv,vr, rr,tdwdx,tdwdy
      real*8  hx,hy, mc, mrho,mhsml,dvx,dvy,dx,dy,Re,rho0
      integer i,j

      etq   = 0.1e0
      rho0=1000.
      alpha = 2.*(real(dim)+2.)
      Re=1000.

      mhsml= (h(i)+h(j))/2.
      vr = 0.
      rr = 0.

      dx=xp(i)-xp(j)

      dy=yp(i)-yp(j)

      vr = dvx*dx+dvy*dy
      rr =dx*dx+dy*dy


      muv = vr/(rr + mhsml*mhsml*etq*etq)
      piv  =  alpha*muv             

      hx=piv*tdwdx
      hy=piv*tdwdy

!$omp atomic
      avdvxdt(i)=avdvxdt(i)+mass(j)/rho(j)/rho(i)*rho0*(1./Re)*hx
!$omp atomic
      avdvydt(i)=avdvydt(i)+mass(j)/rho(j)/rho(i)*rho0*(1./Re)*hy
!$omp atomic
      avdvxdt(j)=avdvxdt(j)-mass(i)/rho(i)/rho(j)*rho0*(1./Re)*hx
!$omp atomic
      avdvydt(j)=avdvydt(j)-mass(i)/rho(i)/rho(j)*rho0*(1./Re)*hy


      return
      end


