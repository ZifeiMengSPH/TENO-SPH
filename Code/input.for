      subroutine input
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
!     Subroutine to create model and set initial conditions 
      use global
      implicit none     
      include 'param.inc'

      integer i,j,k,m
      real*8 length_x,length_y,length_z,dis0_p,tnt_r2,x0,y0,z0
      integer geshu_x,geshu_y,geshu_z,nn
      real*8 tw,tdwdx,tdwdy,tdwdz,dx,dy,dz,r,dr,mhsml
      real*8 xp_temp,yp_temp
      real*8 s1,s2,sum,p_temp

      allocate(xp(maxn))
      allocate(yp(maxn))
      allocate(zp(maxn))
      allocate(vx(maxn))
      allocate(vy(maxn))
      allocate(vz(maxn))
      allocate(mass(maxn))
      allocate(rho(maxn))
      allocate(p(maxn))
      allocate(itype(maxn))
      allocate(h(maxn))
      allocate(u(maxn)) 
      allocate(ntotal)
      allocate(xmax)
      allocate(xmin)
      allocate(ymax)
      allocate(ymin)
      allocate(zmax)
      allocate(zmin)
      allocate(detd)
      allocate(x_maxgeom)
      allocate(x_mingeom)
      allocate(y_maxgeom)
      allocate(y_mingeom)
      allocate(z_maxgeom)
      allocate(z_mingeom)
      allocate(h_max)
      allocate(nvirt)
      allocate(avdudt(maxn))
      allocate(indudt(maxn))
      allocate(indvxdt(maxn))
      allocate(avdvxdt(maxn))
      allocate(indvydt(maxn))
      allocate(avdvydt(maxn))
      allocate(indvzdt(maxn))
      allocate(avdvzdt(maxn))
      allocate(drhodt(maxn))
      allocate(nstart_water)
      allocate(B)
      allocate(cs0)
      allocate(cs1)
      allocate(nbfm)
      allocate(nair)
      allocate(nb_FB)

      allocate(dt)
      allocate(num_threads)
      allocate(nei1_x)
      allocate(nei1_y)
      allocate(dim)

      allocate(pbackground)
      allocate(grx)
      allocate(gry)

      !the artificial speed of sound
      cs0=10.   
      !The parameter of Tait eos
      B=cs0*cs0*1000./7. 
      cs1 =cs0
      !The particle spacing
      detd=0.0025

      !gravity acceleration 
      grx = 0.0   
      gry = 0.0

      !The lenth of the computational domain in x,y,z direction
      length_x=(0.5+0.5*detd)*2.0
      length_y=(0.5+0.5*detd)*2.0
      !The particle number in x,y,z direction
      geshu_x=nint(length_x/detd)
      geshu_y=nint(length_y/detd)

      !The size of time step
      dt= detd/cs0/4.
      !Dimension
      dim=2 
      !Background pressure
      pbackground=0.0  

      !The initial smoothing length
      para_h2dx=1.2
      h_max=para_h2dx*detd


      write(*,*)'  *****************************************************************'
      write(*,*) 'Number of CPU threads is ? '
      write(*,*)'  *****************************************************************'
      read (*,*) num_threads 

      !!**************parameters of creating model
      i=1   

      !The particle label of periodic boundary or dynamic boundary
      nvirt=0
      !The first particle label of water
      nstart_water=1   

      x0=0.0
      y0=0.0
      z0=0.0
      nstart_water=i
      !The first particle label of air. It is not considered here, because Taylor-Green flow is a single phase flow.
      nair=i 

      !!!!create model and set initial conditions, refer to the reference below
      !Taylor, G.I., Green, A.E., 1937. Mechanism of the production of small eddies from large ones. Proc. R. Soc. A 158 (895), 499¨C521
      do j=1,geshu_y  
          do k=1,geshu_x  

          xp(i)=x0
          x0=x0+detd
          yp(i)=y0
          zp(i)=z0

          if(xp(i)>=0.0.and.xp(i)<=1.0 .and. yp(i)>=0.0 .and. yp(i)<=1.0 )then !
              itype(i)=3
              vy(i)=-cos(2.0*pi*xp(i))*sin(2.0*pi*yp(i))
              vx(i)= sin(2.0*pi*xp(i))*cos(2.0*pi*yp(i))
              vz(i)=0.0 
              p(i)=0.25*1000.*(cos(4.*pi*xp(i))+cos(4.*pi*yp(i))) 
              rho(i)=((p(i)-pbackground)/B+1.0_8)**(1.0_8/7.0_8)*1000.0_8 
              mass(i)=rho(i)*detd**real(dim)
              h(i)=para_h2dx*detd
              c(i)=cs0
              i=i+1         
          endif    
          enddo
          y0=y0+detd
          x0=0.0 
      enddo

      !The total number of particles in the computational domain except the periodic boundary
      ntotal=i-1
      ! call the subroutine of the periodic boundary
      call periodic_BC

      !***Delete the previous output data files before the new simulation
      
      open(323,file="./output/ek.dat")
      close(323,status="delete")

      open(423,file="./output/ep.dat")
      close(423,status="delete")

      open(523,file="./output/et.dat")      
      close(523,status="delete")

      !**********Record initial energy      
      ek0=0.
      ep0=0.
      et0=0.

      do i=1,ntotal!+nvirt
          vp(i)=sqrt(vx(i)**2+vy(i)**2)
      enddo
      do i=nstart_water,ntotal!+nvirt 
          ek0=ek0+0.5*mass(i)*vp(i)**2
          ep0=ep0+mass(i)*abs(gry)*yp(i)
      enddo
      et0=ek0+ep0

      open(123,file="./output/ek0.dat")
      write(123,*)ek0
      close(123)

      open(223,file="./output/ep0.dat")
      write(223,*)ep0
      close(223)

      open(623,file="./output/et0.dat")
      write(623,*)et0
      close(623)

      !!!******************************calculate w(¦¤x) for particle shifting technique  

      r=detd
      dx=detd/sqrt(2.0)
      dy=detd/sqrt(2.0)
      mhsml=para_h2dx*detd
      call kernel(r,dx,dy,mhsml,tw,tdwdx,tdwdy)
      wdetd=tw
      if(abs(wdetd)<1.0e-6)then
          write(*,*)'wdetd is 0'
          pause
      endif



      write(*,*)'  ****************************************************************'
      write(*,*)'      Initial particle configuration generated   '   
      write(*,*)'      Total number of particles   ', ntotal+nvirt    	
      write(*,*)'  ****************************************************************'
      return
      end              

