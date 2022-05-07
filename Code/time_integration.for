      subroutine time_integration
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
!     Subroutine to implement the predictor¨Ccorrector time-integration scheme 
      use global
      implicit none     
      include 'param.inc'

      real*8 DD,dt2,dis1,dis2
      integer i,ii,j,k,iiiii,iii,tnt_zongshu
      real*8 xx,xx0,xx1,xx2,kk0,kk1,kk2
      real*8 yy0(6),yy1(6),yy2(6),yy(6)
      real*8 vp_max,shift_xishu
      real*8 x_xiu,y_xiu,factor
      real*8 juli0_tnt2water, juli0_water2tnt  
      real*8  tw,V_j,V_i
      real*8  xleft_dummy, xright_dummy, ydown_dummy
      real max_x, min_x, max_y, min_y

      !*****max_x, min_x, max_y, min_y are boundary lines or mirror lines. 
      max_x=1.0
      min_x=0.0
      max_y=1.0
      min_y=0.0     

      do itimestep = nstart+1, nstart+maxtimestep 

      current_ts=current_ts+1
      if (mod(itimestep,print_step).eq.0) then
          write(*,*)'______________________________________________'
          write(*,*)'  current number of time step =',
     &    itimestep,'     current time=', real(time+dt)
          write(*,*)'______________________________________________'
      endif     

      dt2=dt/2.0 

      time = time + dt2

      call single_step     

!$omp parallel
!$omp do private(i)  
      do i=nair,ntotal!+nvirt          

      xp_min(i)  = xp(i)
      yp_min(i)  = yp(i)
      vx_min(i)  = vx(i)
      vy_min(i)  = vy(i)
      rho_min(i) = rho(i)
      enddo
!$omp end do
!$omp end parallel

!$omp parallel
!$omp do private(i,factor)  
      do i=nair,ntotal!+nvirt
          rho(i) = rho_min(i) + dt2*drhodt(i) 
          xp(i) = xp_min(i) + dt2*vx(i)
          yp(i) = yp_min(i) + dt2*vy(i)
          vx(i) = vx_min(i) + dt2*dvxdt(i)
          vy(i) = vy_min(i) + dt2*dvydt(i)
      enddo
!$omp end do
!$omp end parallel 

!!!!!!!!!!!!!to advoid particle penetration---->periodic boundary 
!$omp parallel
!$omp do private(i)      
      do i=1,ntotal!+nvirt
          if(xp(i)<min_x)then
              xp(i)=xp(i)+max_x
          elseif(xp(i)>max_x)then
              xp(i)=xp(i)-max_x
          elseif(yp(i)>max_y)then
              yp(i)=yp(i)-max_y
          elseif(yp(i)<min_y)then
              yp(i)=yp(i)+max_y
          endif

      enddo
!$omp end do
!$omp end parallel       


!$omp parallel
!$omp do private(i)      
      do i=1,ntotal!+nvirt
          vp(i)=sqrt(vx(i)**2+vy(i)**2)
      enddo
!$omp end do
!$omp end parallel 


      time = time + dt2

      call single_step

!$omp parallel
!$omp do private(i,factor)      
      do i=nair,ntotal!+nvirt
          rho(i) = rho_min(i) + dt*drhodt(i) 
          xp(i) = xp_min(i) + dt*vx(i)
          yp(i) = yp_min(i) + dt*vy(i)
          vx(i) = vx_min(i) + dt*dvxdt(i)
          vy(i) = vy_min(i) + dt*dvydt(i)  
      enddo
!$omp end do
!$omp end parallel 

!!!!!!!!!!!!!to advoid particle penetration---->periodic boundary
!$omp parallel
!$omp do private(i)      
      do i=1,ntotal!+nvirt
          if(xp(i)<min_x)then
              xp(i)=xp(i)+max_x
          elseif(xp(i)>max_x)then
              xp(i)=xp(i)-max_x
          elseif(yp(i)>max_y)then
              yp(i)=yp(i)-max_y
          elseif(yp(i)<min_y)then
              yp(i)=yp(i)+max_y
          endif

      enddo
!$omp end do
!$omp end parallel  


!$omp parallel
!$omp do private(i)      
      do i=1,ntotal!+nvirt
          vp(i)=sqrt(vx(i)**2+vy(i)**2)
      enddo
!$omp end do
!$omp end parallel 

!!!!!!!!!!!!!!!!!!!particle shifting technique
      call PST

!$omp parallel
!$omp do private(i,x_xiu,y_xiu)          
      do i=1,ntotal!+nvirt
          xp(i)=xp(i)+shift_x(i)
          yp(i)=yp(i)+shift_y(i)
      enddo
!$omp end do
!$omp end parallel


!!!**********************************output results
!!!!!!!!set total step of the ensight data
      save_step=nint(time_out/dt) 
      save_step=nint(save_step/200.) 

      if (mod(itimestep,save_step).eq.0) then
          !!output the ensight data in ASCII form
              !call output_ASCII  
                  !!output the ensight data in binary form
                      call output_Binary 
                  endif 

                  !******************************************output time histories

                  ek=0.
                  ep=0.
                  et=0.

                  if (mod(itimestep,100).eq.0) then


                  do i=nstart_water,ntotal!+nvirt 
                      ek=ek+0.5*mass(i)*vp(i)**2.
                      !ep=ep+mass(i)*9.81*yp(i)
                  enddo

                  et=ek+ep

                  open(323,file="./output/ek.dat",position='append')
                  write(323,*)time,ek/et0
                  close(323)

                  open(423,file="./output/ep.dat",position='append')
                  write(423,*)time,ep/et0
                  close(423)

                  open(523,file="./output/et.dat",position='append')
                  write(523,*)time,et/et0
                  close(523)


                  endif



              enddo

              nstart=current_ts

              return

          end
