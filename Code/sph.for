      program SPH
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

      use global
      implicit none     
      include 'param.inc'

      integer yesorno,i
      real*8 s1, s2
      character(len=100) :: inquire_path,create_path

      allocate(vx_min(maxn))
      allocate(vy_min(maxn))
      allocate(vz_min(maxn))
      allocate(u_min(maxn))
      allocate(rho_min(maxn))
      allocate(niac0(maxn))
      allocate(niac(maxn))
      allocate(xp_min(maxn))
      allocate(yp_min(maxn))
      allocate(grid(maxngx,maxngy))
      allocate(dvxdt(maxn))
      allocate(dvydt(maxn))      
      allocate(dvzdt(maxn))
      allocate(dudt(maxn))
      allocate(xgcell(maxn))
      allocate(ygcell(maxn))
      allocate(zgcell(maxn))
      allocate(celldata(maxn))
      allocate(c(maxn))    
      allocate(vcc(maxn))    
      allocate(vp(maxn))

      call time_print
      call time_elapsed(s1) 
      allocate(current_ts)
      allocate(nstart)
      allocate(mm)
      allocate(mm2)
      allocate(maxtimestep)
      allocate(itimestep)
      allocate(time)
      allocate(time_out)
      allocate(save_step)



      !!!**************************free_surface
      allocate(lmat(2,2,maxn))
      allocate(bmat(2,2,maxn))
      allocate(fs(maxn))
      allocate(vtemp(2,1,maxn)) 
      allocate(lambda(maxn))
      allocate(norx(maxn))
      allocate(nory(maxn))
      allocate(norz(maxn))      
      allocate(tangx(maxn))      
      allocate(tangy(maxn))      
      allocate(tpx(maxn))      
      allocate(tpy(maxn))      
      allocate(fs_2h(maxn))      
      allocate(marker(maxn))      
      allocate(div_r(maxn))      

      !!!**********************************shifting
      allocate(shift_x(maxn))      
      allocate(shift_y(maxn))      
      allocate(shift_z(maxn))      
      allocate(wdetd)      
      allocate(rrr3_x(maxn))      
      allocate(rrr3_y(maxn))      
      allocate(r_he(maxn))    


      !!!**********************************filter
      allocate(filter_fenzi(maxn))      
      allocate(filter_fenmu(maxn))      
      allocate(w_total(maxn)) 
      allocate(mw_total(maxn))  
      allocate(niac_total)      
      allocate(gradx_sebiao(maxn))      
      allocate(grady_sebiao(maxn))

      allocate(part_id(maxn))
      allocate(nongdu(maxn))
      allocate(niac_onetype(maxn))

      !!!**********************************IPST

      allocate(unevenness_x(maxn))
      allocate(unevenness_y(maxn))
      allocate(unevenness(maxn))     

      allocate(para_h2dx)
      allocate(p_cedian(maxn))
      allocate(dismin_fs(maxn)) 
      allocate(closest_fs(maxn))

      allocate(Matrix_B11(maxn))
      allocate(Matrix_B12(maxn))
      allocate(Matrix_B21(maxn))
      allocate(Matrix_B22(maxn))

      allocate(Matrix_B11_b(maxn))
      allocate(Matrix_B12_b(maxn))
      allocate(Matrix_B21_b(maxn))
      allocate(Matrix_B22_b(maxn))

      allocate(grad_c_x(maxn))
      allocate(grad_c_y(maxn))

      allocate(id_012(100))
      allocate(id_016(100))

      !!!**********************************DSPH
      allocate(c1_p(maxn))
      allocate(c2_p(maxn))
      allocate(c1_v(maxn))
      allocate(c2_v(maxn))
      allocate(dvxdx(maxn))
      allocate(dvydy(maxn))
      allocate(dpdx(maxn))
      allocate(dpdy(maxn))

      !!!**********************************UNDEX
      allocate(grad_grad_P_x(maxn))
      allocate(grad_grad_P_y(maxn))
      allocate(grad_grad_P(maxn))

      allocate(grad_px(maxn))
      allocate(grad_py(maxn))
      allocate(grad_ux(maxn))
      allocate(grad_uy(maxn))
      allocate(grad_ux_b(maxn))
      allocate(grad_uy_b(maxn))
      allocate(grad_Ex(maxn))
      allocate(grad_Ey(maxn))
      allocate(grad_rhox(maxn))
      allocate(grad_rhoy(maxn))
      allocate(grad_vxx(maxn))
      allocate(grad_vxy(maxn))
      allocate(grad_vyx(maxn))
      allocate(grad_vyy(maxn))

      allocate(grad_px_b(maxn))
      allocate(grad_py_b(maxn))
      allocate(grad_Ex_b(maxn))
      allocate(grad_Ey_b(maxn))
      allocate(grad_rhox_b(maxn))
      allocate(grad_rhoy_b(maxn))
      allocate(grad_vxx_b(maxn))
      allocate(grad_vxy_b(maxn))
      allocate(grad_vyx_b(maxn))
      allocate(grad_vyy_b(maxn))

      allocate(log_rhob_he(maxn))
      allocate(volume(maxn))
      allocate(vx0(maxn))
      allocate(vy0(maxn))
      allocate(e_total(maxn))
      allocate(dvoldt(maxn))
      allocate(dmdt(maxn))
      allocate(dmvxdt(maxn))
      allocate(dmvydt(maxn))
      allocate(dmedt(maxn))
      allocate(dedt(maxn))
      allocate(delta_vx(maxn))
      allocate(delta_vy(maxn))
      allocate(mass_min(maxn))
      allocate(vof(maxn))
      allocate(vof_min(maxn))
      allocate(dvofdt(maxn))
      allocate(dim)
      allocate(vx_smooth(maxn))
      allocate(vy_smooth(maxn))

      allocate(xp_1(maxn))
      allocate(yp_1(maxn))
      allocate(vx_1(maxn))
      allocate(vy_1(maxn))
      allocate(rho_1(maxn))
      allocate(u_1(maxn))
      allocate(xp_2(maxn))
      allocate(yp_2(maxn))
      allocate(vx_2(maxn))
      allocate(vy_2(maxn))
      allocate(rho_2(maxn))
      allocate(u_2(maxn))


      allocate(dummy_gaodu)
      allocate(radius1)
      allocate(x0_water)
      allocate(x1_water)
      allocate(y0_water)
      allocate(y1_water)
      allocate(y_air)
      allocate(width)
      allocate(sponge_x_start)
      allocate(sponge_z_start)
      allocate(sponge_houdu)
      allocate(y_water_p0)
      allocate(accel_x(maxn))
      allocate(accel_y(maxn))
      allocate(txx(maxn))
      allocate(tzz(maxn))
      allocate(bigUdot)
      allocate(bigVdot)

      !*********energy
      allocate(ek)
      allocate(ep)
      allocate(ec)
      allocate(et)
      allocate(ek0)
      allocate(ep0)
      allocate(ep1)
      allocate(et0)

      current_ts=0         
      nstart=0
      mm=0
      time=0.
      
      write(*,*)'  *****************************************************************'
      write(*,*) 'Check if the needed folders (data, output) exist'
      !!!!!!!!create a folder named 'data' to save the field data that can be opened by ensight software
      inquire_path='./data'
      create_path='data'
      call Createfolder(inquire_path,create_path)    
            
      !!!!!!!!create a folder named 'output' to save the time histories data
      inquire_path='./output'
      create_path='output'
      call Createfolder(inquire_path,create_path)
      
      !!!!The final time of this numerical case.
      time_out=2.0 
      call input

      do i=1,ntotal+nvirt
          dismin_fs(i)=0.0
          closest_fs(i)=0

          Matrix_B11(i)=1.0
          Matrix_B12(i)=0.0
          Matrix_B21(i)=0.0
          Matrix_B22(i)=1.0

          grad_px(i)=0.0
          grad_py(i)=0.0
          grad_rhox(i)=0.0
          grad_rhoy(i)=0.0
          grad_vxx(i)=0.0
          grad_vxy(i)=0.0
          grad_vyx(i)=0.0
          grad_vyy(i)=0.0
          accel_x(i)=0.0
          accel_y(i)=0.0
          dvzdt(i)=0.0
          p_cedian(i)=0.0
      enddo


      !!!**************************free_surface

      do i=1,ntotal+nvirt
          fs(i)=0
          lambda(i)=0.0
          !norx(i)=0.0
          !nory(i)=0.0
          !norz(i)=0.0
          tangx(i)=0.0
          tangy(i)=0.0
          div_r(i)=0.0
      enddo


      !!!**********************************shifting
      do i=1,ntotal+nvirt
          shift_x(i)=0.0
          shift_y(i)=0.0
          shift_z(i)=0.0
          vp(i)=0.0
          mw_total(i)=0.0
      enddo


      !!!**********************************DSPH
      do i=1,ntotal+nvirt
          c1_p(i)=0.0
          c2_p(i)=0.0
          c1_v(i)=0.0
          c2_v(i)=0.0
          dvxdx(i)=0.0
          dvydy(i)=0.0
          dpdx(i)=0.0
          dpdy(i)=0.0
      enddo


      !!!**********************************UNDEX

      call check_limits
      !call output_ASCII
      call output_Binary
      maxtimestep=nint(time_out/dt)      


      call time_integration
      call time_print
      call time_elapsed(s2)      
      write (*,*)'        Elapsed CPU time = ', s2-s1

      end

