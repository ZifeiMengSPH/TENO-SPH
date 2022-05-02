      module global
!----------------------------------------------------------------------
!     A module to define global variables 
      real*8, pointer:: vx(:),vy(:),vz(:)
      real*8, pointer:: xp(:),yp(:),zp(:),detdR(:)
      real*8, pointer::p(:) 
      real*8, pointer::mass(:) 
      real*8, pointer::rho(:)
      real*8, pointer::h(:)
      real*8, pointer::u(:),E(:)
      integer,pointer:: itype(:)
      real*8,pointer::dvxdt(:),dvydt(:),dvzdt(:),dudt(:)
      real*8, pointer::avdvxdt(:),avdvydt(:),avdvzdt(:)
      real*8, pointer::indvxdt(:),indvydt(:),indvzdt(:) 
      real*8, pointer::avdudt(:),indudt(:)
      real*8, pointer::drhodt(:),c(:)
      real*8, pointer::time,momentum_x,momentum_y,time_out         
      integer,pointer::maxtimestep,itimestep,ntotal,ntotal_0,save_step
      integer,pointer::current_ts, nstart,mm,mm2
      real*8,pointer::vx_min(:),vy_min(:),vz_min(:),u_min(:),rho_min(:)
      real*8,pointer::xp_min(:),yp_min(:)
      integer,pointer::grid(:,:),xgcell(:),ygcell(:),zgcell(:)
      integer,pointer::celldata(:),niac(:)
      real*8,pointer::ymin,ymax,xmax,xmin,zmin,zmax  
      real*8,pointer::x_maxgeom,x_mingeom,y_maxgeom,y_mingeom
      real*8,pointer::z_mingeom,z_maxgeom
      real*8,pointer::detd,h_max,cs0,B
      integer,pointer:: nvirt,nstart_water,num_threads
      integer,pointer::niac0(:)
      real*8,pointer::vcc(:)
      real*8,pointer::dt
      real*8,pointer::vp(:)
      real*8,pointer::nei1_x,nei1_y


      !!!***************************APR
      integer,pointer::activity(:)
      integer,pointer::is_split(:)
      integer,pointer::id_region(:)
      integer,pointer::id_set(:)
      real*8,pointer::vx_apr(:)
      real*8,pointer::vy_apr(:)
      real*8,pointer::rho_apr(:)
      real*8,pointer::p_apr(:)
      real*8,pointer::u_apr(:)
      real*8,pointer::c_apr(:)
      real*8,pointer::w_apr(:)

      !!!***************************bgm
      real*8,pointer::xp_0(:)
      real*8,pointer::yp_0(:)
      real*8,pointer::zp_0(:)
      integer,pointer::nvirt_0,nvirt_nobgm

      !!!**************************free_surface
      real*8,pointer::vtemp(:,:,:),lmat(:,:,:),bmat(:,:,:)
      real*8,pointer::lambda(:),norx(:),nory(:),norz(:)
      real*8,pointer::tangx(:),tangy(:),tpx(:),tpy(:)
      real*8,pointer::marker(:)
      integer,pointer::fs(:)
      integer,pointer::fs_2h(:)
      real*8,pointer::div_r(:)

      !!!**********************************shifting
      real*8,pointer::wdetd
      real*8,pointer::shift_x(:),shift_y(:),shift_z(:)
      real*8,pointer::rrr3_x(:),rrr3_y(:),r_he(:)

      !!!**********************************filter
      real*8,pointer::filter_fenzi(:),filter_fenmu(:)
      real*8,pointer::w_total(:), mw_total(:)

      !!!**********************************packing
      real*8,pointer::pack_duxdt(:)
      real*8,pointer::pack_duydt(:)
      real*8,pointer::pack_ux(:)
      real*8,pointer::pack_uy(:)

      !!!**********************************particle pair
      integer,pointer::pair_i(:),pair_j(:)
      real*8,pointer ::w(:),dwdx(:),dwdy(:)

      real*8,pointer::p_background
      integer,pointer::niac_total
      real*8,pointer::gradx_sebiao(:),grady_sebiao(:)

      integer,pointer::part_id(:)

      real*8,pointer::nongdu(:)
      integer,pointer::niac_onetype(:)

      !!!**********************************IPST

      real*8,pointer::unevenness_x(:)
      real*8,pointer::unevenness_y(:)
      real*8,pointer::unevenness(:)

      real*8,pointer::air_y
      real*8,pointer::para_h2dx
      real*8,pointer::p_cedian(:)

      real*8,pointer::dismin_fs(:)
      integer,pointer::closest_fs(:)

      real*8,pointer::Matrix_B11(:)
      real*8,pointer::Matrix_B12(:)
      real*8,pointer::Matrix_B21(:)
      real*8,pointer::Matrix_B22(:)
      real*8,pointer::Matrix_B11_b(:)
      real*8,pointer::Matrix_B12_b(:)
      real*8,pointer::Matrix_B21_b(:)
      real*8,pointer::Matrix_B22_b(:)

      real*8,pointer::grad_c_x(:)
      real*8,pointer::grad_c_y(:)


      integer,pointer::id_012(:)
      integer,pointer::id_016(:)

      real*8,pointer::grx
      real*8,pointer::gry
      integer,pointer::nstart_air
      real*8,pointer::cs1

      !!!**********************************DSPH
      real*8,pointer::c1_p(:)
      real*8,pointer::c2_p(:)
      real*8,pointer::c1_v(:)
      real*8,pointer::c2_v(:)
      real*8,pointer::dvxdx(:)
      real*8,pointer::dvydy(:)
      real*8,pointer::dpdx(:)
      real*8,pointer::dpdy(:)

      !!!**********************************UNDEX
      real*8,pointer::grad_grad_P_x(:)
      real*8,pointer::grad_grad_P_y(:)
      real*8,pointer::grad_grad_P(:)

      real*8,pointer::grad_px(:),grad_Ex(:),grad_ux(:)
      real*8,pointer::grad_py(:),grad_Ey(:),grad_uy(:)
      real*8,pointer::grad_rhox(:)
      real*8,pointer::grad_rhoy(:)
      real*8,pointer::grad_vxx(:)
      real*8,pointer::grad_vxy(:)
      real*8,pointer::grad_vyx(:)
      real*8,pointer::grad_vyy(:)

      real*8,pointer::grad_px_b(:),grad_Ex_b(:),grad_ux_b(:)
      real*8,pointer::grad_py_b(:),grad_Ey_b(:),grad_uy_b(:)
      real*8,pointer::grad_rhox_b(:)
      real*8,pointer::grad_rhoy_b(:)
      real*8,pointer::grad_vxx_b(:)
      real*8,pointer::grad_vxy_b(:)
      real*8,pointer::grad_vyx_b(:)
      real*8,pointer::grad_vyy_b(:)

      real*8,pointer::log_rhob_he(:)
      real*8,pointer::volume(:)
      real*8,pointer::vx0(:),vy0(:),e_total(:)
      real*8,pointer::dvoldt(:),dmdt(:),dmvxdt(:),dmvydt(:),dmedt(:),dedt(:)
      real*8,pointer::delta_vx(:),delta_vy(:)
      real*8,pointer::mass_min(:)
      real*8,pointer::vof(:),vof_min(:),dvofdt(:)

      real*8,pointer::vx_smooth(:),vy_smooth(:)

      integer,pointer::dim

      real*8,pointer::xp_1(:),yp_1(:),vx_1(:),vy_1(:),rho_1(:),u_1(:)
      real*8,pointer::xp_2(:),yp_2(:),vx_2(:),vy_2(:),rho_2(:),u_2(:)

      integer,pointer::record_id(:)
      integer,pointer::record_id2(:)

      !****energy
      real*8,pointer:: ek,ep,ec,et,ek0,ep0,et0,ep1
      !******rigid body
      real*8,pointer:: dummy_gaodu,radius1,x0_water,x1_water,y0_water,y1_water,y_air,width
     &,sponge_x_start,sponge_z_start,sponge_houdu,y_water_p0
      real*8,pointer::bigMass,bigInertiaXX,bigInertiaYY,bigInertiaZZ
      real*8,pointer::Box_XC,Box_YC,Box_ZC,bigU,bigV,bigW,bigOmegaX
      real*8,pointer::bigOmegaY,bigOmegaZ,bigUdot,bigVdot,bigWdot,accel_x(:),accel_y(:),txx(:),tzz(:)
      real*8,pointer::pbackground


      integer,pointer::nbfm,nb_FB,nair

      end module 


