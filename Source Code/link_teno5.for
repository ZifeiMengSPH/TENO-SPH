      subroutine link_teno5
*################################################################################
*-------------------------------------------------------------------------------*
** TENO-SPH CODE
** Name of program: TENO-SPH
** Purpose: to simulate 2D Taylor-Green vortex flow by TENO-SPH
** Input file: input.for
** Output file: output_Binary.for or output_ASCII.for
** Programer: Zifei Meng and Pingping Wang
** Last revision: May 05 2024
** An introduction to TENO-SPH is available in the following references:
** [1] Meng, Z. F., Sun, P. N., Xu, Y., Wang, P. P., & Zhang, A. M. (2024). 
** High-order Eulerian SPH scheme through W/TENO reconstruction based on 
** primitive variables for simulating incompressible flows. CMAME, 427, 117065.
** [2] Meng, Z. F., Zhang, A. M., Wang, P. P., Ming, F. R., & Khoo, B. C. (2022). 
** A targeted essentially non-oscillatory (TENO) SPH method and its applications 
** in hydrodynamics. Ocean Engineering, 243, 110100.
** [3] Wang, P. P., Zhang, A. M., Meng, Z. F., Ming, F. R., & Fang, X. L. (2021). 
** A new type of WENO scheme in SPH for compressible flows with discontinuities. CMAME, 381, 113770.   
** Please cite these references if you use this code in a paper.
*-------------------------------------------------------------------------------*
*################################################################################
      use global
      implicit none
      include 'param.inc'

      integer  i,j,ii
      integer  sumiac, maxiac, noiac, miniac, maxp,minp    
      integer xcell,ycell,zcell,xxcell,yycell,zzcell
      integer ngridx,ngridy,ngridz
      real*8 dgeomx,dgeomy,dgeomz,mhsml
      integer dnxgcell,dnygcell,dnzgcell,dpxgcell,dpygcell,dpzgcell,
     &minxcell,minycell,minzcell,maxxcell,maxycell,maxzcell 
      real*8 tw,tdwdx,tdwdy,tdwdz,dx,dy,dz,r,dr,hvcc,dvx,dvy,dvz
      real*8 hx,hy,hz,hex,hey,hez,rhoij,he
      real*8 eij_x,eij_y
      real*8 p_avg,Lee_temp
      real*8 pv1,pv2,pv
      real*8 u_star, vx_star, vy_star, p_star, temp_1, temp_x, temp_y,ul,ur,Rie_beta,crl
      real*8 psi_temp,deltasph_para1
      real*8 eta,vxy_dwdxy,vyx_dwdxy

      real*8 dvoldt_temp,vx0_star,vy0_star,rho_star,dmdt_temp,dmvxdt_temp,dmvydt_temp,dmedt_temp,e_star,dvofdt_temp,vof_star

      !TENO
      real*8 pl_TENO,pr_TENO,ul_TENO,ur_TENO,rhol_TENO,rhor_TENO,rhol_star,rhor_star
      real*8 pk_1,rhok_1,uk_1,pk_2,rhok_2,uk_2,pk2,rhok2,uk2,pk3,rhok3,uk3,beta_p0,beta_p1,beta_p2,beta_rho0,
     &beta_rho1,beta_rho2,beta_u0,beta_u1,beta_u2,alpha_p0,alpha_p1,alpha_p2,alpha_rho0,alpha_rho1,alpha_rho2,alpha_u0
     &,alpha_u1,alpha_u2,gamma_p0,gamma_p1,gamma_p2,gamma_rho0,gamma_rho1,gamma_rho2,gamma_u0,gamma_u1,
     &gamma_u2,tao_5,wu0,wu1,wu2,wp0,wp1,wp2,wrho0,wrho1,wrho2,wc0,wc1,wc2,CT
      integer id_temp
      real*8 xp_temp,yp_temp,r_temp
      real*8 B11_temp,B12_temp,B21_temp,B22_temp,det,det1,det2
      real*8 ul_temp,ur_temp
      integer sig
      real*8 crl_E,Vrl,SL,SR,SM,ML,MR,EL,ER,p_xing

      real*8 vxk_1,vxk_2,vxk2,vxk3,beta_vx0,beta_vx1,beta_vx2,alpha_vx0,alpha_vx1,alpha_vx2
      real*8 gamma_vx0,gamma_vx1,gamma_vx2,vxl_TENO,vxr_TENO,vxl_rela,vxr_rela

      real*8 ck_1,ck_2,ck2,ck3,beta_c0,beta_c1,beta_c2,alpha_c0,alpha_c1,alpha_c2
      real*8 gamma_c0,gamma_c1,gamma_c2,cl_TENO,cr_TENO,factor

      real*8 vluk_1,vluk_2,vluk2,vluk3,beta_vlu0,beta_vlu1,beta_vlu2,alpha_vlu0,alpha_vlu1,alpha_vlu2
      real*8 gamma_vlu0,gamma_vlu1,gamma_vlu2,vlul_TENO,vlur_TENO
      real*8 dx_temp,dy_temp,p_temp1,p_temp2,rho_temp1,rho_temp2,vx_temp1,vx_temp2,vy_temp1,vy_temp2
      real*8 vx_temp,vy_temp,x_temp,y_temp,q,theta

      q=12.
      CT=1.0e-5
      theta=1.0e-15

!$omp parallel
!$omp do private(i)       
      do i=1,ntotal+nvirt
          niac(i) = 0
          niac_onetype(i)=0
          avdudt(i) = 0.
          indudt(i)=0.      
          indvxdt(i) = 0.
          avdvxdt(i) = 0.
          indvydt(i) = 0.
          avdvydt(i) = 0.
          drhodt(i)=0.
          vcc(i)=0.
          txx(i)=0.
          tzz(i)=0.        
          p_cedian(i)=0.0
          rrr3_x(i)=0.0
          rrr3_y(i)=0.0
          r_he(i)=0.0
          div_r(i)=0.0
          gradx_sebiao(i)=0.0
          grady_sebiao(i)=0.0

          Matrix_B11_b(i)=Matrix_B11(i)
          Matrix_B12_b(i)=Matrix_B12(i)
          Matrix_B21_b(i)=Matrix_B21(i)
          Matrix_B22_b(i)=Matrix_B22(i)

          Matrix_B11(i)=0.0
          Matrix_B12(i)=0.0
          Matrix_B21(i)=0.0
          Matrix_B22(i)=0.0

          grad_c_x(i)=0.0
          grad_c_y(i)=0.0
          mw_total(i)=0.0
          filter_fenzi(i)=0.0
          filter_fenmu(i)=0.0
          w_total(i)=0.0
          unevenness_x(i)=0.0
          unevenness_y(i)=0.0
          unevenness(i)=0.0

          grad_px_b(i)=grad_px(i)
          grad_py_b(i)=grad_py(i)
          grad_rhox_b(i)=grad_rhox(i)
          grad_rhoy_b(i)=grad_rhoy(i)
          grad_vxx_b(i)=grad_vxx(i)
          grad_vxy_b(i)=grad_vxy(i)
          grad_vyx_b(i)=grad_vyx(i)
          grad_vyy_b(i)=grad_vyy(i)

          grad_px(i)=0.0
          grad_py(i)=0.0
          grad_rhox(i)=0.0
          grad_rhoy(i)=0.0
          grad_vxx(i)=0.0
          grad_vxy(i)=0.0
          grad_vyx(i)=0.0
          grad_vyy(i)=0.0

      enddo
!$omp end do
!$omp end parallel


      call check_limits

      call init_grid(ngridx,ngridy,dgeomx,dgeomy)

      do i=1,ntotal+nvirt
          call grid_geom(i,xp(i),yp(i),ngridx,ngridy,
     &    dgeomx,dgeomy,xxcell,yycell)  
          xgcell(i)=xxcell
          ygcell(i)=yycell

          celldata(i) = grid(xxcell,yycell) 
          grid(xxcell,yycell) = i 
      enddo

!!!***************************************dummy boundary (In this case, it is not considered)
!      do i=1,nair-1  
!          p(i)=0.
!          p_cedian(i)=0.
!          minxcell = 1
!          maxxcell = 1
!          minycell = 1
!          maxycell = 1
!           
!	    dnxgcell = xgcell(i) - 1         
!	    dnygcell = ygcell(i) - 1
!	    dpxgcell = xgcell(i) + 1         
!	    dpygcell = ygcell(i) + 1
!          
!          minxcell = max(dnxgcell,1)	    
!	    minycell = max(dnygcell,1)
!          maxxcell = min(dpxgcell,ngridx)	    
!	    maxycell = min(dpygcell,ngridy)
!
!!     Search grid:
!
!          !do zcell=minzcell,maxzcell  
!              do ycell=minycell,maxycell
!                  do xcell=minxcell,maxxcell
!                      j = grid(xcell,ycell)
!
!2                     if(j.gt.i) then 
!                          if(itype(j).eq.3 .or. itype(j).eq.1) then
!                              dx= xp(i) - xp(j)                
!                              dy= yp(i) - yp(j)
!                              dz= zp(i) - zp(j)
!                              dr=dx*dx+dy*dy+dz*dz
!                              r = sqrt(dr)   
!                              mhsml = (h(i)+h(j))/2.
!                              if (r.lt.scale_k*mhsml) then  
!                                  call kernel(r,dx,dy,mhsml,tw,tdwdx,tdwdy)   
!                                  
!                                  p(i)=p(i)+p(j)*tw+rho(j)*((-9.81-accel_y(i))*dy+(-accel_x(i))*dx)*tw 
!                                  vcc(i)=vcc(i)+tw
!                                  p_cedian(i)=p_cedian(i)+p(j)*tw 	            
!                                  tzz(i)=tzz(i)+tw
!                              endif !r<kh 
!
!                      endif
!                      j = celldata(j) 
!                      goto 2  
!                  endif !j>i 
!
!              enddo !xcell
!          enddo
!
!      !enddo  !zcell
!      if(vcc(i).gt.1.e-8) then
!          p(i)=p(i)/vcc(i) 
!      else
!          p(i)=0.
!      endif 
!
!      if(tzz(i).gt.1.e-8) then
!          p_cedian(i)=p_cedian(i)/tzz(i)
!      else
!          p_cedian(i)=0.
!      endif 
!      if(p_cedian(i).lt.pbackground) p_cedian(i)=pbackground
!
!      if(p(i).lt.pbackground) p(i)=pbackground
!      rho(i)=((p(i)-pbackground)/B+1)**(1./7.)*1000.
!      c(i) = cs0*((rho(i)/1000.)**3.)
!      mass(i)=rho(i)*detd*detd
!      enddo 	

      do i=1,ntotal+nvirt
          volume(i)=mass(i)/rho(i) 
      enddo      

!!!************************************************************************************************
!!!!!!!!!!!!!!!!!!call OpenMP parallel to accelerate calculation
      call omp_set_num_threads(num_threads)
!$omp parallel
!$omp do private(i,j,dnxgcell,dnygcell,dpxgcell,dpygcell,
!$omp+  minxcell,minycell,maxxcell,maxycell,xcell,ycell,dx,dy,dvx,dvy,dr,r,mhsml,tw,tdwdx,tdwdy,hvcc,
!$omp+  hx,hy,hex,hey,rhoij,he,p_avg,
!$omp+  Lee_temp,pv1,pv2,pv,
!$omp+  eij_x,eij_y,ul,ur,crl,u_star, vx_star, vy_star,
!$omp+  Rie_beta, p_star, temp_1, temp_x, temp_y, 
!$omp+  psi_temp, deltasph_para1,eta,vxy_dwdxy,vyx_dwdxy,

!$omp+  dvoldt_temp,vx0_star,vy0_star,rho_star,dmdt_temp,dmvxdt_temp,dmvydt_temp,dmedt_temp,e_star,dvofdt_temp,vof_star,
!$omp+  id_temp,xp_temp,yp_temp,r_temp,pl_TENO,pr_TENO,ul_TENO,ur_TENO,
!$omp+  pk_1,rhok_1,uk_1, pk_2,rhok_2,uk_2, pk2,rhok2,uk2, pk3,rhok3,uk3,
!$omp+  beta_p0,beta_p1,beta_p2,beta_rho0,beta_rho1,beta_rho2,beta_u0,beta_u1,beta_u2,
!$omp+  alpha_p0,alpha_p1,alpha_p2,alpha_rho0,alpha_rho1,alpha_rho2,alpha_u0,alpha_u1,alpha_u2,
!$omp+  gamma_p0,gamma_p1,gamma_p2,gamma_rho0,gamma_rho1,gamma_rho2,gamma_u0,gamma_u1,gamma_u2,
!$omp+ rhol_TENO,rhor_TENO,rhol_star,rhor_star,sig,crl_E,Vrl,SL,SR,SM,ML,MR,EL,ER,p_xing,
!$omp+ ck_1,ck_2,ck2,ck3,beta_c0,beta_c1,beta_c2,alpha_c0,alpha_c1,alpha_c2,
!$omp+ gamma_c0,gamma_c1,gamma_c2,cl_TENO,cr_TENO,tao_5,wu0,wu1,wu2,wp0,wp1,wp2,wrho0,wrho1,wrho2,wc0,wc1,wc2,
!$omp+ dx_temp,dy_temp,p_temp1,p_temp2,rho_temp1,rho_temp2,vx_temp1,vx_temp2,vy_temp1,vy_temp2,
!$omp+ vx_temp,vy_temp)


      do i=1,ntotal+nvirt

      minxcell = 1
      maxxcell = 1
      minycell = 1
      maxycell = 1

      dnxgcell = xgcell(i) - 1         
      dnygcell = ygcell(i) - 1
      dpxgcell = xgcell(i) + 1         
      dpygcell = ygcell(i) + 1

      minxcell = max(dnxgcell,1)	    
      minycell = max(dnygcell,1)
      maxxcell = min(dpxgcell,ngridx)	    
      maxycell = min(dpygcell,ngridy)
      do ycell=minycell,maxycell
          do xcell=minxcell,maxxcell
              j = grid(xcell,ycell)

 11           if (j.gt.i) then    

              dx= xp(i) - xp(j)
              dy= yp(i) - yp(j)
              dz= zp(i) - zp(j)
              dvx=vx(i) - vx(j)
              dvy=vy(i) - vy(j)
              dvz=vz(i) - vz(j)
              dr=dx*dx+dy*dy+dz*dz
              r = sqrt(dr)
              mhsml = (h(i)+h(j))/2.

              if (r<scale_k*mhsml) then 

!$omp atomic
              niac(i)=niac(i)+1
!$omp atomic
              niac(j)=niac(j)+1

              call kernel(r,dx,dy,mhsml,tw,tdwdx,tdwdy)   

              !!!**************Riemann-SPH based on TENO reconstruction

              eij_x=-dx/r
              eij_y=-dy/r

              ul=vx(i)*eij_x+vy(i)*eij_y
              ur=vx(j)*eij_x+vy(j)*eij_y 

              !!!!search a particle closest to the missing stencil point
              xp_temp=2*xp(i)-xp(j)
              yp_temp=2*yp(i)-yp(j)
              sig=-1
              call closet(i,j,sig,xp_temp,yp_temp,id_temp,ngridx,ngridy)
              !!!!obtain the values of the missing point

              if(id_temp/=0)then
                  dx_temp= xp_temp-xp(id_temp)
                  dy_temp= yp_temp-yp(id_temp)

                  pk_1=p(id_temp)+( grad_px_b(id_temp)*dx_temp + grad_py_b(id_temp)*dy_temp )
                  vx_temp= vx(id_temp)+( grad_vxx_b(id_temp)*dx_temp + grad_vxy_b(id_temp)*dy_temp )
                  vy_temp= vy(id_temp)+( grad_vyx_b(id_temp)*dx_temp + grad_vyy_b(id_temp)*dy_temp )
                  uk_1=( vx_temp*(-dx) + vy_temp*(-dy) )/r 
              else !!!if it cannot find the closet particle, then the order is reduced. Just in case.
                  pk_1=p(i)
                  uk_1=ul 
              endif    
              !!!!End 


              xp_temp=3*xp(i)-2*xp(j)
              yp_temp=3*yp(i)-2*yp(j)
              sig=-1
              call closet(i,j,sig,xp_temp,yp_temp,id_temp,ngridx,ngridy)
              if(id_temp/=0)then
                  pk_2=p(id_temp)
                  uk_2=( vx(id_temp)*(-dx) + vy(id_temp)*(-dy) )/r
              else
                  pk_2=p(i)
                  uk_2=ul 
              endif

              if(id_temp/=0)then
                  dx_temp= xp_temp-xp(id_temp)
                  dy_temp= yp_temp-yp(id_temp)

                  pk_2=p(id_temp)+( grad_px_b(id_temp)*dx_temp + grad_py_b(id_temp)*dy_temp )
                  vx_temp= vx(id_temp)+( grad_vxx_b(id_temp)*dx_temp + grad_vxy_b(id_temp)*dy_temp )
                  vy_temp= vy(id_temp)+( grad_vyx_b(id_temp)*dx_temp + grad_vyy_b(id_temp)*dy_temp )
                  uk_2=( vx_temp*(-dx) + vy_temp*(-dy) )/r 
              endif   




              xp_temp=2*xp(j)-xp(i)
              yp_temp=2*yp(j)-yp(i)
              sig=1
              call closet(i,j,sig,xp_temp,yp_temp,id_temp,ngridx,ngridy)
              if(id_temp/=0)then
                  pk2=p(id_temp)
                  uk2=( vx(id_temp)*(-dx) + vy(id_temp)*(-dy) )/r
              else
                  pk2=p(j)
                  uk2=ur 
              endif

              if(id_temp/=0)then
                  dx_temp= xp_temp-xp(id_temp)
                  dy_temp= yp_temp-yp(id_temp)

                  pk2=p(id_temp)+( grad_px_b(id_temp)*dx_temp + grad_py_b(id_temp)*dy_temp )
                  vx_temp= vx(id_temp)+( grad_vxx_b(id_temp)*dx_temp + grad_vxy_b(id_temp)*dy_temp )
                  vy_temp= vy(id_temp)+( grad_vyx_b(id_temp)*dx_temp + grad_vyy_b(id_temp)*dy_temp )
                  uk2=( vx_temp*(-dx) + vy_temp*(-dy) )/r 
              endif    


              !******to implement TENO reconstruction for left states
              !***********calculate smoothness indicators
              beta_p2=13/12.*(pk_2-2.*pk_1+p(i))**2+1/4.*(pk_2-4.*pk_1+3.*p(i))**2.
              beta_p0=13/12.*(pk_1-2.*p(i)+p(j))**2+1/4.*(pk_1-p(j))**2.
              beta_p1=13/12.*(p(i)-2.*p(j)+pk2)**2+1/4.*(3.*p(i)-4.*p(j)+pk2)**2.

              beta_u2=13/12.*(uk_2-2.*uk_1+ul)**2+1/4.*(uk_2-4.*uk_1+3.*ul)**2.
              beta_u0=13/12.*(uk_1-2.*ul+ur)**2+1/4.*(uk_1-ur)**2.
              beta_u1=13/12.*(ul-2.*ur+uk2)**2+1/4.*(3.*ul-4.*ur+uk2)**2.

              !***********calculate Scale separation
              tao_5=abs(beta_p1-beta_p2)
              alpha_p0=(1.0+tao_5/(beta_p0+theta))**q 
              alpha_p1=(1.0+tao_5/(beta_p1+theta))**q 
              alpha_p2=(1.0+tao_5/(beta_p2+theta))**q 

              !***********calculate ENO-like stencil selection
              gamma_p0=alpha_p0/(alpha_p0+alpha_p1+alpha_p2)
              gamma_p1=alpha_p1/(alpha_p0+alpha_p1+alpha_p2)
              gamma_p2=alpha_p2/(alpha_p0+alpha_p1+alpha_p2)

              if(gamma_p0<CT)then
                  gamma_p0=0.
              else
                  gamma_p0=1.0
              endif

              if(gamma_p1<CT)then
                  gamma_p1=0.
              else
                  gamma_p1=1.0
              endif

              if(gamma_p2<CT)then
                  gamma_p2=0.
              else
                  gamma_p2=1.0
              endif


              !***********calculate non-linear weights
              wp0=10./16.*gamma_p0/(10./16.*gamma_p0+5./16.*gamma_p1+1./16.*gamma_p2)
              wp1=5./16.*gamma_p1/(10./16.*gamma_p0+5./16.*gamma_p1+1./16.*gamma_p2)
              wp2=1./16.*gamma_p2/(10./16.*gamma_p0+5./16.*gamma_p1+1./16.*gamma_p2)

              tao_5=abs(beta_u1-beta_u2)
              alpha_u0=(1.0+tao_5/(beta_u0+theta))**q
              alpha_u1=(1.0+tao_5/(beta_u1+theta))**q
              alpha_u2=(1.0+tao_5/(beta_u2+theta))**q

              gamma_u0=alpha_u0/(alpha_u0+alpha_u1+alpha_u2)
              gamma_u1=alpha_u1/(alpha_u0+alpha_u1+alpha_u2)
              gamma_u2=alpha_u2/(alpha_u0+alpha_u1+alpha_u2)

              if(gamma_u0<CT)then
                  gamma_u0=0.
              else
                  gamma_u0=1.0
              endif

              if(gamma_u1<CT)then
                  gamma_u1=0.
              else
                  gamma_u1=1.0
              endif

              if(gamma_u2<CT)then
                  gamma_u2=0.
              else
                  gamma_u2=1.0
              endif

              wu0=10./16.*gamma_u0/(10./16.*gamma_u0+5./16.*gamma_u1+1./16.*gamma_u2)
              wu1=5./16.*gamma_u1/(10./16.*gamma_u0+5./16.*gamma_u1+1./16.*gamma_u2)
              wu2=1./16.*gamma_u2/(10./16.*gamma_u0+5./16.*gamma_u1+1./16.*gamma_u2)


              !******calculate the reconstructed left states, such as pl_TENO,ul_TENO
              pl_TENO=1.0/8.0*( wp0*(-pk_1+6.0*p(i)+3.0*p(j))+wp1*(3.0*p(i)+6.0*p(j)-pk2)+
     &        wp2*(3.0*pk_2-10.0*pk_1+15.0*p(i) ))

              ul_TENO=1.0/8.0*( wu0*(-uk_1+6.0*ul+3.0*ur)+wu1*(3.0*ul+6.0*ur-uk2)+
     &        wu2*(3.0*uk_2-10.0*uk_1+15.0*ul) )


              xp_temp=3*xp(j)-2*xp(i)
              yp_temp=3*yp(j)-2*yp(i)
              sig=1
              call closet(i,j,sig,xp_temp,yp_temp,id_temp,ngridx,ngridy)
              if(id_temp/=0)then
                  pk3=p(id_temp)
                  uk3=( vx(id_temp)*(-dx) + vy(id_temp)*(-dy) )/r
              else
                  pk3=p(j)
                  uk3=ur 
              endif

              if(id_temp/=0)then
                  dx_temp= xp_temp-xp(id_temp)
                  dy_temp= yp_temp-yp(id_temp)

                  pk3=p(id_temp)+( grad_px_b(id_temp)*dx_temp + grad_py_b(id_temp)*dy_temp )
                  vx_temp= vx(id_temp)+( grad_vxx_b(id_temp)*dx_temp + grad_vxy_b(id_temp)*dy_temp )
                  vy_temp= vy(id_temp)+( grad_vyx_b(id_temp)*dx_temp + grad_vyy_b(id_temp)*dy_temp )
                  uk3=( vx_temp*(-dx) + vy_temp*(-dy) )/r 
              endif    

              !******to implement TENO reconstruction for right states
              !***********calculate smoothness indicators
              beta_p2=13/12.*(pk_1-2.*p(i)+p(j))**2+1/4.*(pk_1-4.*p(i)+3.*p(j))**2.
              beta_p0=13/12.*(p(i)-2.*p(j)+pk2)**2+1/4.*(pk2-p(i))**2.
              beta_p1=13/12.*(p(j)-2.*pk2+pk3)**2+1/4.*(3.*p(j)-4.*pk2+pk3)**2.

              beta_u2=13/12.*(uk_1-2.*ul+ur)**2+1/4.*(uk_1-4.*ul+3.*ur)**2.
              beta_u0=13/12.*(ul-2.*ur+uk2)**2+1/4.*(uk2-ul)**2.
              beta_u1=13/12.*(ur-2.*uk2+uk3)**2+1/4.*(3.*ur-4.*uk2+uk3)**2.

              !***********calculate Scale separation
              tao_5=abs(beta_p1-beta_p2)
              alpha_p0=(1.0+tao_5/(beta_p0+theta))**q
              alpha_p1=(1.0+tao_5/(beta_p1+theta))**q
              alpha_p2=(1.0+tao_5/(beta_p2+theta))**q

              !***********calculate ENO-like stencil selection
              gamma_p0=alpha_p0/(alpha_p0+alpha_p1+alpha_p2)
              gamma_p1=alpha_p1/(alpha_p0+alpha_p1+alpha_p2)
              gamma_p2=alpha_p2/(alpha_p0+alpha_p1+alpha_p2)

              if(gamma_p0<CT)then
                  gamma_p0=0.
              else
                  gamma_p0=1.0
              endif

              if(gamma_p1<CT)then
                  gamma_p1=0.
              else
                  gamma_p1=1.0
              endif

              if(gamma_p2<CT)then
                  gamma_p2=0.
              else
                  gamma_p2=1.0
              endif

              !***********calculate non-liner weights
              wp0=10./16.*gamma_p0/(10./16.*gamma_p0+1./16.*gamma_p1+5./16.*gamma_p2)
              wp1=1./16.*gamma_p1/(10./16.*gamma_p0+1./16.*gamma_p1+5./16.*gamma_p2)
              wp2=5./16.*gamma_p2/(10./16.*gamma_p0+1./16.*gamma_p1+5./16.*gamma_p2)


              tao_5=abs(beta_u1-beta_u2)
              alpha_u0=(1.0+tao_5/(beta_u0+theta))**q
              alpha_u1=(1.0+tao_5/(beta_u1+theta))**q
              alpha_u2=(1.0+tao_5/(beta_u2+theta))**q

              gamma_u0=alpha_u0/(alpha_u0+alpha_u1+alpha_u2)
              gamma_u1=alpha_u1/(alpha_u0+alpha_u1+alpha_u2)
              gamma_u2=alpha_u2/(alpha_u0+alpha_u1+alpha_u2)

              if(gamma_u0<CT)then
                  gamma_u0=0.
              else
                  gamma_u0=1.0
              endif

              if(gamma_u1<CT)then
                  gamma_u1=0.
              else
                  gamma_u1=1.0
              endif

              if(gamma_u2<CT)then
                  gamma_u2=0.
              else
                  gamma_u2=1.0
              endif

              wu0=10./16.*gamma_u0/(10./16.*gamma_u0+1./16.*gamma_u1+5./16.*gamma_u2)
              wu1=1./16.*gamma_u1/(10./16.*gamma_u0+1./16.*gamma_u1+5./16.*gamma_u2)
              wu2=5./16.*gamma_u2/(10./16.*gamma_u0+1./16.*gamma_u1+5./16.*gamma_u2)


              !******calculate the reconstructed right states, such as pr_TENO,ur_TENO
              pr_TENO=1.0/8.0*(wp0*(3.0*p(i)+6.0*p(j)-pk2)+wp1*(15.0*p(j)-10.0*pk2+3.0*pk3)+
     &        wp2*(-pk_1+6.0*p(i)+3.0*p(j)))

              ur_TENO=1.0/8.0*(wu0*(3.0*ul+6.0*ur-uk2)+wu1*(15.0*ur-10.0*uk2+3.0*uk3)+
     &        wu2*(-uk_1+6.0*ul+3.0*ur))

              !********************************Roe approximate Riemann solver
              crl=(c(i)*rho(i)**1.5+c(j)*rho(j)**1.5)/(sqrt(rho(i))+sqrt(rho(j)))
              u_star=0.5*(ur_TENO+ul_TENO+(pl_TENO-pr_TENO)/crl)
              vx_star=u_star*(-dx)/r+0.5*(vx(i)+vx(j))-0.5*(ur_TENO+ul_TENO)*(-dx)/r
              vy_star=u_star*(-dy)/r+0.5*(vy(i)+vy(j))-0.5*(ur_TENO+ul_TENO)*(-dy)/r
              p_star=0.5*(pl_TENO+pr_TENO)+0.5*crl*(ul_TENO-ur_TENO)

              call real_visc(i,j,tdwdx,tdwdy,dvx,dvy)

              !********************************Riemann-SPH discrete governing equations
              vxy_dwdxy=(vx(i)-vx_star)*tdwdx+(vy(i)-vy_star)*tdwdy
              vyx_dwdxy=(vx_star-vx(j))*tdwdx+(vy_star-vy(j))*tdwdy

!$omp atomic
              drhodt(i) = drhodt(i) + 2.0_8*rho(i)*vxy_dwdxy*mass(j)/rho(j) 
!$omp atomic
              drhodt(j) = drhodt(j) + 2.0_8*rho(j)*vyx_dwdxy*mass(i)/rho(i)                         

              temp_1=p_star/rho(i)/rho(j)
              temp_x=temp_1*tdwdx
              temp_y=temp_1*tdwdy

!$omp atomic
              indvxdt(i)=indvxdt(i)-2.0_8*mass(j)*temp_x
!$omp atomic
              indvxdt(j)=indvxdt(j)+2.0_8*mass(i)*temp_x
!$omp atomic
              indvydt(i)=indvydt(i)-2.0_8*mass(j)*temp_y
!$omp atomic
              indvydt(j)=indvydt(j)+2.0_8*mass(i)*temp_y

              !!!Matrix of Gradient correction 
              Matrix_B11(i) = Matrix_B11(i)+ (-dx)*tdwdx*volume(j)
              Matrix_B12(i) = Matrix_B12(i)+ (-dx)*tdwdy*volume(j)
              Matrix_B21(i) = Matrix_B21(i)+ (-dy)*tdwdx*volume(j)
              Matrix_B22(i) = Matrix_B22(i)+ (-dy)*tdwdy*volume(j)

              Matrix_B11(j) = Matrix_B11(j)+ dx*(-tdwdx)*volume(i)
              Matrix_B12(j) = Matrix_B12(j)+ dx*(-tdwdy)*volume(i)
              Matrix_B21(j) = Matrix_B21(j)+ dy*(-tdwdx)*volume(i)
              Matrix_B22(j) = Matrix_B22(j)+ dy*(-tdwdy)*volume(i)

              p_temp1=(p(j)-p(i))*tdwdx
              p_temp2=(p(j)-p(i))*tdwdy
              rho_temp1=(rho(j)-rho(i))*tdwdx
              rho_temp2=(rho(j)-rho(i))*tdwdy
              vx_temp1=(vx(j)-vx(i))*tdwdx
              vx_temp2=(vx(j)-vx(i))*tdwdy
              vy_temp1=(vy(j)-vy(i))*tdwdx
              vy_temp2=(vy(j)-vy(i))*tdwdy
              !!!!Some gradients or divergences of pressure, density, velocity
              grad_px(i) = grad_px(i)+ p_temp1*volume(j)
              grad_py(i) = grad_py(i)+ p_temp2*volume(j)
              grad_rhox(i) = grad_rhox(i)+ rho_temp1*volume(j)
              grad_rhoy(i) = grad_rhoy(i)+  rho_temp2*volume(j)
              grad_vxx(i)=grad_vxx(i)+vx_temp1*volume(j)
              grad_vxy(i)=grad_vxy(i)+vx_temp2*volume(j)
              grad_vyx(i)=grad_vyx(i)+vy_temp1*volume(j)
              grad_vyy(i)=grad_vyy(i)+vy_temp2*volume(j)

              grad_px(j) = grad_px(j)+ p_temp1*volume(i)
              grad_py(j) = grad_py(j)+ p_temp2*volume(i)
              grad_rhox(j) = grad_rhox(j)+ rho_temp1*volume(i)
              grad_rhoy(j) = grad_rhoy(j)+  rho_temp2*volume(i)
              grad_vxx(j)=grad_vxx(j)+vx_temp1*volume(i)
              grad_vxy(j)=grad_vxy(j)+vx_temp2*volume(i)
              grad_vyx(j)=grad_vyx(j)+vy_temp1*volume(i)
              grad_vyy(j)=grad_vyy(j)+vy_temp2*volume(i)


              endif  !r<kh
              j = celldata(j)
              goto 11  
          endif !j>i 

      enddo
      enddo
      enddo !ntotal   
!$omp end do
!$omp end parallel


!!!!!!!!!!!implement gradient correction procedure
      do i=nair,ntotal+nvirt
          B11_temp=Matrix_B11(i)
          B12_temp=Matrix_B12(i)
          B21_temp=Matrix_B21(i)
          B22_temp=Matrix_B22(i)

          !!!*****************************************for 1D problem
          !det=B11_temp
          !
          !if(abs(det)<1.0e-8)then
          !    write(*,*)'det==0 in Link_CSPH'
          !    pause
          !endif
          !Matrix_B11(i)=1.0/det
          !
          !grad_px(i)= grad_px(i)*Matrix_B11(i)
          !grad_rhox(i)= grad_rhox(i)*Matrix_B11(i)
          !grad_vxx(i)= grad_vxx(i)*Matrix_B11(i)

          !!!*****************************************for 2D problem
          det=B11_temp*B22_temp-B12_temp*B21_temp

          if(abs(det)<1.0e-8)then
              write(*,*)'det==0 in Link_CSPH',i
              pause
          endif

          Matrix_B11(i)= B22_temp/det
          Matrix_B12(i)=-B12_temp/det
          Matrix_B21(i)=-B21_temp/det
          Matrix_B22(i)= B11_temp/det

          x_temp=grad_px(i)
          y_temp=grad_py(i)
          grad_px(i)   = Matrix_B11(i)*x_temp + Matrix_B12(i)*y_temp
          grad_py(i)   = Matrix_B21(i)*x_temp + Matrix_B22(i)*y_temp

          x_temp=grad_rhox(i)
          y_temp=grad_rhoy(i)
          grad_rhox(i) = Matrix_B11(i)*x_temp + Matrix_B12(i)*y_temp
          grad_rhoy(i) = Matrix_B21(i)*x_temp + Matrix_B22(i)*y_temp

          x_temp=grad_vxx(i)
          y_temp=grad_vxy(i)
          grad_vxx(i) = Matrix_B11(i)*x_temp + Matrix_B12(i)*y_temp
          grad_vxy(i) = Matrix_B21(i)*x_temp + Matrix_B22(i)*y_temp

          x_temp=grad_vyx(i)
          y_temp=grad_vyy(i)
          grad_vyx(i) = Matrix_B11(i)*x_temp + Matrix_B12(i)*y_temp
          grad_vyy(i) = Matrix_B21(i)*x_temp + Matrix_B22(i)*y_temp
      enddo 


      return
      end
