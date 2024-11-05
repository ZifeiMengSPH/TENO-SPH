      subroutine closet(id_i,id_j,sig_temp,xg,yg,id_WENO,ngridx,ngridy)
!----------------------------------------------------------------------
!     Subroutine to search a particle closest to the missing stencil point (see Meng, 2022) 
!     Meng, Z. F., Zhang, A. M., Wang, P. P., Ming, F. R., & Khoo, B. C. (2022). 
!     A targeted essentially non-oscillatory (TENO) SPH method and its applications in hydrodynamics. Ocean Engineering, 243, 110100.
      use global
      implicit none
      include 'param.inc'   

      real*8 xg,yg,zg
      integer xxcell,yycell,zzcell,xcell,ycell
      integer ngridx,ngridy,ngridz
      real*8 dgeomx,dgeomy,dgeomz
      integer dnxgcell,dnygcell,dnzgcell,dpxgcell,dpygcell,dpzgcell,
     &minxcell,minycell,minzcell,maxxcell,maxycell,maxzcell 
      integer i,j,id_WENO,xgcell_temp,ygcell_temp
      real*8 dis_min,dx,dy,dr
      integer id_i,id_j,sig_temp


      id_WENO=0
      dis_min=100000.0

      xgcell_temp=int((xg-x_mingeom)/(scale_k*h_max) + 1.e0)         
      ygcell_temp=int((yg-y_mingeom)/(scale_k*h_max) + 1.e0)

      minxcell = 1
      maxxcell = 1
      minycell = 1
      maxycell = 1

      dnxgcell = xgcell_temp - 1         
      dnygcell = ygcell_temp - 1
      dpxgcell = xgcell_temp + 1         
      dpygcell = ygcell_temp + 1

      minxcell = max(dnxgcell,1)	    
      minycell = max(dnygcell,1)
      maxxcell = min(dpxgcell,ngridx)	    
      maxycell = min(dpygcell,ngridy)


      do ycell=minycell,maxycell 
          do xcell=minxcell,maxxcell
              j = grid(xcell,ycell)
61            if(j==0) cycle  !!!This is a global search method to find the particle closet to the missing point

              dx=xg-xp(j)  
              dy=yg-yp(j)                     
              dr=dx*dx+dy*dy

              if(dr<dis_min)then
                  dis_min=dr
                  id_WENO=j
              endif

22            j = celldata(j)
              goto 61 
          enddo
      enddo


      end