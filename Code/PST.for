      subroutine PST  
!----------------------------------------------------------------------
!     Subroutine of particle shifting technique (see Sun, 2017) 
!     Sun, P. N., Colagrossi, A., Marrone, S., & Zhang, A. M. (2017). 
!     The ¦Äplus-SPH model: Simple procedures for a further improvement of the SPH scheme. Computer Methods in Applied Mechanics and Engineering, 315, 25-49.

      use global      
      implicit none
      include 'param.inc' 

      integer i,ii,j
      real*8  dx,dy,r
      real*8  tw,tdwdx,tdwdy
      integer ngridx,ngridy,ngridz
      real*8  dgeomx,dgeomy,dgeomz
      integer minxcell,minycell,minzcell,maxxcell,maxycell,maxzcell 
      integer dnxgcell,dnygcell,dnzgcell,dpxgcell,dpygcell,dpzgcell
      integer xcell,ycell,zcell,xxcell,yycell,zzcell
      real*8  r1,mhsml,dr,V_j,V_i
      real*8  b11_c,b21_c
      real*8  temp,factor,x_xiu,y_xiu,norx_temp,nory_temp


      call init_grid(ngridx,ngridy,dgeomx,dgeomy)

      do i=1,ntotal+nvirt

      call grid_geom(i,xp(i),yp(i),ngridx,ngridy,
     &dgeomx,dgeomy,xxcell,yycell)  
      xgcell(i)=xxcell          
      ygcell(i)=yycell          
      celldata(i) = grid(xxcell,yycell) 
      grid(xxcell,yycell) = i 

      enddo

      !!********calculate the shifting value for all particles except the boundary particles.
      do i=1,ntotal

      b11_c=0.0
      b21_c=0.0

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
111           if(j==0) cycle  
              dx=xp(i)-xp(j)  
              dy=yp(i)-yp(j)
              dr=dx*dx+dy*dy
              r=sqrt(dr)
              mhsml=h(i)
              r1=scale_k*mhsml

              if(r<r1)then
                  v_j=mass(j)/rho(j)
                  call kernel(r,dx,dy,mhsml,tw,tdwdx,tdwdy) 

                  b11_c =b11_c + tdwdx*v_j  *(tw/wdetd)
                  b21_c =b21_c + tdwdy*v_j  *(tw/wdetd)

              endif

222           j = celldata(j)
              goto 111 
          enddo
      enddo


      shift_x(i) =-b11_c
      shift_y(i) =-b21_c
      shift_z(i) =0.0

      shift_x(i)=shift_x(i)*0.1* 2.0*(para_h2dx*detd)**2   
      shift_y(i)=shift_y(i)*0.1* 2.0*(para_h2dx*detd)**2 

      enddo

      end