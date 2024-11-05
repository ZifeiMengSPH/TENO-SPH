      subroutine grid_geom (i,xg,yg,ngridx,ngridy,dgeomx,dgeomy,xxcell,yycell)
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