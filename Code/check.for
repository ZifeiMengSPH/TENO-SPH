      subroutine check_limits
!----------------------------------------------------------------------
!     Subroutine to update the limits of grids to prepare for linked list search method (see Liu & Liu 2003) 
!     Liu, G.R., Liu, M.B., 2003. Smoothed Particle Hydrodynamics: A Meshfree Particle Method. World Scientific
      use global
      implicit none     
      include 'param.inc'
      integer i
      write(*,*) '      check',ntotal+nvirt 

      xmax=maxval(xp(1:ntotal+nvirt))
      xmin=minval(xp(1:ntotal+nvirt))
      
      ymax=maxval(yp(1:ntotal+nvirt))
      ymin=minval(yp(1:ntotal+nvirt))
      
      x_maxgeom=xmax+20.5*h_max
      x_mingeom=xmin-20.5*h_max  
                     
      y_maxgeom=ymax+20.5*h_max
      y_mingeom=ymin-20.5*h_max 


      end