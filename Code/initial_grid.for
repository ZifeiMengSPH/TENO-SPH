      subroutine init_grid (ngridx,ngridy,
     &dgeomx,dgeomy)
!----------------------------------------------------------------------
!     Subroutine to establish the initial grid for linked list search method (see Liu & Liu 2003) 
!     Liu, G.R., Liu, M.B., 2003. Smoothed Particle Hydrodynamics: A Meshfree Particle Method. World Scientific

      use global
      implicit none
      include 'param.inc'
      integer ngridx,ngridy,ngridz
      real*8 dgeomx,dgeomy,dgeomz,hh2
      integer i,j,k

      ngridx = 1
      ngridy = 1    

      dgeomx = x_maxgeom- x_mingeom
      dgeomy = y_maxgeom- y_mingeom

      hh2=scale_k*h_max
      ngridx=min(int(dgeomx/hh2) + 1,maxngx)      
      ngridy=min(int(dgeomy/hh2) + 1,maxngy)


      do i=1,ngridx
          do j=1,ngridy
              grid(i,j)=0
          enddo
      enddo

      return
      end
