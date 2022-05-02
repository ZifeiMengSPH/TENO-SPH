      subroutine single_step
!----------------------------------------------------------------------
!     Subroutine to calculate the pressure and update the state of governing equations 
      use global
      implicit none
      include 'param.inc'

      integer i
      real*8 det_r
      integer num_tnt
      real*8 mass_total, energy_total, volume_total
      real*8 vp_max
      real*8 pmax_water

      call periodic_BC

      do i=nstart_water,ntotal!+nvirt !calculate the pressure of water
          call water_Tait(p(i),rho(i),c(i))
      enddo
      
      !!!!call TENO-SPH scheme to update the governing equations
      call link_teno5

      do i=1,ntotal!+nvirt
          dvxdt(i) = 0.
          dvydt(i) = 0. 
          dvzdt(i) = 0. 
      enddo


      do i=nair,ntotal!+nvirt
          dvxdt(i) = indvxdt(i) + avdvxdt(i)+grx         
          dvydt(i) = indvydt(i) + avdvydt(i)+gry  
      enddo

      return
      end
