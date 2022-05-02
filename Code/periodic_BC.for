      subroutine periodic_BC
!----------------------------------------------------------------------
!     Subroutine to implement the periodic boundary 
      use global
      implicit none
      include 'param.inc'
      integer i
      real max_x, min_x, max_y, min_y

      !*****max_x, min_x, max_y, min_y are boundary lines or mirror lines. 
      max_x=1.0
      min_x=0.0
      max_y=1.0
      min_y=0.0

      nvirt=0
      do i=1,ntotal           
          if(xp(i).le. (min_x+10*h_max)) then !10*h_max is the domain of the periodic boundary. The number 10 can be changed.
              nvirt=nvirt+1                     
              xp(ntotal+nvirt)=xp(i)+max_x+detd
              yp(ntotal+nvirt)=yp(i)
              zp(ntotal+nvirt)=zp(i)

              p(ntotal+nvirt)=p(i)
              rho(ntotal+nvirt)=rho(i)
              mass(ntotal+nvirt)=mass(i)
              c(ntotal+nvirt)=c(i)
              vx(ntotal+nvirt)=vx(i)
              vy(ntotal+nvirt)=vy(i)
              vz(ntotal+nvirt)=vz(i)
              vp(ntotal+nvirt)=vp(i)
              h(ntotal+nvirt)=h(i)

              Matrix_B11(ntotal+nvirt)=0.0_8
              Matrix_B12(ntotal+nvirt)=0.0_8
              Matrix_B21(ntotal+nvirt)=0.0_8
              Matrix_B22(ntotal+nvirt)=0.0_8
              itype(ntotal+nvirt)=2
          elseif(xp(i).ge. (max_x-10*h_max))then
              nvirt=nvirt+1                     
              xp(ntotal+nvirt)=xp(i)-max_x-detd
              yp(ntotal+nvirt)=yp(i)
              zp(ntotal+nvirt)=zp(i)

              p(ntotal+nvirt)=p(i)
              rho(ntotal+nvirt)=rho(i)
              mass(ntotal+nvirt)=mass(i)
              c(ntotal+nvirt)=c(i)
              vx(ntotal+nvirt)=vx(i)
              vy(ntotal+nvirt)=vy(i)
              vz(ntotal+nvirt)=vz(i)
              vp(ntotal+nvirt)=vp(i)
              h(ntotal+nvirt)=h(i)

              Matrix_B11(ntotal+nvirt)=0.0_8
              Matrix_B12(ntotal+nvirt)=0.0_8
              Matrix_B21(ntotal+nvirt)=0.0_8
              Matrix_B22(ntotal+nvirt)=0.0_8
              itype(ntotal+nvirt)=2

          endif
      enddo

      do i=1,ntotal+nvirt           
          if(yp(i).le.(min_y+10*h_max)) then
              nvirt=nvirt+1                     
              xp(ntotal+nvirt)=xp(i)
              yp(ntotal+nvirt)=yp(i)+max_y+detd
              zp(ntotal+nvirt)=zp(i)

              p(ntotal+nvirt)=p(i)
              rho(ntotal+nvirt)=rho(i)
              mass(ntotal+nvirt)=mass(i)
              c(ntotal+nvirt)=c(i)
              vx(ntotal+nvirt)=vx(i)
              vy(ntotal+nvirt)=vy(i)
              vz(ntotal+nvirt)=vz(i)
              vp(ntotal+nvirt)=vp(i)
              h(ntotal+nvirt)=h(i)

              Matrix_B11(ntotal+nvirt)=0.0_8
              Matrix_B12(ntotal+nvirt)=0.0_8
              Matrix_B21(ntotal+nvirt)=0.0_8
              Matrix_B22(ntotal+nvirt)=0.0_8
              itype(ntotal+nvirt)=2
          elseif(yp(i).ge.(max_y-10*h_max)) then
              nvirt=nvirt+1                     
              xp(ntotal+nvirt)=xp(i)
              yp(ntotal+nvirt)=yp(i)-max_y-detd
              zp(ntotal+nvirt)=zp(i)

              p(ntotal+nvirt)=p(i)
              rho(ntotal+nvirt)=rho(i)
              mass(ntotal+nvirt)=mass(i)
              c(ntotal+nvirt)=c(i)
              vx(ntotal+nvirt)=vx(i)
              vy(ntotal+nvirt)=vy(i)
              vz(ntotal+nvirt)=vz(i)
              vp(ntotal+nvirt)=vp(i)
              h(ntotal+nvirt)=h(i)

              Matrix_B11(ntotal+nvirt)=0.0_8
              Matrix_B12(ntotal+nvirt)=0.0_8
              Matrix_B21(ntotal+nvirt)=0.0_8
              Matrix_B22(ntotal+nvirt)=0.0_8
              itype(ntotal+nvirt)=2
          endif
      enddo


      end