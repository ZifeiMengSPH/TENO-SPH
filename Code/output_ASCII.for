      subroutine output_ASCII
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
      !!subroutine of the ensight data in ASCII form
      use global
      implicit none     
      include 'param.inc'

      character(len=3)::ctemp
      integer k1,k2,k3,k4,k5,i,k


      k=1
      k1=0
      k2=0
      k3=0
      k4=0
      k5=0

      do i=1,ntotal+nvirt
          if(itype(i).eq.-2.or.itype(i).eq.2)then      
              k1=k1+1
              part_id(i)=1
          elseif(itype(i)==3)then
              k2=k2+1
              part_id(i)=2
              !elseif(itype(i)==-4)then 
              !    k3=k3+1
              !    part_id(i)=3
          elseif(itype(i)==4)then
              k4=k4+1
              part_id(i)=4
          elseif(itype(i).eq.1)then    
              k5=k5+1
              part_id(i)=5
          endif
      enddo

      mm=mm+1
      write(ctemp,'(I3.3)') mm

! ***************************     geometry node    ******************************
      if(geometry_node)  then
          k=1
          open(mm,file='./data/geometry.geo_'//trim(adjustl(ctemp)))
          write(mm,"(A4,A10,e12.5e2)")"****",'totaltime=', time
          write(mm,"(A25)")'this is the geometry file'
          write(mm,"(A13)")'node id given'
          write(mm,"(A16)")'element id given'
          write(mm,"(A7)")'extents'
          write(mm,13)x_mingeom,x_maxgeom,y_mingeom,y_maxgeom,z_mingeom,z_maxgeom
   13     format(2(e12.5e2,e12.5e2/),(e12.5e2,e12.5e2))


          if(k1.gt.0)then
              write(mm,"(A4)")'part'
              write(mm,10) k
              write(mm,"(A17)") 'this is for part1'
              write(mm,"(A11)") 'coordinates'
              write(mm,10) k1 

              do i=1,ntotal+nvirt
                  if(part_id(i)==1) write(mm,10) i
              enddo    	

              do i=1,ntotal+nvirt
                  if(part_id(i)==1) write(mm,11) xp(i)
              enddo

              do i=1,ntotal+nvirt
                  if(part_id(i)==1) write(mm,11) yp(i)	
              enddo

              do i=1,ntotal+nvirt
                  if(part_id(i)==1) write(mm,11) zp(i)
              enddo


          endif

          if(k2.gt.0)then
              k=k+1
              write(mm,"(A4)")'part'
              write(mm,10) k
              write(mm,"(A17)") 'this is for part2'
              write(mm,"(A11)") 'coordinates'
              write(mm,10) k2      
              do i=1,ntotal+nvirt     
                  if(part_id(i)==2) write(mm,10) i
              enddo    	

              do i=1,ntotal+nvirt
                  if(part_id(i)==2) write(mm,11) xp(i)	
              enddo
              do i=1,ntotal+nvirt
                  if(part_id(i)==2) write(mm,11) yp(i)	
              enddo
              do i=1,ntotal+nvirt
                  if(part_id(i)==2) write(mm,11) zp(i)
              enddo				
          endif

          if(k3.gt.0)then
              k=k+1
              write(mm,"(A4)")'part'
              write(mm,10) k
              write(mm,"(A17)") 'this is for part3'
              write(mm,"(A11)") 'coordinates'
              write(mm,10) k3

              do i=1,ntotal+nvirt     
                  if(part_id(i)==3)	
     &            write(mm,10) i
              enddo    	

              do i=1,ntotal+nvirt     
                  if(part_id(i)==3) 
     &            write(mm,11) xp(i)	
              enddo
              do i=1,ntotal+nvirt     
                  if(part_id(i)==3) 
     &            write(mm,11) yp(i)	
              enddo
              do i=1,ntotal+nvirt     
                  if(part_id(i)==3)
     &            write(mm,11) zp(i)
              enddo


          endif

          if(k4.gt.0)then
              k=k+1
              write(mm,"(A4)")'part'
              write(mm,10) k
              write(mm,"(A17)") 'this is for part4'
              write(mm,"(A11)") 'coordinates'
              write(mm,10) k4

              do i=1,ntotal+nvirt     
                  if(part_id(i)==4)	
     &            write(mm,10) i
              enddo    	

              do i=1,ntotal+nvirt     
                  if(part_id(i)==4) 
     &            write(mm,11) xp(i)	
              enddo
              do i=1,ntotal+nvirt     
                  if(part_id(i)==4) 
     &            write(mm,11) yp(i)	
              enddo
              do i=1,ntotal+nvirt     
                  if(part_id(i)==4)
     &            write(mm,11) zp(i)
              enddo

          endif

          if(k5.gt.0)then
              k=k+1
              write(mm,"(A4)")'part'
              write(mm,10) k
              write(mm,"(A17)") 'this is for part5'
              write(mm,"(A11)") 'coordinates'
              write(mm,10) k5

              do i=1,ntotal+nvirt     
                  if(part_id(i)==5)	
     &            write(mm,10) i
              enddo    	

              do i=1,ntotal+nvirt     
                  if(part_id(i)==5) 
     &            write(mm,11) xp(i)	
              enddo
              do i=1,ntotal+nvirt     
                  if(part_id(i)==5) 
     &            write(mm,11) yp(i)	
              enddo
              do i=1,ntotal+nvirt     
                  if(part_id(i)==5)
     &            write(mm,11) zp(i)
              enddo

          endif

          close(mm)
      endif		

      ! ***************************     geometry ele    ******************************

      if(geometry_ele)then 
          k=1
          open(mm,file='./data/geometry.geo_'//trim(adjustl(ctemp)))
          write(mm,"(A4,A10,e12.5e2)")"****",'totaltime=', time
          write(mm,"(A25)")'this is the geometry file'
          write(mm,"(A13)")'node id given'
          write(mm,"(A16)")'element id given'
          write(mm,"(A7)")'extents'
          write(mm,13)x_mingeom,x_maxgeom,z_mingeom,z_maxgeom  


          if(k1.gt.0)then
              write(mm,"(A4)")'part'
              write(mm,10) k
              write(mm,"(A17)") 'this is for part1'
              write(mm,"(A11)") 'coordinates'
              write(mm,10) k1

              do i=1,ntotal+nvirt
                  if(part_id(i)==1) write(mm,10) i
              enddo    	

              do i=1,ntotal+nvirt
                  if(part_id(i)==1) write(mm,11) xp(i)
              enddo

              do i=1,ntotal+nvirt
                  if(part_id(i)==1)  write(mm,11) yp(i)	
              enddo

              do i=1,ntotal+nvirt
                  if(part_id(i)==1) write(mm,11) zp(i)	
              enddo

              write(mm,"(A5)") 'point'
              write(mm,10) k1

              do i=1,ntotal+nvirt
                  if(part_id(i)==1) write(mm,10) i
              enddo  

              do i=1,k1
                  write(mm,10) i
              enddo				
          endif

          if(k2.gt.0)then
              k=k+1
              write(mm,"(A4)")'part'
              write(mm,10) k
              write(mm,"(A17)") 'this is for part2'
              write(mm,"(A11)") 'coordinates'
              write(mm,10) k2

              do i=1,ntotal+nvirt     
                  if(part_id(i)==2) write(mm,10) i
              enddo    	

              do i=1,ntotal+nvirt
                  if(part_id(i)==2) write(mm,11) xp(i)	
              enddo

              do i=1,ntotal+nvirt
                  if(part_id(i)==2) write(mm,11) yp(i)	
              enddo
              do i=1,ntotal+nvirt
                  if(part_id(i)==2) write(mm,11) zp(i)
              enddo

              write(mm,"(A5)") 'point'
              write(mm,10) k2
              do i=1,ntotal+nvirt
                  if(part_id(i)==2) write(mm,10) i
              enddo
              do i=1,k2
                  write(mm,10) i
              enddo				
          endif

          if(k3.gt.0)then
              k=k+1
              write(mm,"(A4)")'part'
              write(mm,10) k
              write(mm,"(A17)") 'this is for part3'
              write(mm,"(A11)") 'coordinates'
              write(mm,10) k3

              do i=1,ntotal+nvirt     
                  if(part_id(i)==3)
     &            write(mm,10) i
              enddo    	

              do i=1,ntotal+nvirt     
                  if(part_id(i)==3) 
     &            write(mm,11) xp(i)	
              enddo
              do i=1,ntotal+nvirt     
                  if(part_id(i)==3)
     &            write(mm,11) yp(i)	
              enddo
              do i=1,ntotal+nvirt     
                  if(part_id(i)==3) 
     &            write(mm,11) zp(i)
              enddo

              write(mm,"(A5)") 'point'
              write(mm,10) k3
              do i=1,ntotal+nvirt     
                  if(part_id(i)==3) 
     &            write(mm,10) i
              enddo
              do i=1,k3
                  write(mm,10) i
              enddo				
          endif

          if(k4.gt.0)then
              k=k+1
              write(mm,"(A4)")'part'
              write(mm,10) k
              write(mm,"(A17)") 'this is for part4'
              write(mm,"(A11)") 'coordinates'
              write(mm,10) k4

              do i=1,ntotal+nvirt     
                  if(part_id(i)==4)
     &            write(mm,10) i
              enddo    	

              do i=1,ntotal+nvirt     
                  if(part_id(i)==4) 
     &            write(mm,11) xp(i)	
              enddo
              do i=1,ntotal+nvirt     
                  if(part_id(i)==4)
     &            write(mm,11) yp(i)	
              enddo
              do i=1,ntotal+nvirt     
                  if(part_id(i)==4) 
     &            write(mm,11) zp(i)
              enddo

              write(mm,"(A5)") 'point'
              write(mm,10) k4
              do i=1,ntotal+nvirt     
                  if(part_id(i)==4) 
     &            write(mm,10) i
              enddo
              do i=1,k4
                  write(mm,10) i
              enddo				
          endif

          if(k5.gt.0)then
              k=k+1
              write(mm,"(A4)")'part'
              write(mm,10) k
              write(mm,"(A17)") 'this is for part5'
              write(mm,"(A11)") 'coordinates'
              write(mm,10) k5

              do i=1,ntotal+nvirt     
                  if(part_id(i)==5)
     &            write(mm,10) i
              enddo    	

              do i=1,ntotal+nvirt     
                  if(part_id(i)==5) 
     &            write(mm,11) xp(i)	
              enddo
              do i=1,ntotal+nvirt     
                  if(part_id(i)==5)
     &            write(mm,11) yp(i)	
              enddo
              do i=1,ntotal+nvirt     
                  if(part_id(i)==5) 
     &            write(mm,11) zp(i)
              enddo

              write(mm,"(A5)") 'point'
              write(mm,10) k5
              do i=1,ntotal+nvirt     
                  if(part_id(i)==5) 
     &            write(mm,10) i
              enddo
              do i=1,k5
                  write(mm,10) i
              enddo				
          endif

          close(mm)		
      endif

!***********************     pressure  element   ******************************	

      if(pressure_ele_out) then
          k=1	  
          open(mm+3,file='./data/pressure.Esca_'//trim(adjustl(ctemp)))
          write(mm+3,"(A14,A10,e12.5e2)")'*pressurefile*','totaltime=', time
          if(k1.gt.0)then	 
              write(mm+3,"(A4)")'part' 
              write(mm+3,10) k	  
              write(mm+3,"(A5)") 'point'
              do i=1,ntotal+nvirt
                  if(part_id(i)==1) write(mm+3,11) p(i)              
              enddo	
          endif

          if(k2.gt.0)then
              k=k+1	 
              write(mm+3,"(A4)")'part' 
              write(mm+3,10) k	  
              write(mm+3,"(A5)") 'point'
              do i=1,ntotal+nvirt
                  if(part_id(i)==2) write(mm+3,11)  p(i)
              enddo	
          endif

          if(k3.gt.0)then
              k=k+1	 
              write(mm+3,"(A4)")'part' 
              write(mm+3,10) k	  
              write(mm+3,"(A5)") 'point'
              do i=1,ntotal+nvirt
                  if(part_id(i)==3) 
     &            write(mm+3,11)  p(i)
              enddo	
          endif

          if(k4.gt.0)then
              k=k+1	 
              write(mm+3,"(A4)")'part' 
              write(mm+3,10) k	  
              write(mm+3,"(A5)") 'point'
              do i=1,ntotal+nvirt
                  if(part_id(i)==4) 
     &            write(mm+3,11)  p(i)
              enddo	
          endif

          if(k5.gt.0)then
              k=k+1	 
              write(mm+3,"(A4)")'part' 
              write(mm+3,10) k	  
              write(mm+3,"(A5)") 'point'
              do i=1,ntotal+nvirt
                  if(part_id(i)==5) 
     &            write(mm+3,11)  p(i)
              enddo	
          endif

          close(mm+3)
      endif	   
!*************************   velocity  element   *********************************

      if(velocity_ele_out) then
          k=1			       
          open(mm+4,file='./data/velocity.Evec_'//trim(adjustl(ctemp))) 
          write(mm+4,"(A14,A10,e12.5e2)")'*velocityfile*','totaltime=', time
          if(k1.gt.0)then
              write(mm+4,"(A4)")'part'
              write(mm+4,10) k
              write(mm+4,"(A5)") 'point'

              do i=1,ntotal+nvirt
                  if(part_id(i)==1) write(mm+4,11) vx(i)
              enddo
              do i=1,ntotal+nvirt
                  if(part_id(i)==1) write(mm+4,11) vy(i)
              enddo
              do i=1,ntotal+nvirt
                  if(part_id(i)==1) write(mm+4,11) vz(i)	
              enddo
          endif	 

          if(k2.gt.0)then
              k=k+1
              write(mm+4,"(A4)")'part'
              write(mm+4,10) k
              write(mm+4,"(A5)") 'point'

              do i=1,ntotal+nvirt
                  if(part_id(i)==2) write(mm+4,11) vx(i)
              enddo
              do i=1,ntotal+nvirt
                  if(part_id(i)==2) write(mm+4,11) vy(i)
              enddo
              do i=1,ntotal+nvirt
                  if(part_id(i)==2) write(mm+4,11) vz(i)
              enddo
          endif

          if(k3.gt.0)then
              k=k+1
              write(mm+4,"(A4)")'part'
              write(mm+4,10) k
              write(mm+4,"(A5)") 'point'

              do i=1,ntotal+nvirt
                  if(part_id(i)==3) 
     &            write(mm+4,11) vx(i)
              enddo
              do i=1,ntotal+nvirt
                  if(part_id(i)==3) 
     &            write(mm+4,11) vy(i)
              enddo
              do i=1,ntotal+nvirt
                  if(part_id(i)==3) 
     &            write(mm+4,11) vz(i)	
              enddo
          endif	  

          if(k4.gt.0)then
              k=k+1
              write(mm+4,"(A4)")'part'
              write(mm+4,10) k
              write(mm+4,"(A5)") 'point'

              do i=1,ntotal+nvirt
                  if(part_id(i)==4) 
     &            write(mm+4,11) vx(i)
              enddo
              do i=1,ntotal+nvirt
                  if(part_id(i)==4) 
     &            write(mm+4,11) vy(i)
              enddo
              do i=1,ntotal+nvirt
                  if(part_id(i)==4) 
     &            write(mm+4,11) vz(i)	
              enddo
          endif

          if(k5.gt.0)then
              k=k+1
              write(mm+4,"(A4)")'part'
              write(mm+4,10) k
              write(mm+4,"(A5)") 'point'

              do i=1,ntotal+nvirt
                  if(part_id(i)==5) 
     &            write(mm+4,11) vx(i)
              enddo
              do i=1,ntotal+nvirt
                  if(part_id(i)==5) 
     &            write(mm+4,11) vy(i)
              enddo
              do i=1,ntotal+nvirt
                  if(part_id(i)==5) 
     &            write(mm+4,11) vz(i)	
              enddo
          endif

          close(mm+4) 
      endif      	  	     


!**************************     rho  element   ******************************	

      if(1) then
          k=1	  
          open(mm+10,file='./data/rho.Esca_'//trim(adjustl(ctemp)))
          write(mm+10,"(A14,A10,e12.5e2)")'*rhofile*','totaltime=', time
          if(k1.gt.0)then	 
              write(mm+10,"(A4)")'part' 
              write(mm+10,10) k	  
              write(mm+10,"(A5)") 'point'
              do i=1,ntotal+nvirt
                  if(part_id(i)==1) write(mm+10,11) real(rho(i))             
              enddo	
          endif

          if(k2.gt.0)then
              k=k+1	 
              write(mm+10,"(A4)")'part' 
              write(mm+10,10) k	  
              write(mm+10,"(A5)") 'point'
              do i=1,ntotal+nvirt
                  if(part_id(i)==2) write(mm+10,11)  real(rho(i))
              enddo	
          endif

          if(k3.gt.0)then
              k=k+1	 
              write(mm+10,"(A4)")'part' 
              write(mm+10,10) k	  
              write(mm+10,"(A5)") 'point'
              do i=1,ntotal+nvirt
                  if(part_id(i)==3) 
     &            write(mm+10,11)  real(rho(i))
              enddo	
          endif

          if(k4.gt.0)then
              k=k+1	 
              write(mm+10,"(A4)")'part' 
              write(mm+10,10) k	  
              write(mm+10,"(A5)") 'point'
              do i=1,ntotal+nvirt
                  if(part_id(i)==4) 
     &            write(mm+10,11)  real(rho(i))
              enddo	
          endif

          if(k5.gt.0)then
              k=k+1	 
              write(mm+10,"(A4)")'part' 
              write(mm+10,10) k	  
              write(mm+10,"(A5)") 'point'
              do i=1,ntotal+nvirt
                  if(part_id(i)==5) 
     &            write(mm+10,11)  real(rho(i))
              enddo	
          endif

          close(mm+10)
      endif	

!**************************     mass  element   ******************************	

      if(1) then
          k=1	  
          open(mm+13,file='./data/mass.Esca_'//trim(adjustl(ctemp)))
          write(mm+13,"(A14,A10,e12.5e2)")'*massfile*','totaltime=', time
          if(k1.gt.0)then	 
              write(mm+13,"(A4)")'part' 
              write(mm+13,10) k	  
              write(mm+13,"(A5)") 'point'
              do i=1,ntotal+nvirt
                  if(part_id(i)==1) write(mm+13,11) real(mass(i))             
              enddo	
          endif

          if(k2.gt.0)then
              k=k+1	 
              write(mm+13,"(A4)")'part' 
              write(mm+13,10) k	  
              write(mm+13,"(A5)") 'point'
              do i=1,ntotal+nvirt
                  if(part_id(i)==2) write(mm+13,11)  real(mass(i))
              enddo	
          endif

          if(k3.gt.0)then
              k=k+1	 
              write(mm+13,"(A4)")'part' 
              write(mm+13,10) k	  
              write(mm+13,"(A5)") 'point'
              do i=1,ntotal+nvirt
                  if(part_id(i)==3) 
     &            write(mm+13,11)  real(mass(i))
              enddo	
          endif

          if(k4.gt.0)then
              k=k+1	 
              write(mm+13,"(A4)")'part' 
              write(mm+13,10) k	  
              write(mm+13,"(A5)") 'point'
              do i=1,ntotal+nvirt
                  if(part_id(i)==4) 
     &            write(mm+13,11)  real(mass(i))
              enddo	
          endif

          if(k5.gt.0)then
              k=k+1	 
              write(mm+13,"(A4)")'part' 
              write(mm+13,10) k	  
              write(mm+13,"(A5)") 'point'
              do i=1,ntotal+nvirt
                  if(part_id(i)==5) 
     &            write(mm+13,11)  real(mass(i))
              enddo	
          endif

          close(mm+13)
      endif	

!!*************************   shifting direction   *********************************
!
      if(1) then
          k=1			       
          open(mm+18,file='./data/shift.Evec_'//trim(adjustl(ctemp))) 
          write(mm+18,"(A14,A10,e12.5e2)")'*shiftfile*','totaltime=', time
          if(k1.gt.0)then
              write(mm+18,"(A4)")'part'
              write(mm+18,10) k
              write(mm+18,"(A5)") 'point'

              do i=1,ntotal+nvirt
                  if(part_id(i)==1) write(mm+18,11) shift_x(i)
              enddo
              do i=1,ntotal+nvirt
                  if(part_id(i)==1) write(mm+18,11) shift_y(i)
              enddo
              do i=1,ntotal+nvirt
                  if(part_id(i)==1) write(mm+18,11) shift_z(i)	
              enddo
          endif	 

          if(k2.gt.0)then
              k=k+1
              write(mm+18,"(A4)")'part'
              write(mm+18,10) k
              write(mm+18,"(A5)") 'point'

              do i=1,ntotal+nvirt
                  if(part_id(i)==2) write(mm+18,11) shift_x(i)
              enddo
              do i=1,ntotal+nvirt
                  if(part_id(i)==2) write(mm+18,11) shift_y(i)
              enddo
              do i=1,ntotal+nvirt
                  if(part_id(i)==2) write(mm+18,11) shift_z(i)
              enddo
          endif

          if(k3.gt.0)then
              k=k+1
              write(mm+18,"(A4)")'part'
              write(mm+18,10) k
              write(mm+18,"(A5)") 'point'

              do i=1,ntotal+nvirt
                  if(part_id(i)==3) 
     &            write(mm+18,11) shift_x(i)
              enddo
              do i=1,ntotal+nvirt
                  if(part_id(i)==3) 
     &            write(mm+18,11) shift_y(i)
              enddo
              do i=1,ntotal+nvirt
                  if(part_id(i)==3) 
     &            write(mm+18,11) shift_z(i)	
              enddo
          endif	

          if(k4.gt.0)then
              k=k+1
              write(mm+18,"(A4)")'part'
              write(mm+18,10) k
              write(mm+18,"(A5)") 'point'

              do i=1,ntotal+nvirt
                  if(part_id(i)==4) 
     &            write(mm+18,11) shift_x(i)
              enddo
              do i=1,ntotal+nvirt
                  if(part_id(i)==4) 
     &            write(mm+18,11) shift_y(i)
              enddo
              do i=1,ntotal+nvirt
                  if(part_id(i)==4) 
     &            write(mm+18,11) shift_z(i)	
              enddo
          endif	

          close(mm+18) 
      endif    

!*************************   c    *********************************
      if(1) then
          k=1	  
          open(mm+21,file='./data/c.Esca_'//trim(adjustl(ctemp)))
          write(mm+21,"(A14,A10,e12.5e2)")'*cfile*','totaltime=', time
          if(k1.gt.0)then	 
              write(mm+21,"(A4)")'part' 
              write(mm+21,10) k	  
              write(mm+21,"(A5)") 'point'
              do i=1,ntotal+nvirt
                  if(part_id(i)==1) write(mm+21,11) real(c(i))             
              enddo	
          endif

          if(k2.gt.0)then
              k=k+1	 
              write(mm+21,"(A4)")'part' 
              write(mm+21,10) k	  
              write(mm+21,"(A5)") 'point'
              do i=1,ntotal+nvirt
                  if(part_id(i)==2) write(mm+21,11)  real(c(i))
              enddo	
          endif

          if(k3.gt.0)then
              k=k+1	 
              write(mm+21,"(A4)")'part' 
              write(mm+21,10) k	  
              write(mm+21,"(A5)") 'point'
              do i=1,ntotal+nvirt
                  if(part_id(i)==3) 
     &            write(mm+21,11)  real(c(i))
              enddo	
          endif

          if(k4.gt.0)then
              k=k+1	 
              write(mm+21,"(A4)")'part' 
              write(mm+21,10) k	  
              write(mm+21,"(A5)") 'point'
              do i=1,ntotal+nvirt
                  if(part_id(i)==4) 
     &            write(mm+21,11)  real(c(i))
              enddo	
          endif

          if(k5.gt.0)then
              k=k+1	 
              write(mm+21,"(A4)")'part' 
              write(mm+21,10) k	  
              write(mm+21,"(A5)") 'point'
              do i=1,ntotal+nvirt
                  if(part_id(i)==5) 
     &            write(mm+21,11)  real(c(i))
              enddo	
          endif

          close(mm+21)
      endif	


!!*************************   unevenness   *********************************
      if(1) then
          k=1	  
          open(mm+27,file='./data/unevenness.Esca_'//trim(adjustl(ctemp)))
          write(mm+27,"(A14,A10,e12.5e2)")'*unevennessfile*','totaltime=', time
          if(k1.gt.0)then	 
              write(mm+27,"(A4)")'part' 
              write(mm+27,10) k	  
              write(mm+27,"(A5)") 'point'
              do i=1,ntotal+nvirt
                  if(part_id(i)==1) write(mm+27,11) real(unevenness(i))             
              enddo	
          endif

          if(k2.gt.0)then
              k=k+1	 
              write(mm+27,"(A4)")'part' 
              write(mm+27,10) k	  
              write(mm+27,"(A5)") 'point'
              do i=1,ntotal+nvirt
                  if(part_id(i)==2) write(mm+27,11)  real(unevenness(i))
              enddo	
          endif

          if(k3.gt.0)then
              k=k+1	 
              write(mm+27,"(A4)")'part' 
              write(mm+27,10) k	  
              write(mm+27,"(A5)") 'point'
              do i=1,ntotal+nvirt
                  if(part_id(i)==3) 
     &            write(mm+27,11)  real(unevenness(i))
              enddo	
          endif
          close(mm+27)
      endif	


      !**********************************************create case file
      open(mm+40,file='./data/11output.case')
      write(mm+40,*) 'FORMAT' 
      write(mm+40,*) 'type:			ensight gold'
      write(mm+40,*) 'GEOMETRY'
      write(mm+40,*) 'model:			geometry.geo_***'
      write(mm+40,*) 'VARIABLE'
      write(mm+40,*) 'scalar per element:   1 pressure pressure.Esca_***'
      write(mm+40,*) 'scalar per element:   1 rho rho.Esca_***'
      write(mm+40,*) 'scalar per element:   1 mass mass.Esca_***'
      write(mm+40,*) 'scalar per element:   1 c c.Esca_***'
      write(mm+40,*) 'scalar per element:   1 unevenness unevenness.Esca_***'
      write(mm+40,*) 'vector per element:   1 velocity velocity.Evec_***'
      write(mm+40,*) 'vector per element:   1 shift shift.Evec_***'
      write(mm+40,*) 'TIME'
      write(mm+40,*) 'time set:              1'
      write(mm+40,*) 'number of steps:',mm
      write(mm+40,*) 'filename start number: 1'
      write(mm+40,*) 'filename increment:    1'
      write(mm+40,*) 'time values:           1'
      do i=2,mm
          write(mm+40,*) '              ',i
      enddo
      close(mm+40)


  10  format(I10)      
  11  format(e12.5e2)

      end    