      subroutine output_Binary
      !!subroutine of the ensight data in binary form
      use global
      implicit none     
      include 'param.inc'


      character(len=4)::ctemp
      integer k1,k2,k3,k4,k5,i,k

      character(len=80)::char1,char2,char3,char4,char5,char6,char7,char8,char9,char10,char11,char12,char13,char14,char15,char16
      character(len=80)::char112

      char1="C Binary"
      char2="****"
      char3="geometry file"
      char4="node id given"
      char5="element id given"
      char6="extents"
      char7="part"
      char8="part1"
      char9="coordinates"
      char10="part2"
      char11="part3"
      char112="part4"
      char12="point"
      char13="vortex file"
      char14="velocity file"
      char15="rho file"
      char16="pressure file"



      k=1
      k1=0
      k2=0
      k3=0
      k4=0
      k5=0

      !the label for different parts, such as itype(i)=3 corresponding to water part
      do i=1,ntotal+nvirt
          if(itype(i).eq.-2.or.itype(i).eq.2)then       
              k1=k1+1
              part_id(i)=1
          elseif(itype(i)==3)then
              k2=k2+1
              part_id(i)=2
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

!**********************************************************************************************
      if(1)then       
          k=1
          open(mm,file='./data/geometry.geo_'//trim(adjustl(ctemp)),form='binary',recl=80,access='sequential')
          write(mm) char1     !"C Binary"
          write(mm) char2     !"****"
          write(mm) char3     !"geometry file"
          write(mm) char4     !"node id given"
          write(mm) char5     !"element id given"
          write(mm) char6     !"extents"   
          write(mm) real(x_mingeom),real(x_maxgeom),real(y_mingeom),real(y_maxgeom),0.0,0.0

          if(k1.gt.0)then
              write(mm) char7     !"part"
              write(mm) k         !
              write(mm) char8     !"part1"
              write(mm) char9     !"coordinates"
              write(mm) k1        !total number of particles

              do i=1,ntotal+nvirt
                  if(part_id(i)==1) write(mm) i
              enddo    	

              do i=1,ntotal+nvirt
                  if(part_id(i)==1) write(mm) real(xp(i))
              enddo

              do i=1,ntotal+nvirt
                  if(part_id(i)==1)  write(mm) real(yp(i))	
              enddo

              do i=1,ntotal+nvirt
                  if(part_id(i)==1) write(mm) real(zp(i))	
              enddo

              write(mm) char12    !"point"
              write(mm) k1        

              do i=1,ntotal+nvirt
                  if(part_id(i)==1) write(mm) i  
              enddo  
              do i=1,k1
                  write(mm) i         
              enddo				
          endif


          if(k2.gt.0)then
              k=k+1
              write(mm) char7     !"part"
              write(mm) k         !
              write(mm) char10    !"part2"
              write(mm) char9     !"coordinates"
              write(mm) k2        !total number of particles

              do i=1,ntotal+nvirt     
                  if(part_id(i)==2) write(mm) i
              enddo    	

              do i=1,ntotal+nvirt
                  if(part_id(i)==2) write(mm) real(xp(i))
              enddo
              do i=1,ntotal+nvirt
                  if(part_id(i)==2) write(mm) real(yp(i))	
              enddo
              do i=1,ntotal+nvirt
                  if(part_id(i)==2) write(mm) real(zp(i))
              enddo

              write(mm) char12    !"point"
              write(mm) k2
              do i=1,ntotal+nvirt
                  if(part_id(i)==2) write(mm) i
              enddo
              do i=1,k2
                  write(mm) i
              enddo				
          endif

          if(k3.gt.0)then
              k=k+1
              write(mm)char7          !"part"
              write(mm) k             !
              write(mm) char11        !"part3"
              write(mm) char9         !"coordinates"
              write(mm) k3            !total number of particles

              do i=1,ntotal+nvirt     
                  if(part_id(i)==3) write(mm) i
              enddo    	

              do i=1,ntotal +nvirt    
                  if(part_id(i)==3) write(mm) real(xp(i))
              enddo
              do i=1,ntotal +nvirt    
                  if(part_id(i)==3) write(mm) real(yp(i))
              enddo
              do i=1,ntotal+nvirt     
                  if(part_id(i)==3) write(mm) real(zp(i))	
              enddo

              write(mm) char12  !"point"
              write(mm) k3
              do i=1,ntotal+nvirt     
                  if(part_id(i)==3) write(mm) i
              enddo
              do i=1,k3
                  write(mm) i
              enddo				
          endif

          if(k4.gt.0)then
              k=k+1
              write(mm)char7          !"part"
              write(mm) k             !
              write(mm) char112        !"part3"
              write(mm) char9         !"coordinates"
              write(mm) k4            !total number of particles

              do i=1,ntotal+nvirt     
                  if(part_id(i)==4) write(mm) i
              enddo    	

              do i=1,ntotal +nvirt    
                  if(part_id(i)==4) write(mm) real(xp(i))
              enddo
              do i=1,ntotal +nvirt    
                  if(part_id(i)==4) write(mm) real(yp(i))
              enddo
              do i=1,ntotal+nvirt     
                  if(part_id(i)==4) write(mm) real(zp(i))	
              enddo

              write(mm) char12  !"point"
              write(mm) k4
              do i=1,ntotal+nvirt     
                  if(part_id(i)==4) write(mm) i
              enddo
              do i=1,k4
                  write(mm) i
              enddo				
          endif

          close(mm)		
      endif


!******************************************************************************************
      if(1) then    

      k=1	  
      open(mm+3,file='./data/pressure.Esca_'//trim(adjustl(ctemp)),form='binary',recl=80,access='sequential')
      write(mm+3) char16           !"pressure file"                                                 
      if(k1.gt.0)then	                                                                              
          write(mm+3) char7        !"part"                                                          
          write(mm+3) k	                                                                          
          write(mm+3) char12       !"point"                                                         
          do i=1,ntotal+nvirt                                                                        
              if(part_id(i)==1) write(mm+3) real(p(i))                             
          enddo	
      endif

      if(k2.gt.0)then
          k=k+1	 
          write(mm+3)char7        !"part"
          write(mm+3) k	  
          write(mm+3) char12      !"point" 
          do i=1,ntotal+nvirt
              if(part_id(i)==2) write(mm+3) real(p(i))
          enddo	
      endif

      if(k3.gt.0)then
          k=k+1	 
          write(mm+3)char7        !"part"
          write(mm+3) k	  
          write(mm+3) char12      !"point"
          do i=1,ntotal+nvirt
              if(part_id(i)==3) write(mm+3) real(p(i))
          enddo	
      endif

      if(k4.gt.0)then
          k=k+1	 
          write(mm+3)char7        !"part"
          write(mm+3) k	  
          write(mm+3) char12      !"point"
          do i=1,ntotal+nvirt
              if(part_id(i)==4) write(mm+3) real(p(i))
          enddo	
      endif

      close(mm+3)
      endif	   



      if(1) then 
          k=1			       
          open(mm+4,file='./data/velocity.Evec_'//trim(adjustl(ctemp)),form='binary',recl=80,access='sequential') 
          write(mm+4) char14
          if(k1.gt.0)then
              write(mm+4)char7        !"part"
              write(mm+4) k
              write(mm+4) char12      !"point"

              do i=1,ntotal+nvirt
                  if(part_id(i)==1) write(mm+4) real(vx(i))
              enddo
              do i=1,ntotal+nvirt
                  if(part_id(i)==1) write(mm+4) real(vy(i))
              enddo
              do i=1,ntotal+nvirt
                  if(part_id(i)==1) write(mm+4) real(vz(i))
              enddo
          endif	 

          if(k2.gt.0)then
              k=k+1
              write(mm+4)char7        !"part"
              write(mm+4) k
              write(mm+4)char12       !"point"

              do i=1,ntotal+nvirt
                  if(part_id(i)==2) write(mm+4) real(vx(i))
              enddo
              do i=1,ntotal+nvirt
                  if(part_id(i)==2) write(mm+4) real(vy(i))
              enddo
              do i=1,ntotal+nvirt
                  if(part_id(i)==2) write(mm+4) real(vz(i))
              enddo
          endif

          if(k3.gt.0)then
              k=k+1
              write(mm+4)char7    !"part"
              write(mm+4) k
              write(mm+4) char12  !"point"

              do i=1,ntotal+nvirt
                  if(part_id(i)==3)write(mm+4) real(vx(i))
              enddo
              do i=1,ntotal+nvirt
                  if(part_id(i)==3)write(mm+4) real(vy(i))
              enddo
              do i=1,ntotal+nvirt
                  if(part_id(i)==3)write(mm+4) real(vz(i))
              enddo  
          endif	

          if(k4.gt.0)then
              k=k+1
              write(mm+4)char7    !"part"
              write(mm+4) k
              write(mm+4) char12  !"point"

              do i=1,ntotal+nvirt
                  if(part_id(i)==4)write(mm+4) real(vx(i))
              enddo
              do i=1,ntotal+nvirt
                  if(part_id(i)==4)write(mm+4) real(vy(i))
              enddo
              do i=1,ntotal+nvirt
                  if(part_id(i)==4)write(mm+4) real(vz(i))
              enddo  
          endif	

          close(mm+4) 
      endif    



      if(0) then    

      k=1	  
      open(mm+5,file='./data/vortex.Esca_'//trim(adjustl(ctemp)),form='binary',recl=80,access='sequential')
      write(mm+5) char13                                                           
      if(k1.gt.0)then	                                                                              
          write(mm+5) char7                                                                  
          write(mm+5) k	                                                                          
          write(mm+5) char12                                                                
          do i=1,ntotal+nvirt                                                                        
              if(part_id(i)==1) write(mm+5) real(grad_vyx(i)-grad_vxy(i))                            
          enddo	
      endif

      if(k2.gt.0)then
          k=k+1	 
          write(mm+5)char7        
          write(mm+5) k	  
          write(mm+5) char12       
          do i=1,ntotal+nvirt
              if(part_id(i)==2) write(mm+5) real(grad_vyx(i)-grad_vxy(i))
          enddo	
      endif

      if(k3.gt.0)then
          k=k+1	 
          write(mm+5)char7        
          write(mm+5) k	  
          write(mm+5) char12      
          do i=1,ntotal+nvirt
              if(part_id(i)==3) write(mm+5) real(grad_vyx(i)-grad_vxy(i))
          enddo	
      endif

      if(k4.gt.0)then
          k=k+1	 
          write(mm+5)char7        
          write(mm+5) k	  
          write(mm+5) char12      
          do i=1,ntotal+nvirt
              if(part_id(i)==4) write(mm+5) real(grad_vyx(i)-grad_vxy(i))
          enddo	
      endif

      close(mm+5)
      endif	   


      if(0) then    
          k=1	  
          open(mm+6,file='./data/rho.Esca_'//trim(adjustl(ctemp)),form='binary',recl=80,access='sequential')
          write(mm+6) char15                                                           
          if(k1.gt.0)then	                                                                              
              write(mm+6) char7                                                                  
              write(mm+6) k	                                                                          
              write(mm+6) char12                                                                
              do i=1,ntotal+nvirt                                                                        
                  if(part_id(i)==1) write(mm+6) real(rho(i))                             
              enddo	
          endif

          if(k2.gt.0)then
              k=k+1	 
              write(mm+6)char7        
              write(mm+6) k	  
              write(mm+6) char12       
              do i=1,ntotal+nvirt
                  if(part_id(i)==2) write(mm+6) real(rho(i))
              enddo	
          endif

          if(k3.gt.0)then
              k=k+1	 
              write(mm+6)char7        
              write(mm+6) k	  
              write(mm+6) char12      
              do i=1,ntotal+nvirt
                  if(part_id(i)==3) write(mm+6) real(rho(i))
              enddo	
          endif

          if(k4.gt.0)then
              k=k+1	 
              write(mm+6)char7        
              write(mm+6) k	  
              write(mm+6) char12      
              do i=1,ntotal+nvirt
                  if(part_id(i)==4) write(mm+6) real(rho(i))
              enddo	
          endif

          close(mm+6)
      endif	 


      !**********************************************create case file that is opened by ensight software
      open(mm+40,file='./data/11output.case')
      write(mm+40,*) 'FORMAT' 
      write(mm+40,*) 'type:			ensight gold'
      write(mm+40,*) 'GEOMETRY'
      write(mm+40,*) 'model:			geometry.geo_***'
      write(mm+40,*) 'VARIABLE'
      write(mm+40,*) 'scalar per element:   1 pressure pressure.Esca_***'
      write(mm+40,*) 'vector per element:   1 velocity velocity.Evec_***'
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