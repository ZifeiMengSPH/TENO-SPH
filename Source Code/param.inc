!---------------------------All parameters of this program

!---------------------------maxmum particle numbers and particle interactions
      integer maxn,max_interaction
      parameter ( maxn    = 400000  , 
     &             max_interaction = 50 * maxn )
!---------------------------Parameters of the linked-list search method
	   integer maxngx,maxngy,maxngz
	  parameter ( maxngx  = 4000      , !
     &            maxngy  = 600      ,
     &            maxngz  = 600       )
!------------------------------------------------------------------------
      integer print_step
      parameter ( print_step = 1)
!---------------------------Parameter of the radius of the support domain  
      real*8 scale_k
      parameter(  scale_k=3.0         )  
!------------------------------------------------------------------------        
      real*8,parameter:: pi = 3.14159265358979323846
!---------------------------Parameters of the output_ASCII.for and the output_Binary.for
      logical geometry_node,geometry_ele,pressure_node_out,
     &        velocity_node_out, pressure_ele_out,velocity_ele_out,  
     &         Mises_out,PEEQ_out,xyzp_out,out_Binary
              
      parameter (out_Binary  = .false. )
      parameter (geometry_node  = .false. )
      parameter (geometry_ele   = .true. )
      parameter (pressure_node_out  = .false. ) 
      parameter (velocity_node_out  = .false. ) 
      parameter (pressure_ele_out  = .true. ) 
      parameter (velocity_ele_out  = .true. ) 
      parameter (Mises_out  = .false. ) 
      parameter (PEEQ_out  = .false. )
      parameter (xyzp_out  = .false. ) 
	

