      subroutine kernel(r,dx,dy,mhsml,tw,tdwdx,tdwdy)   
!----------------------------------------------------------------------
!     Subroutine to create the improved Gaussian kernel
      use global
      implicit none
      include 'param.inc'

      real*8 tw
      real*8 q, factor
      real*8 tdwdx,tdwdy,dx,dy,r,mhsml

      q = r/mhsml 
      tw = 0.e0_8
      tdwdx = 0.e0_8
      tdwdy = 0.e0_8

      factor = 1.e0_8 / (mhsml**(dim)*(pi**(dim/2)))
      factor=factor /(1._8-10._8*exp(-9._8)) 

      if(q.ge.0.and.q.le.3.) then
          tw = factor * (exp(-q*q)-exp(-9._8))
          tdwdx = factor *exp(-q*q)* ( -2._8* dx/mhsml/mhsml)
          tdwdy = factor *exp(-q*q)* ( -2._8* dy/mhsml/mhsml)
      else
          tw = 0._8
          tdwdx = 0._8
          tdwdy = 0._8
      endif


      return	
      end
