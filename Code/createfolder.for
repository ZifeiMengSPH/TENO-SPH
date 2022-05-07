      subroutine Createfolder(inquire_path,create_path)
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
!----------------------------------------------------------------------
!     Subroutine to create folder
      use IFPORT
      implicit none
      integer(kind=4) :: status,err
      logical(kind=4) :: ierr
      character(len=100) :: inquire_path,create_path

      inquire(DIRECTORY=trim(adjustl(inquire_path)), EXIST=ierr)
      if(ierr) then
          print*,'The needed directory has existed'
      else
          status=SYSTEM('md '//trim(adjustl(create_path)))
          if(status==-1) then
              err=ierrno()
              print*,'Error=',err,'please see the Fortran help document'
              print*,' '
              stop 'Folder creation failed'
          else
              print*,'The needed directory does not exist and creat it successfully'
          endif
      end if

      return
      end subroutine Createfolder