      subroutine Createfolder(inquire_path,create_path)
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