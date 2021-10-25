!
!***************************************************************************************
!
!> main.f90
!!
!! Description: Main driver routine for CDE (Chemical Discovery Engine).
!!              Can perform a variety of simulations, including NEB, CINEB,
!!              and Graph-driven reaction-path sampling.
!!
!! Lead author: Scott Habershon.
!!
!! Start date: January 2015
!!
!! Usage: cde.x input
!!
!***************************************************************************************
!
Program CDE
  use constants
  use globaldata
  use rpath
  use IO
  use pathopt
  use pes
  use pathfinder
  implicit none
  type(rxp) :: rp


  ! Read the input file.
  !
  Call ReadInput()


  ! Initialize PES calculations.
  !
  Call SetupEnergyCalc(PEStype,PESfile,PESExecutable)
  Call SetupGeomOpt(PESopttype,PESOptfile,PESoptexecutable)

  ! Start random numbers.
  idum = SetRanSeed(irun)

  ! Switch based on calctype.
  !
  select case (calctype)

     !***************************************************
     ! PATH OPTIMIZATION
     !***************************************************
     !
  case('optpath')

     write(logfile,'("* Starting path optimization..."/)')

     ! Initialize the reaction path object.
     !
     Call NewPath(rp, startfrompath, startfile, endfile, pathfile, nimage, &
          pathinit, .FALSE., idum)

     ! Print the initial path.
     !
     Call PrintPathToFile(rp,'initial_path.xyz')
     write(logfile,'("- Initial path written to initial_path.xyz")')

     ! If stripinactive is true, then we strip out the inactive molecules (i.e. those not
     ! involved in the reaction).
     !
     if (stripinactive) then
        write(logfile,'("*** Stripping inactive molecules from path...")')
        Call StripInactiveFromPath( rp, 'strip_path.xyz', FixedDOF, FixedAtom, NDOFconstr, Natomconstr,.true.)
        write(logfile,'("- Stripped initial path written to strip_path.xyz")')

        ! Now delete the original path and read the new one.
        !
        Call DeletePath(rp)

        ! ..and create a new path by reading from file...
        !
        Call NewPath(rp, .TRUE., startfile, endfile, 'strip_path.xyz', nimage, &
             pathinit, .FALSE., idum)

     endif


     ! Initialize constraints in the path.
     !
     Call SetPathConstraints(rp, NDOFconstr, FixedDOF, Natomconstr, Fixedatom)

     ! Now optimize the path.
     !
     Call OptimizePath( rp )

     ! Tidy up....
     Call DeletePath(rp)


     !***************************************************
     ! DOUBLE-ENDED PATH-FINDING CALCULATION
     !***************************************************
     !
  case ('pathfind')

     write(logfile,'("* Starting double-ended mechanism-finding calculation..."/)')
     Call RunPathFinder()

     !***************************************************
     ! SINGLE-ENDED NETWORK GENERATION CALCULATION
     !***************************************************
     !
  case ('netgrow')

     write(logfile,'("* Starting single-ended network generation calculation..."/)')
     Call RunNetGrow()




  case default
    print*, '* Unknown calculation type', calctype
    stop
  end select


  ! Finish things up.
  !
  write(logfile,'(/"================================================================")')
  write(logfile,'("                     CALCULATION COMPLETE                       ")')
  write(logfile,'("================================================================"/)')
  call flush(logfile)

end Program CDE
