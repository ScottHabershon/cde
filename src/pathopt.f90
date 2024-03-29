

!
!***************************************************************************************
!
!> @brief Contains routines for path optimization, including CINEB.
!
!***************************************************************************************
!
Module pathopt

  use globaldata
  use constants
  use pes
  use chemstr
  use rpath
  implicit none

contains

  !************************************************************************
  !> InterpolatePath
  !!
  !! Interpolates a reaction path, without optimisation.
  !!
  !! - rp: Generated interpolation of reaction path.
  !!
  !************************************************************************
  !
  subroutine InterpolatePath(rp)
    type(rxp), intent(out) :: rp

    integer :: nmol(2)
    logical :: ldum

    ! Initialize the reaction path object.
    Call NewPath(rp, startfrompath, startfile, endfile, pathfile, nimage, &
        pathinit, .FALSE., idum)

    ! Print the initial path.
    Call PrintPathToFile(rp, 'linear_path.xyz')
    write(logfile,'("- Initial path written to linear_path.xyz")')

    ! If stripinactive is true, then we strip out the inactive molecules (i.e. those not
    ! involved in the reaction).
    if (stripinactive) then
      write(logfile, '("*** Stripping inactive molecules from path...")')
      Call StripInactiveFromPath(rp, 'linear_path_stripped.xyz', FixedDOF, FixedAtom, NDOFconstr, Natomconstr, .true.)
      write(logfile, '("- Stripped initial path written to linear_path_stripped.xyz")')

      ! Now delete the original path and read the new one.
      Call DeletePath(rp)
      ! ..and create a new path by reading from file...
      Call NewPath(rp, .TRUE., startfile, endfile, 'linear_path_stripped.xyz', nimage, &
            pathinit, .FALSE., idum)
    endif

    if (idpppath) then
      ! Initialize constraints in the path.
      Call SetPathConstraints(rp, NDOFconstr, FixedDOF, Natomconstr, Fixedatom)

      ! Determine initial graphs and molecules for end-points.
      call GetGraph(rp%cx(1))
      call GetMols(rp%cx(1))
      nmol(1) = rp%cx(1)%nmol
      call GetGraph(rp%cx(rp%nimage))
      call GetMols(rp%cx(rp%nimage))
      nmol(2) = rp%cx(rp%nimage)%nmol

      ! Do we optimize the end-points before NEB?
      if (optendsbefore) then
        call AbInitio(rp%cx(1), 'optg', ldum)
        call AbInitio(rp%cx(rp%nimage), 'optg', ldum)
        call GetGraph(rp%cx(1))
        call GetMols(rp%cx(1))
        call GetGraph(rp%cx(rp%nimage))
        call GetMols(rp%cx(rp%nimage))
        if (nmol(1) .ne.  rp%cx(1)%nmol .or. nmol(2) .ne. rp%cx(rp%nimage)%nmol )  then
          write(6,'("Warning: Number of molecules changed during optimization!: was: ",2I2,", is: ", 2I2)') &
          nmol(1:2), rp%cx(1)%nmol , rp%cx(rp%nimage)%nmol
          write(logfile,'("Warning: Number of molecules changed during optimization!: was: ",2I2,", is: ", 2I2)') &
          nmol(1:2), rp%cx(1)%nmol , rp%cx(rp%nimage)%nmol
        endif
      endif

      ! If required, redraw the internal images.
      if (reconnect) then
        rp%coeff(:,:,:) = 0.d0
        call FourierToPath(rp)
      endif

      ! Recalculate graphs and mols.
      call GetGraph(rp%cx(1))
      call GetMols(rp%cx(1))
      call GetGraph(rp%cx(rp%nimage))
      call GetMols(rp%cx(rp%nimage))

      ! Calculate IDPP path.
      call FindIDPPPath(rp, NEBIter, NEBConv, NEBstep, NEBspring)

      ! Print the IDPP path.
      Call PrintPathToFile(rp,'idpp_path.xyz')
      write(logfile,'("- IDPP-refined path written to idpp_path.xyz")')
    endif

    return
  end subroutine InterpolatePath

  !************************************************************************
  !> OptimizePath
  !!
  !! Initiates a path optimization calculation.
  !!
  !! - rp: Initial reaction path.
  !! - nebiter: number of NEB iterations
  !! - neboutfreq: number of NEB iterations between outputs
  !! - nebspring: NEB spring parameter in atomic units
  !! - pathoptmethod: Method used to optimize reaction path.
  !!                  Allowed values are 'NEB'....
  !! - NEBconv - Convergence threshold on average force for NEB calculations.
  !! - CIthresh - Average force modulus when climbing-image should be activated.
  !! - inputfile: The input filename, used to generate further output files.
  !! - logfile: Log file unit nunmber for output
  !!
  !************************************************************************
  !
  Subroutine OptimizePath(rp)
    implicit none
    double precision :: tmpCIthresh, tmpVSthresh
    type(rxp) :: rp

    ! Optimize the path based on the method in pathoptmethod
    !
    select case(pathoptmethod)

    case('cineb')
      if (anebb .gt. 0) then
        tmpCIthresh = CIthresh ; tmpVSthresh = VSthresh
        CIthresh = 0.0d0    ;    VSthresh = 0.0d0
      endif
      call CINEB(rp)

    case default
      stop '* ERROR: Unknown pathoptmethod in OptimizePath in pathopt.f90'

    end select

    return
  end Subroutine OptimizePath


  !
  !************************************************************************
  !> CINEB
  !!
  !! Runs a climbing-image nudged elastic band calculation starting from
  !! an initial reaction-path defined in rp.
  !!
  !! - rp: The reaction-path object to be refined.
  !!
  !************************************************************************
  !
  Subroutine CINEB(rp)

    implicit none
    logical :: minimize, ci_flag, isom
    type(rxp) :: rp, rp1, rp2, rp_init
    integer :: i, l, m, iter, j, k, idof, n1, n2, itmp, imove, nmol(2), tmpiter(rp%nimage), imax
    integer :: isum, isum2
    real(8) :: lambda, fnorm, f0, fp, fm, alpha, alpha1, alpha2, alpha3, Emax
    real(8) :: f1, f2, f3, factor, a, b, c, alpha_min, test, fb, step
    real(8) :: fnorm1, fnorm2, best, ran2, fmax, U1(3,3)
    real(8) :: aa, mx, gg, xe, maxf, dx, dy, dz, distance
    real(8), allocatable :: v(:,:), F(:,:), alph(:)  !Fire variables
    real(8) :: f_inc, f_dec, f_alp, alph_start, pfire, delt !Fire variables
    real(8), allocatable :: pr(:,:,:), pg(:,:,:), halph(:), hbeta(:) ! heavy-ball variables
    integer :: Nmin, Nmin_count !Fire variables
    character (len=5) :: x1, cdum, string
    logical :: success, shim, variablespring, callgather, ldum, debug, bfrozen(rp%cx(1)%na)
    integer  :: MolecAct(rp%cx(1)%na), nactmol, layers
    double precision  :: PP(rp%na,rp%na), D(3,3), rr(rp%nimage,3,rp%na), energy(rp%nimage), similarity(2)
    integer, allocatable           :: fgraph(:,:,:)
    character(len=2), allocatable  :: flabel(:,:)
    integer,          allocatable  :: fmap(:,:,:), fna(:)
    type(cxs)                      :: tcx(2), tcx2(2), cx1,cx2

    debug = .true.

    ! Determine initial graphs and molecules for end-points.
    !
    call GetGraph(rp%cx(1))
    call GetMols(rp%cx(1))
    nmol(1) = rp%cx(1)%nmol
    call GetGraph(rp%cx(rp%nimage))
    call GetMols(rp%cx(rp%nimage))
    nmol(2) = rp%cx(rp%nimage)%nmol
    if (debug) call PrintPathToFile(rp,'interpol1.xyz' )

    ! Do we optimize the end-points before NEB?
    !
    if (optendsbefore) then
      call AbInitio(rp%cx(1), 'optg', ldum)
      call AbInitio(rp%cx(rp%nimage), 'optg', ldum)
      call GetGraph(rp%cx(1))
      call GetMols(rp%cx(1))
      call GetGraph(rp%cx(rp%nimage))
      call GetMols(rp%cx(rp%nimage))
      if (nmol(1) .ne.  rp%cx(1)%nmol .or. nmol(2) .ne. rp%cx(rp%nimage)%nmol )  then
        write(6,'("Warning: Number of molecules changed during optimization!: was: ",2I2,", is: ", 2I2)') &
        nmol(1:2), rp%cx(1)%nmol , rp%cx(rp%nimage)%nmol
        write(logfile,'("Warning: Number of molecules changed during optimization!: was: ",2I2,", is: ", 2I2)') &
        nmol(1:2), rp%cx(1)%nmol , rp%cx(rp%nimage)%nmol
      endif
    endif

    ! If required, redraw the internal images.
    if (reconnect) then
      rp%coeff(:,:,:) = 0.d0
      call FourierToPath(rp)
    endif

    ! define spring constant between images:
    rp%ks = NEBspring
    variablespring = .false.

    ! Recalculate graphs and mols.
    !
    call GetGraph(rp%cx(1))
    call GetMols(rp%cx(1))
    call GetGraph(rp%cx(rp%nimage))
    call GetMols(rp%cx(rp%nimage))
    if (debug) call PrintPathToFile(rp, 'interpol2.xyz')

    ! If required, generate Image Dependent Pair Potential Path interpolation.
    ! (Smidstrup et al J Chem phys 2014, 140):
    !
    if (idpppath) then
      call FindIDPPPath(rp, NEBIter*250, NEBConv*0.1d0, NEBstep, NEBspring)
    endif

    ! Sort out orientation of images.
    !
    ! Call SortOrientation(rp)
    
    ! Output initial reoriented path:
    !
    call PrintPathToFile(rp, 'initial_reoriented_path.xyz')

    ! Set the range of images to be optimized.
    !
    if (optendsduring) then
       n1 = 1
       n2 = rp%nimage
    else
       n1 = 2
       n2 = rp%nimage - 1
    endif

    ! Output initial info.
    !
    write(logfile,'("* Starting CINEB calculation..."/)')
    write(logfile,'("- NEB force-convergence monitor: ",1x,A)')trim(inputfile)//'.nebconv'
    write(logfile,'("- NEB energy profile outputs: ",1x,A/)')trim(inputfile)//'.nebprofile'

    ! Open NEB convergence monitor file.
    !
    open(18, file = trim( inputfile )//'.nebconv', status = 'unknown')
    write(18,'("# Force convergence during optimization (in au - without and with springs")')
    open(17, file = trim( inputfile )//'.nebprofile', status = 'unknown')
    Write(17,'("# NEB energy profile during optimization (in au)")')
    write(17,'("# Lambda    |    Energy / Eh  | Relative energy / Eh")')
    open(101, file = trim( inputfile )//'.conv', status = 'unknown')
    write(101,'("# Force convergence during optimization (in au - RMS, IMAX")')

    ! Calculate the energy along the initial path.
    !
    call GetPathGradients(rp, success, .true.)

    ! Add restraint forces to end-points if required.
    !
    if (nebrestrend .and. optendsduring .and. (.not. optendsbefore) ) then
      ! If optsduring is true but optsbefore is not, then perform the optimization with some graph constrains
      ! if nebrestrend is true, these parameters are used in a tanch function that switches the graph constraints as the
      ! forces on the RP apprroach convergence NEBConv
      maxf = 0.750d0
      mx = 1.0d0-NEBConv/4.0d0
      xe = maxf
      aa = atanh(2.0*mx-1.0)/(xe-xe/2.0)
      call GraphConstraints(rp%cx(1), gdsrestspring, nbstrength, nbrange, kradius)
      call GraphConstraints(rp%cx(rp%nimage), gdsrestspring, nbstrength, nbrange, kradius)
      rp%cx(1)%vcon = 0.0d0 ; rp%cx(rp%nimage)%vcon = 0.0d0
    endif

    ! Output the energies along the path.
    !
    write(logfile, '("- Calculating initial energy along path...BEFORE ENDPOINT OPTIMIZATION"/)')
    write(logfile, '("# Lambda    |    Energy / Eh  | Relative energy / Eh")')
    open(14, file = trim( inputfile )//'.energy-neb-start', status = 'unknown')
    write(14, '(/"# NEB energy profile at start of NEB run (in au)"/)')
    write(14, '("# Lambda    |    Energy / Eh  | Relative energy / Eh")')
    lambda = 0.0
    do i = 1, rp%nimage
    !  lambda = real(i-1)/real(rp%nimage-1)
      write(14, '(3(1x,f14.8))') lambda, rp%cx(i)%vcalc, (rp%cx(i)%vcalc - rp%cx(1)%vcalc)
      write(17, '(3(1x,f14.8))') lambda, rp%cx(i)%vcalc, (rp%cx(i)%vcalc - rp%cx(1)%vcalc)
      write(logfile, '(3(1x,f14.8))') lambda, rp%cx(i)%vcalc, (rp%cx(i)%vcalc - rp%cx(1)%vcalc)
      if (i .lt. rp%nimage) lambda = lambda + norm2(reshape(rp%cx(i+1)%r - rp%cx(i)%r, (/rp%na*3/)))
    enddo
    close(14)
    write(17,*)
    write(logfile,*)

    ! If outputting every iteration, create the full trajectory file.
    !
    if (NEBoutfreq == 1) then
      call PrintPathToFile(rp, 'full_neb_traj.xyz')
    endif

    ! Calculate initial norm of projected forces, assuming no climbing image...
    !
    select case (projforcetype)
      case(1)
        call GetProjForces1(rp, .false., .FALSE., optendsduring)
      case(2)
        call GetProjForces2(rp, .false., .FALSE., optendsduring)
      case(3)
        call GetProjForces3(rp, .false., .FALSE., optendsduring)
    end select
    call GetForceNorm(rp, fnorm1, fmax, n1, n2)
    select case (projforcetype)
      case(1)
        call GetProjForces1(rp, .true., .FALSE., optendsduring)
      case(2)
        call GetProjForces2(rp, .true., .FALSE., optendsduring)
      case(3)
        call GetProjForces3(rp, .true., .FALSE., optendsduring)
    end select
    call GetForceNorm(rp, fnorm2, fmax, n1, n2)
    write(18,*)' 0 ', fnorm2, fmax
    write(logfile, '("*** Initial force norms: ")')
    write(logfile, '("* |Forces| WITHOUT springs   = ",1x,f14.8," au ")') fnorm1
    write(logfile, '("* |Forces| WITH springs      = ",1x,f14.8," au ")') fnorm2
    write(logfile, '("* Maximum force WITH springs = ",1x,f14.8," au "/)') fmax

    ! Now perform optimization of the path based on the method in NEBmethod:
    !
    select case (NEBmethod)
      case('quickmin')  !< Quickmin, with Euler update.
        write(logfile, '("* Running QUICKMIN optimization..."/)')

        ! Zero initial momenta.
        do i = n1, n2
          rp%cx(i)%p = 0.d0
        enddo

      case('steepest')  !< Simple steepest descent, taking NEBstep * Force step each iteration
        write(logfile, '("* Running STEEPEST DESCENT optimization..."/)')

      case('fire')
        write(logfile, '("* Running FIRE optimization..."/)')
        allocate(v(rp%nimage, rp%na*3), F(rp%nimage, rp%na*3), alph(rp%nimage))
        Nmin = 5
        f_inc = 1.10d0 ; f_dec = 0.50d0 ; f_alp = 0.990d0 ; alph_start = 0.10d0
        delt = NEBstep/10.0 ; alph = alph_start ; Nmin_count = 0
        do i = n1, n2
          rp%cx(i)%p = 0.0d0
          do j = 1, rp%na
            rp%cx(i)%p(1:3, j) = rp%cx(i)%p(1:3, j) + delt * rp%cx(i)%force(1:3, j)/2.0d0
          !  v(i,(j-1)*3+1:(j-1)*3+3) = rp%cx(i)%p(1:3,j)/rp%cx(i)%mass(j)
            v(i, (j-1)*3+1:(j-1)*3+3) = rp%cx(i)%p(1:3, j)
          !   v = 0.d0
          enddo
        enddo
      case default
        print *, '* Unknown NEBmethod in CINEB in pathopt.f90: ', NEBmethod
        stop
    end select

    !****************************************************************************
    ! Loop over optimization iterations.
    !
    ci_flag = .FALSE.
    do iter = 1, NEBiter

      ! Update coordinates according to projected force.
      !
      select case (NEBmethod)
        case('quickmin')  !< Quickmin, with Euler update.
          ! Recalculate projected momenta.
          !
          do i = n1, n2
            call GetProjectedMomenta(rp%cx(i))
            idof = 0
            do j = 1, rp%na
              if (.not. rp%cx(i)%fixedatom(j)) then
                do k = 1, 3
                  idof = idof + 1
                  if (.not. rp%cx(i)%FixedDOF(idof)) then
                    rp%cx(i)%r(k, j) = rp%cx(i)%r(k, j) + NEBstep * ( rp%cx(i)%p(k, j) / rp%cx(i)%mass(j) )
                !   rp%cx(i)%r(k,j) = rp%cx(i)%r(k,j) + NEBstep * ( rp%cx(i)%p(k,j) )
                    Rp%cx(i)%p(k, j) = rp%cx(i)%p(k, j) + NEBStep * rp%cx(i)%force(k, j)
                  endif
                enddo
              else
                idof = idof + 3
              endif
            enddo
          enddo

        case('steepest')  !< Simple steepest descent, taking NEBstep * Force step each iteration
          do i = n1, n2
            idof = 0
            do j = 1, rp%na
              if (.not. rp%cx(i)%fixedatom(j)) then
                do k = 1, 3
                  idof = idof + 1
                  if (.not. rp%cx(i)%FixedDOF(idof)) then
                    rp%cx(i)%r(k, j) = rp%cx(i)%r(k, j) + NEBstep * rp%cx(i)%force(k, j)
                  endif
                enddo
              else
                idof = idof + 3
              endif
            enddo
          enddo

        case('fire')  !< FIRE update
          do i = n1, n2
            F(i, :) = reshape(rp%cx(i)%force, (/rp%na*3/))
            pfire = dot_product(v(i, :), F(i, :))
            if (norm2(F(i, :)) .gt. 0.00001) then
              v(i, :) = (1.0d0-alph(i)) * v(i, :) + alph(i) * F(i, :) * norm2(v(i, :))/norm2(F(i, :))
            endif
            if (pfire > 0.0d0 .and. Nmin_count > Nmin) then
              delt = min(delt * f_inc, 1.0 * NEBstep)
              Alph(i) = alph(i) * f_alp
            elseif (pfire <= 0.0d0) then
              Nmin_count = 0
              delt = delt * f_dec
              v(i, :) = 0.0d0
              rp%cx(i)%p = 0.0d0
              alph(i) = alph_start
            endif
            if (pfire > 0.0d0) Nmin_count = Nmin_count + 1
            do j = 1, rp%na
              rp%cx(i)%r(1:3, j) = rp%cx(i)%r(1:3, j) + delt * v(i, (j-1)*3+1:(j-1)*3+3)
              ! rp%cx(i)%p(1:3,j) = v(i,(j-1)*3+1:(j-1)*3+3)*rp%cx(i)%mass(j) + delt * rp%cx(i)%force(1:3,j)
              !v(i,(j-1)*3+1:(j-1)*3+3) = rp%cx(i)%p(1:3,j) / rp%cx(i)%mass(j)
              rp%cx(i)%p(1:3, j) = v(i, (j-1)*3+1:(j-1)*3+3) + delt * rp%cx(i)%force(1:3, j)
              v(i, (j-1)*3+1:(j-1)*3+3) = rp%cx(i)%p(1:3, j)
            enddo
          enddo

        case default

      end select

      ! Recalculate energy and projected forces.
      !
      call GetPathGradients(rp, success, .false.)

      ! Remove overall translation and rotation.
      ! TBD

      ! Add restraint forces to end-points if required.
      !
      if (nebrestrend .and. optendsduring .and. (.not. optendsbefore)) then
        gg = 0.50d0*(1.0d0+tanh((fnorm1-xe/2.0d0)*aa))
        call GraphConstraints(rp%cx(1), gdsrestspring*gg, nbstrength*gg, nbrange, kradius)
        call GraphConstraints(rp%cx(rp%nimage), gdsrestspring*gg, nbstrength*gg, nbrange, kradius)
        rp%cx(1)%vcon = 0.0d0 ; rp%cx(rp%nimage)%vcon = 0.0d0
      endif

      ! Check convergence - calculate perpendicular forces without spring.
      !
      select case (projforcetype)
        case(1)
          call GetProjForces1(rp, .false., ci_flag, optendsduring)
        case(2)
          call GetProjForces2(rp, .false., ci_flag, optendsduring)
        case(3)
          call GetProjForces3(rp, .false., ci_flag, optendsduring)
        case default
      end select
      call GetForceNorm(rp, fnorm1, fmax, n1, n2)
      select case (projforcetype)
        case(1)
          call GetProjForces1(rp, .true., ci_flag, optendsduring)
        case(2)
          call GetProjForces2(rp, .true., ci_flag, optendsduring)
        case(3)
          call GetProjForces3(rp, .true., ci_flag, optendsduring)
        case default
      end select
      call GetForceNorm(rp, fnorm2, fmax, n1, n2)
      write(18, *) iter, fnorm2, fmax
      write(logfile, '(/"*** Optimization progress summary at iteration: ",1x,i5)') iter
      write(logfile, '("* |Forces| WITHOUT springs =   ",1x,f14.8," au ")') fnorm1
      write(logfile, '("* |Forces| WITH springs =      ",1x,f14.8," au ")') fnorm2
      write(logfile, '("* Maximum force WITH springs = ",1x,f14.8," au "/)') fmax

      ! Determine variable springs, if being used:
      !
      if (variablespring) then
        call VariableSprings(rp)
      endif

      ! Do we turn on CINEB forces?
      !
      if (fnorm2 <= CIthresh) then
        write(logfile, '("* |Forces| < CIthresh :: TURNING ON CLIMBING-IMAGE")')
        ci_flag = .TRUE.
      endif

      ! Do we turn on variable spring strength?
      !
      if (fnorm2 <= VSthresh .and. VSthresh /= 0.0d0 ) then
        write(logfile, '("* |Forces| < VSthresh :: TURNING ON VARIABLE SPRING")')
        variablespring = .TRUE.
      endif

      ! Output data every NEBoutfreq steps.
      !
      if (mod(iter,NEBoutfreq) == 0) then
        lambda = 0.0
        do i = 1, rp%nimage
          write(17, '(3(1x,f14.8))') lambda, rp%cx(i)%vcalc, (rp%cx(i)%vcalc - rp%cx(1)%vcalc)
          if (i .lt. rp%nimage) lambda = lambda + norm2(reshape(rp%cx(i+1)%r-rp%cx(i)%r, (/rp%na*3/)))
        enddo
        write(17, *)
        write (x1,fmt4) iter
        ! If writing every step, just append to a centralised trajectory instead of multiple files.
        if (NEBoutfreq == 1) then
          call PrintPathToFile(rp, 'full_neb_traj.xyz', .true.)
        else
          call PrintPathToFile(rp, trim( inputfile ) //'_'//trim(x1)//'.xyz')
        endif
        write(logfile, '(/"*** Current path output: ITERATION =",1x,i5,1x,":: OUTPUT FILE =",1x,A/)') iter, &
        trim(inputfile) //'_'//trim(x1)//'.xyz'
      endif

      ! Are we converged?
      !
!      if (fnorm1 <= NEBconv .and. fmax <= NEBmaxconv ) then             ! changed to fnorm 2 !
      if (fnorm2 <= NEBconv) then             ! changed to fnorm 2 !
        write(logfile, '("* |Forces| < NEBconv :: CONVERGED")')
        exit
      endif

      ! New convergence, assuming CI-NEB.
      !
      ! Emax = -1d6
      ! do i = n1, n2
      !    if (rp%cx(i)%vcalc > Emax ) then
      !       Emax = rp%cx(i)%vcalc
      !       imax = i
      !    endif
      ! enddo
      ! fnorm1 = 0.d0
      ! isum = 0
      ! isum2 = 0
      ! fnorm2 = 0.d0
      ! do i = n1, n2
      !    idof = 0
      !    do j = 1, rp%na
      !       if (.not.rp%cx(i)%Fixedatom(j)) then
      !          do k = 1, 3
      !             idof = idof + 1
      !             if (.not.rp%cx(i)%FixedDof(idof)) then
      !               if (i /= imax) then
      !                isum = isum + 1
      !                fnorm1 = fnorm1 + rp%cx(i)%force(k,j)**2
      !              else
      !                isum2 = isum2 + 1
      !                fnorm2 = fnorm2 + rp%cx(i)%force(k,j)**2
      !              endif
      !             endif
      !          enddo
      !       else
      !          idof = idof + 3
      !       endif
      !    enddo
      ! enddo
      ! fnorm1 = dsqrt(fnorm1 / dble(isum))
      ! fnorm2 = dsqrt( fnorm2 / dble(isum2) )
      ! write(101,*)iter,fnorm1,fnorm2
      ! if (fnorm2 < 5d-4 .and. fnorm1 < 5d-3 ) then
      !   write(logfile,'("* |Forces| < NEBconv :: CONVERGED")')
      !   exit
      ! endif

    enddo


    ! Output the RMS displacement from the initial path, ignoring fixed end-points.
    !
    call NewPath(rp_init, .TRUE., startfile, endfile, 'initial_reoriented_path.xyz', nimage, pathinit, &
      .FALSE., rp%cx(1)%na)
    distance = 0.d0
    do i = 2, rp%nimage-1
      do j = 1, rp%cx(i)%na
        dx = rp%cx(i)%r(1, j) - rp_init%cx(i)%r(1, j)
        dy = rp%cx(i)%r(2, j) - rp_init%cx(i)%r(2, j)
        dz = rp%cx(i)%r(3, j) - rp_init%cx(i)%r(3, j)
        distance = distance + dsqrt(dx*dx + dy*dy + dz*dz)
      enddo
    enddo
    distance = distance * bohr_to_ang
    write(logfile, '("* Sum of displacements from initial path =",1x,f14.8,1x,"Angstrom")') distance

    ! Output the energy along the final path.
    !
    open(15, file = trim( inputfile )//'.energy-neb-end', status = 'unknown')
    write(15, '("# NEB energy profile at END of (CI)NEB run (in au)")')
    write(15, '("# Lambda    |    Energy / Eh  | Relative energy / Eh")')
    lambda = 0.0
    imax = 1
    do i = 1, rp%nimage
      write(15, '(3(1x, f14.8))') lambda, rp%cx(i)%vcalc, (rp%cx(i)%vcalc - rp%cx(1)%vcalc)
      if (i .lt. rp%nimage) lambda = lambda + norm2(reshape(rp%cx(i+1)%r-rp%cx(i)%r, (/rp%na*3/)))
      if (rp%cx(i)%vcalc > rp%cx(imax)%vcalc) imax = i
    enddo
    close(15)

    ! Print the TS structure.
    call PrintCXSToFile(rp%cx(imax), trim(inputfile)//'_opt_ts.xyz', 0.0d0)

    ! Print the final path to inputfile_nebfinal.xyz
    !
    call PrintPathToFile(rp, trim( inputfile ) //'_nebfinal.xyz')

    ! Close up files.
    !
    close(15)
    close(16)
    close(18)

    return
  end Subroutine CINEB


end Module pathopt
