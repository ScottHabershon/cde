!
!***************************************************************************************
!
!> @brief Contains routines for automatic determination of a reaction-mechanism
!! connecting two distinct (user-input) structures.
!!
!! The method employed here has been reported in:
!!
!! Automatic Proposal of Multistep Reaction Mechanisms using a Graph-Driven Search,
!! I. Ismail and H. B. V. A. Stuttaford-Fowler and C. Ochan Ashok and C. Robertson and S. Habershon,
!! J Phys Chem A  123  3407-3417  (2019)
!!
!! Fast screening of homogeneous catalysis mechanisms using graph-driven searches and
!! approximate quantum chemistry, C Robertson, S Habershon, Catalysis Science & Technology,
!! 9, 6357 (2019)
!!
!! Bibtex records:
!!
!! @article{Habershon:2019ab,
!!	Author = {Ismail, Idil and Stuttaford-Fowler, Holly B V A and Ochan Ashok, Curtis and Robertson, Christopher and Habershon, Scott},
!!	Journal = {J Phys Chem A},
!! Month = {Apr},
!! Number = {15},
!! Pages = {3407-3417},
!! Title = {Automatic Proposal of Multistep Reaction Mechanisms using a Graph-Driven Search},
!! Volume = {123},
!! Year = {2019}}
!!
!! @article{Habershon:2019ee,
!!	Author = {Robertson, Christopher and Habershon, Scott},
!!	Date-Added = {2019-12-03 15:41:20 +0000},
!!	Date-Modified = {2019-12-03 15:41:29 +0000},
!!	Doi = {10.1039/C9CY01997A},
!!	Issue = {22},
!!	Journal = {Catal. Sci. Technol.},
!!	Pages = {6357-6369},
!!	Publisher = {The Royal Society of Chemistry},
!!	Title = {Fast screening of homogeneous catalysis mechanisms using graph-driven searches and approximate quantum chemistry},
!!	Url = {http://dx.doi.org/10.1039/C9CY01997A},
!!	Volume = {9},
!!	Year = {2019},
!!  }
!!
!!---------------------------
!! DESCRIPTION
!!---------------------------
!!
!! Given input xyz files for target reactants and products, this routine will
!! try to find a reaction mechanism connecting the two by using graph-moves defined
!! in the "movefile".
!!
!! After determination of a path which definitively connects reactants and products,
!! we then generate xyz-files for each reaction-path along the reaction-mechanism
!! using the idea of optimization under the graph-restraining potential.
!!
!! If optaftermove = .TRUE. in the input file, this code will also perform geometry
!! optimization for all of the intermediate structures and report the results in the
!! .log file.
!!
!!--------------------------
!! OUTPUTS
!!--------------------------
!!
!! After completion, the following outputs are worth looking at:
!! - mcopt.dat: contains the graph-error function as a function of number of MC iterations
!!              in the calculation. This should reach zero for a successful calculation.
!!
!! - final_path.xyz: Contains the intermediate structures (xyz format) along the located reaction
!!                   mechanism. If optaftermove = .true., the structures are geometry-optimized using
!!                   whichever electronic structure methods is specified in the input file.
!!
!! - final_path_rx_*.xyz: This series of xyz files, one for each reaction-step in the mechanism,
!!                        contain initial nimage-long paths for each intermediate reaction.
!!
!! - .log: The .log file will contain useful information about the reaction, including the
!!         calculated Q-value and the maximum dE value.
!!
!! In addition, if igfunc = 4 (so we are targetting formation of a single molecule or
!! subset of molecules, rather than an entire structure), the further outputs are
!! produced which have the prefix 'adjusted_'. This are files which have been determined
!! for the reaction path after removing molecules and reactions which are not relevant
!! to the formation of the product.
!
!***************************************************************************************
!

Module pathfinder

  use globaldata
  use constants
  use pes
  use chemstr
  use rpath
  use io
  implicit none

contains

  !************************************************************************
  !> RunPathFinder
  !!
  !! Runs automatic determination of a reaction-path connecting
  !! two structures.
  !!
  !************************************************************************
  !
  Subroutine RunPathFinder( )
    implicit none
    real(8) :: error, error_old, rx, vbest, ir
    real(8) :: x,y,z, gsum, Qvalue, Qvalue2,dE, Qmax,error_temp,bsum
    integer :: ii, imove, irx, itype,k,it,ncount,imcount,nsuccess,ia,nn,ncharge,TotalCharge
    integer :: i, j, imc, ix(NAMOVEMAX), na,isum,isum2,i1,i2,im1(NRXNMAX),im2(NRXNMAX)
    integer :: iend
    integer, allocatable ::  movenum(:), moveatoms(:,:)
    integer, allocatable ::  msmovenum(:,:), msmoveatoms(:,:,:)
    integer, allocatable :: chargemove(:,:), chargemove_store(:,:)
    type (cxs) :: cx_start, cx_end
    type(cxs), allocatable :: cx(:), cxtemp(:)
    logical :: errflag, fail, success
    logical :: printflag,cyc, errx,ChangeCharges
    logical, allocatable :: atomchange(:), bondchange(:,:)
    integer :: rxindex(NAMAX), nrx
    integer, allocatable :: movenum_store(:), moveatoms_store(:,:)
    integer :: itry, idof, irxn
    real(8) :: beta, fac, mctstore, vbe, TotalError, grapherror,maxbarrier
    type(rxp) :: rp
    character (len=7) :: x1
    character (len=8) :: fmt
    character (len=20) :: printfile, file_root

    ! Open output file mcopt.dat to print graph-error function.
    !
    write(logfile,'("* Running pathfinder calculation...."/)')
    write(logfile,'("* Optimization function information will be printed to mcopt.dat")')
    open(72,file='mcopt.dat',status='unknown')
    write(72,'("#    Step number     |     TOTAL error    |    GRAPH error    |    VBE term ")')

    ! Assign the end-point structures from the input files, startfile and endfile.
    ! Note that startfile and endfile are read from the main input file.
    !
    write(logfile,'(/"* Reading reactant structure...")')
    Call ReadCXS( cx_start, startfile )
    Call SetMass(cx_start)
    Call GetGraph( cx_start )
    Call Getmols(cx_start)
    Call PrintCXSGraphInfo(cx_start,logfile,"Reactant structure")

    write(logfile,'(/"* Reading product structure...")')
    Call ReadCXS( cx_end, endfile )
    Call SetMass(cx_end)
    Call GetGraph( cx_end )
    Call Getmols(cx_end)
    Call PrintCXSGraphInfo(cx_end,logfile,"Product structure")


    ! Allocate space for atomchange and bondchange arrays - these will indicate at each MC search
    ! step which atoms and bonds are allowed to change.
    !
    na = max(cx_start%na,cx_end%na)
    allocate( atomchange(na) )
    allocate( bondchange(na,na) )


    ! Read the graphmoves.
    !
    Call ReadGraphMoves( movefile )


    ! Based on the movefile, decide whether or not we're also going to have to
    ! consider changes in charge states too.
    !
    ChangeCharges = .FALSE.
    do i = 1, ngmove
      if (namove(i) == 0) then
        ChangeCharges = .TRUE.
      endif
    enddo


    ! Create space for the TOTAL reaction-string. This is defined by
    ! nrxn in the input file. Here, nrxn defines the (maximum) number of reactions
    ! which can occur. Here, cx(1) is the structure formed after application of reaction
    ! 1 to the cx_start, cx(2) is the structure formed after application of reaction 2 to
    ! cx(1), and so on....
    !
    allocate( cx(nrxn) )


    ! Copy all aspects of the end-point structures to the intermediates for consistency:
    !
    do i = 1, nrxn
      Call CopytoNewCXS( cx_start, cx(i) )
      Call SetMass(cx(i))
    enddo

    ! Allocate space for the moves:
    !
    allocate( movenum(nrxn) )
    allocate( moveatoms(nrxn,NAMOVEMAX) )
    allocate( movenum_store(nrxn) )
    allocate( moveatoms_store(nrxn,NAMOVEMAX) )
    allocate( chargemove(nrxn,NMOLMAX) )
    allocate( chargemove_store(nrxn,NMOLMAX))

    ! Set all the initial moves to be null moves.
    !
    write(logfile,'("* Setting all initial moves to null...")')
    do i = 1, nrxn
      movenum(i) = 0
      moveatoms(i,1:NAMOVEMAX) = 0
    enddo


    ! Temp - test for Noha reaction
!    movenum(1) = 5
!    moveatoms(1,1) = 68
!    moveatoms(1,2) = 1
!    moveatoms(1,3) = 18
!    movenum(2) = 5
!    moveatoms(2,1) = 21
!    moveatoms(2,2) = 2
!    moveatoms(2,3) = 20


    ! Zero the charges for all molecules along the starting path...
    !
    cx_start%molcharge(:) = 0
    do i = 1, nrxn
      cx(i)%molcharge(:) = 0
      chargemove(i,:) = 0
    enddo

    ! Evaluate the initial error, defined as the graph-distance between the
    ! graph resulting from applying nrxn reactions and the graph for the
    ! end-point:
    !
    Call GetPathFitness( cx_start, cx_end, cx, nrxn, movenum, moveatoms, errflag, &
    GraphError, TotalError, vbe )

!     print*,'ERROR = ',GraphError

!  Call GraphsToCoords(cx_start, cx, nrxn, .TRUE., 'final_path.xyz' )
! stop

    ! If an error flag is returned here, something is wrong with the initial path....
    !
    if (errflag)Stop '* Initial errflag for path is TRUE - something weird going &
    on with initial path....go check it out....'


    ! Output initial errors to logfile.
    !
    write(logfile,'(/"* Initial Graph Error = ",2x,f10.5/)')GraphError
    write(logfile,'("* Initial Total Error = ",2x,f10.5/)')TotalError
    imc = 0
    write(72,41)imc,totalerror,grapherror,vbe,dble(imc)
    call flush(logfile)


    ! Loop over MC moves, trying to anneal down to zero error:
    !
    write(logfile,'("*** Starting reaction-path search ***"/)'); call flush(logfile)
    mctstore = mcrxntemp
    nsuccess = 0
    outer: do imc = 1, nmcrxn

      if (grapherror <= GCONV .and. imc == 1) then
        write(72,41)imc,totalerror,grapherror,vbe,dble(nsuccess)/dble(imc)
        call flush(72)
        write(logfile,'("================================================================")')
        write(logfile,'("            ZERO_ERROR MECHANISM LOCATED AT FIRST STEP...       ")')
        write(logfile,'("================================================================"/)')
        call flush(logfile)
        exit outer
      endif

      ! Linear decrease of temperature...
      !
      mcrxntemp = (1.d0 - (dble(imc)/dble(nmcrxn)) ) * mctstore

      ! Store old error value.
      !
      error_old = totalerror

      ! Make a change to the current reaction-mechanism...
      !
      Call UpdateMechanism(nrxn,movenum,moveatoms,bondchange,atomchange,na,cx_start,cx,rxindex,cyc,&
      movenum_store, moveatoms_store)
      if (cyc) cycle outer


      ! Propagate the graph and evaluate the new error
      !
      Call GetPathFitness( cx_start, cx_end, cx, nrxn, movenum, moveatoms, errflag, &
      GraphError, TotalError, vbe )


      ! Check for convergence, based on GRAPHERROR!!!
      !
      if (grapherror <= GCONV) then
        write(72,41)imc,totalerror,grapherror,GCONV,dble(nsuccess)/dble(imc)
        call flush(72)
        write(logfile,'("================================================================")')
        write(logfile,'("            SUCCESSFULLY LOCATED ZERO_ERROR MECHANISM           ")')
        write(logfile,'("================================================================"/)')
        call flush(logfile)
        exit outer
      endif

      ! Do we change the charges on the molecules along the path?
      !
      if (ChangeCharges) Call UpdateCharges(nrxn,cx,chargemove, chargemove_store, errflag)


      ! Metropolis reject/accept criterion.
      !
      if (errflag) then
        movenum(:) = movenum_store(:)
        moveatoms(:,:) = moveatoms_store(:,:)
        Totalerror = error_old
        if (changecharges) then
          chargemove = chargemove_store
        endif
      else
        beta = 1.d0 / (KBOLTZ * mcrxntemp)
        fac = exp(- beta * (TotalError - error_old) )
        Call random_number(ir)
        if ( ir > fac ) then
          movenum(:) = movenum_store(:)
          moveatoms(:,:) = moveatoms_store(:,:)
          Totalerror = error_old
          if (changecharges) then
            chargemove = chargemove_store
          endif
        else
          nsuccess = nsuccess + 1
        endif
      endif

      ! Output current error every 1000 steps.
      !
      if (mod(imc,1000) == 0) then
        write(72,41)imc,totalerror,grapherror,GCONV,dble(nsuccess)/dble(imc);call flush(72)
      endif

    enddo outer
    close(72)

    ! Print a warning and stop if the calculation ran but didn't converge.....
    !
    if (GraphError > GCONV) then
      write(logfile,'("NOPE! ==========================================================")')
      write(logfile,'("                FAILED TO FIND A ZERO_ERROR MECHANISM           ")')
      write(logfile,'("               xyz output files WILL NOT be printed out!        ")')
      write(logfile,'("NOPE! =========================================================="/)')
      call flush(logfile)
      stop
    endif

    ! After finishing optimization, propagate the final string, evaluate the error
    ! and output details of final string to logfile.
    !
    Call GetPathFitness( cx_start, cx_end, cx, nrxn, movenum, moveatoms, errflag, &
    GraphError, TotalError, vbe )

    ! At this point, we find the first instance of the target molecule(s)
    ! being formed then set the remainder of movenum to be zero....
    !
    Call TrimPath(cx_start, cx_end, cx, nrxn, movenum, moveatoms, iend )

    ! Print molecules generated along reaction:
    !
    Call PrintMolsAlongPath( nrxn, cx_start, cx_end, cx, movenum,chargemove, &
    changecharges,moveatoms )


    ! At this point, we need to turn graphs into coordinates....
    !
    write(logfile,'("* Generating structures from connectivity matrices..."/)')
    write(logfile,'("* Reaction end-points printed to final_path.xyz"/)')
    Call GraphsToCoords(cx_start, cx, nrxn, .TRUE., 'final_path.xyz' )
    write(logfile,'("* Finished structure generation."/)')


    ! Now we need to do error correction - if the end-point structures have changed due to geometry
    ! optimization, such that the structures no longer represent the zero-error reaction path, we need to
    ! work out which reactions "did" happen, then modify the movenum(:) and moveatoms(:,:). We also need to
    ! check that the reaction-string is still valid.
    !
    ! write(6,*)'OPT = ',optaftermove
    ! if (optaftermove) then
    !    Call ReactionErrorCorrection(cx_start, cx, nrxn, movenum, moveatoms,errx)
    !    if (errx) then
    !      write(logfile,'(/"*** ERROR CORRECTION ATTEMPT FAILED ***"/)')
    !    else
    !      write(logfile,'(/"*** ERROR CORRECTION ATTEMPT SUCCESSFUL ***"/)')
    !      write(logfile,'("\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/")')
    !      write(logfile,'("*** CORRECTED MECHANISM ***")')
    !      Call PrintMolsAlongPath( nrxn, cx_start, cx_end, cx, movenum,chargemove,changecharges,moveatoms )
    !
    !     ! Recalcualte energy.
    !      write(6,*)'INTO AbInitio here...',cx_start%nmol;call flush(6)
    !     Call AbInitio( cx_start, 'ener', success )
    !     do i = 1, nrxn
    !       print*,'INTO AbInitio here...',cx(i)%nmol,cx(i)%molcharge(1:cx(i)%nmol)
    !       Call AbInitio( cx(i), 'ener', success)
    !     enddo
    !
    !   endif
    ! endif

    ! Calculate and output the energy information, including the calculated Q-value:
    !
    if (optaftermove) then
      open(94,file = 'final_path_energy.dat', status='unknown')
      write(94,'("# Reaction index | Absolute energy in au | Relative energy of reaction in kJ/mol")')
      write(94,'("# Energy relative to reactants in kJ/mol")')
      write(logfile,'("===============        ENERGY INFORMATION         ==============="/)')
      write(logfile,'("* Initial energy =",3x,f14.8,2x,"Eh")')cx_start%vcalc
      call flush(logfile)
      write(94,*)'0  ',cx_start%vcalc, '  0.00 ','   0.00 '
      Qvalue = 0.d0
      Qvalue2 = 0.d0
      Qmax = -1d6
      ncount = 0
      do irxn = 1, nrxn
        if (irxn == 1) then
          dE = (cx(irxn)%vcalc - cx_start%vcalc) * au_to_kjmol
        else
          dE = (cx(irxn)%vcalc - cx(irxn-1)%vcalc) * au_to_kjmol
        endif
        write(94,*)irxn,cx(irxn)%vcalc, dE,(cx(irxn)%vcalc-cx_start%vcalc)*au_to_kjmol
        if (abs(de) > Qmax)Qmax = abs(de)
        Qvalue2 = Qvalue2 + abs(dE)
        if (movenum(irxn) /= 0) then
          ncount = ncount + 1
          Qvalue = Qvalue + abs(dE)
        endif
        write(logfile,'("* Reaction #",1x,i4,1x,":",2x,"dE = ",2x,f14.8)')irxn,dE
      enddo
      Qvalue = Qvalue / dble(ncount)
      write(logfile,'(/"* Final Q-value (NOT averaged over non-zero rxns) =",2x,f16.8,2x,"kJ/mol"/)')Qvalue2
      write(logfile,'(/"* Number of non-zero reactions =",2x,i4,2x,"kJ/mol"/)')ncount
      write(logfile,'(/"* Final Q-value (averaged over non-zero rxns) =",2x,f16.8,2x,"kJ/mol"/)')Qvalue
      write(logfile,'(/"* Maximum dE =",2x,f16.8,2x,"kJ/mol"/)')Qmax
      write(logfile,'("================================================================"/)')
      close(94)
    endif

    ! For each reaction along the path, print nimage intermediate snapshots to file.
    ! These snapshots can be used for subsequenct NEB calculations.
    !
    write(logfile,'("* Reaction paths printed to final_path_rx_*.xyz")'); call flush(logfile)
    file_root = 'final_path'
    Call PrintMechanismPaths(nrxn, cx_start, cx, file_root, maxbarrier,bsum,movenum )
    write(logfile,*)"* Maximum barrier = ",maxbarrier * au_to_kjmol,"kJ/mol";call flush(logfile)
    write(logfile,*)"* Barrier sum = ",bsum,"au";call flush(logfile)


    ! If we're running a calculation with igfunc=4, so we're trying to target formation
    ! of a specific molecule (or set of molecules), then we need to only focus on
    ! outputting mechanism and energy data for the set of molecules leading to the target.
    !
    if (igfunc == 4)Call AdjustPaths(cx_start, cx_end, cx, nrxn, movenum, moveatoms,iend )


    ! Clean up memory.
    !
    deallocate( movenum )
    deallocate( moveatoms )
    deallocate( movenum_store )
    deallocate( moveatoms_store )

    41 Format(i7,3x,f16.7,3x,f16.7,3x,f16.7,3x,f16.7)

    return
  end Subroutine RunPathFinder


  !************************************************************************
  !> AdjustPaths
  !!
  !! This routine performs a series of functions to
  !! simplify the post-analysis of reaction-path simulations.
  !! (i) If pesfull = .FALSE., we print out the energies and Q-values
  !!     calculated only for those molecules involved in formation of the
  !!     product.
  !! (ii) We print out a new set of reaction paths which only involve
  !!      the molecules containing atoms ending up in the final structure.
  !! (iii) We strip the adjusted paths of non-reactive molecules.
  !!
  !! The final set of xyz path files should then be most appropriate
  !! for post-analysis.
  !!
  !! Input parameters:
  !! - cx_start: CXS object for reactants.
  !! - cx_end: CXS object for target product.
  !! - cx(nrxn): Intermediate CXS objects along mechanism path.
  !! - nrxn: Maximum number of active reactions.
  !! - movenum(:): Graph-move number at each reaction step.
  !! - moveatoms(:,:): Atoms involved in graph-moves at each step.
  !! - ifound: Reaction step-number at which target products are first
  !!           found.
  !!
  !************************************************************************
  Subroutine AdjustPaths(cx_start, cx_end, cx, nrxn, movenum, moveatoms, ifound )
    implicit none
    integer :: i,j,k,nrxn, ifound,itarg,irxn,ncount,vtot_last,iend,jj
    type(cxs) :: cx_start, cx_end, cx(nrxn),cx1,cx2
    integer :: movenum(nrxn), moveatoms(nrxn,NAMOVEMAX),im,imove
    integer :: m, ii, n1
    real(8) :: grapherror,Qvalue,Qvalue2,dE,vtot,vtot_start,Qmax,v1,v2
    real(8) :: drsq1, drsq2, xtemp, ytemp, ztemp,dx,dy,dz,de_sum
    logical :: errflag, active(0:NRXNMAX,NMOLMAX),barrierless
    character (len=2) :: ctemp
    character (len=7) :: x1
    character (len=8) :: fmt
    character (len=20) :: file_root

    write(logfile,'(/"================================================================"/)')
    write(logfile,'(/"* Adjusting paths to focus on target formation..."/)')

    ! Propagate the graph from the start to generate the final graph:
    !
    !  Call PropagateGraphs( cx_start, cx, nrxn, movenum, moveatoms, errflag, GraphError )

    ! Determine which reaction-step has first appearance of products.
    !
    write(logfile,'("* Target found in reaction:",2x,i4)')ifound
    write(logfile,'("* Molecule ID of target is:",2x,i4)')cx(ifound)%itargetmol
    call flush(logfile)
    itarg = cx(ifound)%itargetmol


    ! We now need to step backwards from the frame in which the target molecules
    ! were found and keep track of all the molecules containing atoms..
    !
    ! First, set the active flag to true for the product molecule (itarg) in the frame
    ! when it is found (ifound).
    !
    active(:,:) = .FALSE.
    active(ifound,itarg) = .true.

    ! Now for the ifound-1 frame, work out the active molecules - this is neccesary,
    ! so when can then recheck the ifound frame to account for active molecules
    ! formed by e.g. dissociation, which we might have missed if we only set the
    ! target molecule flag to "active".
    !
    do jj = 1, nrxn

      irxn = ifound - 1
      do i = 1, cx(irxn)%nmol
        do j = 1, cx(irxn+1)%nmol
          if (active(irxn+1,j)) then
            if (ContainsTargetAtom(cx(irxn),i,cx(irxn+1),j)) then
              active(irxn,i) = .TRUE.
            endif
          endif
        enddo
      enddo

      ! Now, based on the active molecules in ifound-1, make sure that everything
      ! is activated in frame ifound.
      !
      do i = 1, cx(ifound)%nmol
        do j = 1, cx(ifound-1)%nmol
          if (active(ifound-1,j)) then
            if (ContainsTargetAtom(cx(ifound),i,cx(ifound-1),j)) then
              active(ifound,i) = .TRUE.
            endif
          endif
        enddo
      enddo


      ! Now that we've identified molecules which are active in the ifound
      ! frame, go back and identify active molecules in other frames.
      !
      do irxn = ifound-1, 0, -1

        ! Run through the molecules in this frame and check whether they
        ! contain any of the active molecules in the previous frame.
        !
        if (irxn > 0 ) then

          do i = 1, cx(irxn)%nmol
            do j = 1, cx(irxn+1)%nmol
              if (active(irxn+1,j)) then
                if (ContainsTargetAtom(cx(irxn),i,cx(irxn+1),j)) then
                  active(irxn,i) = .TRUE.
                endif
              endif
            enddo
          enddo

        else if (irxn == 0) then

          do i = 1, cx_start%nmol
            do j = 1, cx(irxn+1)%nmol
              if (active(irxn+1,j))then
                if (ContainsTargetAtom(cx_start,i,cx(irxn+1),j)) then
                  active(irxn,i) = .TRUE.
                endif
              endif
            enddo
          enddo

        endif
      enddo

      ! Now go forwards, activating molecules.
      !
      do irxn = 1, ifound

        ! Run through the molecules in this frame and check whether they
        ! contain any of the active molecules in the previous frame.
        !
        if (irxn == 1 ) then

          do i = 1, cx(irxn)%nmol
            do j = 1, cx_start%nmol
              if (active(irxn-1,j)) then
                if (ContainsTargetAtom(cx(irxn),i,cx_start,j)) then
                  active(irxn,i) = .TRUE.
                endif
              endif
            enddo
          enddo

        else if (irxn > 1) then

          do i = 1, cx(irxn)%nmol
            do j = 1, cx(irxn-1)%nmol
              if (active(irxn-1,j))then
                if (ContainsTargetAtom(cx(irxn),i,cx(irxn-1),j)) then
                  active(irxn,i) = .TRUE.
                endif
              endif
            enddo
          enddo

        endif
      enddo

    enddo

    ! Q - do we have to keep going through above?

    ! Now calculate the energy along the reaction path using only active molecules
    !
    if (optaftermove .eqv. .TRUE. .and. pesfull .eqv. .FALSE. ) then
      open(94,file = 'adjusted_path_energy.dat', status='unknown')
      write(94,'("# Reaction index | Absolute energy in au | Relative energy of reaction in kJ/mol")')
      write(94,'("# Energy relative to reactants in kJ/mol")')
      write(logfile,'("=============== ADJUSTED ENERGY INFORMATION ===============")')

      ! Start at reactants and, for each reaction, calculate the energy
      ! change for the active molecules.
      !
      ncount = 0
      Qvalue = 0.d0
      Qvalue2 = 0.d0
      Qmax = -1d6
      de_sum = 0.d0
      do i = 0, ifound-1
        im = i + 1

        write(logfile,'(/"- REACTION NUMBER:",2x,i4/)')im

        v1 = 0.d0
        v2 = 0.d0

        if (i == 0) then

          outer2: do j = 1, cx_start%nmol
            if (active(i,j)) then

              ! Check if this molecule contains any of the reactive atoms.
              !
              imove = movenum(im)
              do k = 1, namove(imove)
                do m = 1, cx_start%namol(j)
                  if (cx_start%molid(j,m) == moveatoms(im,k))then
                    write(logfile,*)"- REACTANT:",j,cx_start%molen(j),MolecularFormula( cx_start, j )
                    v1 = v1 + cx_start%molen(j)
                    cycle outer2
                  endif
                enddo
              enddo
            endif
          enddo outer2

        else if (i /= 0) then

          outer3: do j = 1, cx(i)%nmol
            if (active(i,j)) then

              ! Check if this molecule contains any of the reactive atoms.
              !
              imove = movenum(im)
              do k = 1, namove(imove)
                do m = 1, cx(i)%namol(j)
                  if (cx(i)%molid(j,m) == moveatoms(im,k))then
                    write(logfile,*)"- REACTANT:",j,cx(i)%molen(j),MolecularFormula(cx(i),j)
                    v1 = v1 + cx(i)%molen(j)
                    cycle outer3
                  endif
                enddo
              enddo
            endif
          enddo outer3

        endif

        ! Now calculate the energy after reaction.
        outer4: do j = 1, cx(i+1)%nmol
          if (active(i+1,j)) then

            ! Check if this molecule contains any of the reactive atoms.
            !
            imove = movenum(im)
            do k = 1, namove(imove)
              do m = 1, cx(i+1)%namol(j)
                if (cx(i+1)%molid(j,m) == moveatoms(im,k))then
                  write(logfile,*)"+ PRODUCT:",j,cx(i+1)%molen(j),MolecularFormula( cx(i+1), j )
                  v2 = v2 + cx(i+1)%molen(j)
                  cycle outer4
                endif
              enddo
            enddo
          endif
        enddo outer4

        dE = (v2 - v1) * au_to_kjmol
        de_sum = de_sum + de
        if (abs(de) > Qmax)Qmax = abs(de)
        Qvalue2 = Qvalue2 + abs(dE)
        if (movenum(im) /= 0) then
          ncount = ncount + 1
          Qvalue = Qvalue + abs(dE)
        endif

        write(logfile,'("* Reaction Energy change (au,au,kJ/mol):",1x,i3,3(f14.8))')im,v1,v2,de
        write(94,*)im, de

        ! Store the reactant energy...
        !
        if (de_sum > 0.d0) then
          write(logfile,'("*** BARRIER REACTION DETECTED ***")')
        endif

      enddo
      Qvalue = Qvalue / dble(ncount)
      write(logfile,'(/"* ADJUSTED Final Q-value (NOT averaged over non-zero rxns) =",2x,f16.8,2x,"kJ/mol"/)')Qvalue2
      write(logfile,'(/"* ADJUSTED Number of non-zero reactions =",2x,i4,2x,"kJ/mol"/)')ncount
      write(logfile,'(/"* ADJUSTED Final Q-value (averaged over non-zero rxns) =",2x,f16.8,2x,"kJ/mol"/)')Qvalue
      write(logfile,'(/"* ADJUSTED Maximum dE =",2x,f16.8,2x,"kJ/mol"/)')Qmax
      write(logfile,'(/"* Is it barrierless?",2x,L/)')barrierless

      write(logfile,'("================================================================"/)')
      close(94)
    endif

    ! Now output the structures along the reaction-path, but only using the
    ! active molecules.
    !
    ! NOTE: This might mean that different paths have different numbers of
    ! atoms or molecules! Might be tought to visualize in e.g. VMD....
    !
    ! NOTE: We can't guarantee that the atom list in each molecule will be the same
    !       order! This will have to be sorted out before we perform e.g. NEB!
    !
    file_root = 'adjusted_path'
    do i = 0, ifound-1
      im = i + 1

      if (im <= 9) then
        fmt = '(I1.1)'
        write (x1,fmt)im
      else
        fmt = '(I2.2)'
        write (x1,fmt)im
      endif

      ! Print out reactant structure for this reaction.
      !
      open(21,file=trim(file_root)//'_reactants_'//trim(x1)//'.xyz')

      ! First, count number of atoms in active molecules...
      !
      n1 = 0
      if (i == 0) then
        do j = 1, cx_start%nmol
          if (active(0,j)) n1 = n1 + cx_start%namol(j)
        enddo
      else
        do j = 1, cx(i)%nmol
          if (active(i,j)) n1 = n1 + cx(i)%namol(j)
        enddo
      endif

      ! Write correct number of atoms to file.
      !
      write(21,*)n1
      write(21,*)

      ! Now print out atomic coordinates.
      !
      if (i == 0) then
        do j = 1, cx_start%nmol
          if (active(0,j)) then
            do k = 1, cx_start%namol(j)
              ii = cx_start%molid(j,k)
              write(21,*)cx_start%atomlabel(ii),cx_start%r(1:3,ii)*bohr_to_ang
            enddo
          endif
        enddo
      else
        do j = 1, cx(i)%nmol
          if (active(i,j)) then
            do k = 1, cx(i)%namol(j)
              ii = cx(i)%molid(j,k)
              write(21,*)cx(i)%atomlabel(ii),cx(i)%r(1:3,ii)*bohr_to_ang
            enddo
          endif
        enddo
      endif
      close(21)

      ! Now we've printed reactants, let's print out the products too.
      !
      open(21,file=trim(file_root)//'_products_'//trim(x1)//'.xyz')

      ! Again, work out the number of atoms.
      !
      n1 = 0
      if (i == nrxn) then
        do j = 1, cx_end%nmol
          if (active(i,j)) n1 = n1 + cx_end%namol(j)
        enddo
      else
        do j = 1, cx(i+1)%nmol
          if (active(i+1,j)) n1 = n1 + cx(i+1)%namol(j)
        enddo
      endif

      ! Write out number of atoms.
      !
      write(21,*)n1
      write(21,*)

      ! Write out the product coordinates.
      !
      if (i == nrxn) then
        do j = 1, cx_end%nmol
          if (active(i,j)) then
            do k = 1, cx_end%namol(j)
              ii = cx_end%molid(j,k)
              write(21,*)cx_end%atomlabel(ii),cx_end%r(1:3,ii)*bohr_to_ang
            enddo
          endif
        enddo
      else
        do j = 1, cx(i+1)%nmol
          if (active(i+1,j)) then
            do k = 1, cx(i+1)%namol(j)
              ii = cx(i+1)%molid(j,k)
              write(21,*)cx(i+1)%atomlabel(ii),cx(i+1)%r(1:3,ii)*bohr_to_ang
            enddo
          endif
        enddo
      endif
      close(21)
    enddo

    ! Finally, for each pair of reactants and products, we reorder the atoms in the product
    ! structure to best match those in the reactants.
    !
    file_root = 'adjusted_path'
    do i = 0, ifound-1
        im = i + 1

        if (im <= 9) then
          fmt = '(I1.1)'
          write (x1,fmt)im
        else
          fmt = '(I2.2)'
          write (x1,fmt)im
        endif

        Call ReadCXS( cx1, trim(file_root)//'_reactants_'//trim(x1)//'.xyz' )
        Call ReadCXS( cx2, trim(file_root)//'_products_'//trim(x1)//'.xyz' )

        if (cx1%na /= cx2%na)stop 'ERROR: inconsistent numbers of atoms in path-reordering &
        in AdjustPaths...'

        ! Reorder atoms....first use bubble sort to put atoms in atomic number...
        !
        do k = 1, cx1%na ! Loop over run-throughs.
          do j = 1, cx1%na-1 ! Loop over atoms in cx1
            do m = j+1, cx1%na
              if ( LabelToNumber(cx1%atomlabel(j)) <= LabelToNumber(cx1%atomlabel(m) ) ) then
                ctemp = cx1%atomlabel(j)
                xtemp = cx1%r(1,j)
                ytemp = cx1%r(2,j)
                ztemp = cx1%r(3,j)
                cx1%atomlabel(j) = cx1%atomlabel(m)
                cx1%r(1:3,j) = cx1%r(1:3,m)
                cx1%atomlabel(m) = ctemp
                cx1%r(1,m) = xtemp
                cx1%r(2,m) = ytemp
                cx1%r(3,m) = ztemp
              endif
            enddo
          enddo
        enddo

        ! Do the same for the products.
        !
        do k = 1, cx2%na ! Loop over run-throughs.
          do j = 1, cx2%na-1 ! Loop over atoms in cx1
            do m = j+1, cx2%na
              if ( LabelToNumber(cx2%atomlabel(j)) <= LabelToNumber(cx2%atomlabel(m) ) ) then
                ctemp = cx2%atomlabel(j)
                xtemp = cx2%r(1,j)
                ytemp = cx2%r(2,j)
                ztemp = cx2%r(3,j)
                cx2%atomlabel(j) = cx2%atomlabel(m)
                cx2%r(1:3,j) = cx2%r(1:3,m)
                cx2%atomlabel(m) = ctemp
                cx2%r(1,m) = xtemp
                cx2%r(2,m) = ytemp
                cx2%r(3,m) = ztemp
              endif
            enddo
          enddo
        enddo

        ! Now loop through the products and permute atoms if they are "closer" to the
        ! reactant structure...
        do m = 1, cx1%na
        do k = 1, cx1%na
          do j = 1, cx2%na ! Loop over atoms in cx1

            if (cx1%atomlabel(k) == cx2%atomlabel(j)) then

              dx = cx2%r(1,j) - cx1%r(1,k)
              dy = cx2%r(2,j) - cx1%r(2,k)
              dz = cx2%r(3,j) - cx1%r(3,k)
              drsq1 = dx*dx + dy*dy + dz*dz
              dx = cx2%r(1,j) - cx1%r(1,j)
              dy = cx2%r(2,j) - cx1%r(2,j)
              dz = cx2%r(3,j) - cx1%r(3,j)
              drsq2 = dx*dx + dy*dy + dz*dz

              ! if drsq1 < drsq2, j is closer to k, so swap!
              if ( drsq1 < drsq2 ) then
                ctemp = cx2%atomlabel(j)
                xtemp = cx2%r(1,j)
                ytemp = cx2%r(2,j)
                ztemp = cx2%r(3,j)
                cx2%atomlabel(j) = cx2%atomlabel(k)
                cx2%r(1:3,j) = cx2%r(1:3,k)
                cx2%atomlabel(k) = ctemp
                cx2%r(1,k) = xtemp
                cx2%r(2,k) = ytemp
                cx2%r(3,k) = ztemp
              endif
            endif
            enddo
          enddo
        enddo

        ! Finally, print the two reordered structures to file.
        !
        Call PrintCXSToFile(cx1,'reorder_reactants_'//trim(x1)//'.xyz',0.d0)
        Call PrintCXSToFile(cx2,'reorder_products_'//trim(x1)//'.xyz',0.d0)

        CAll DeleteCXS(cx1)
        Call DeleteCXS(cx2)

      enddo

    return
  end Subroutine AdjustPaths


  !************************************************************************
  !> ContainsTargetAtom
  !!
  !! Logical function to check whether cx1 contains any of the atoms in
  !! molecule number "itarget" in cx2.
  !!
  !************************************************************************
  !
  Function ContainsTargetAtom(cx1,mol1,cx2,itarget) result(flag)
    implicit none
    type(cxs) :: cx1,cx2
    integer :: mol1,j,k,itarget
    logical :: flag

    flag = .FALSE.
    outer1: do k = 1, cx1%namol(mol1)
      do j = 1, cx2%namol(itarget)
        if ( cx1%molid(mol1,k) == cx2%molid(itarget,j)) then
          flag = .TRUE.
          exit outer1
        endif
      enddo
    enddo outer1

    return
  end Function ContainsTargetAtom

  !************************************************************************
  !> EvaluateGraphError
  !!
  !! Calculates the error measure between the final graph cx and the
  !! target graph cx_end. The error is defines as the sum of squared
  !! differences for each graph element. If errflag = .TRUE., it means
  !! that the reaction path itself is not valid and we assign it
  !! and arbitrarily large error.
  !!
  !! cx - Final graph in reaction string.
  !! cx_end - Target graph for reaction string.
  !! error - Calculated graph error.
  !! errflag - Logical flag indicating whether the reaction string is valid
  !!           or not.
  !************************************************************************
  !
  Subroutine EvaluateGraphError( cx, cx_end, error, errflag )
    implicit none
    integer :: na, i, j, lwork, info,ic,nc,id1,id2
    type(cxs) :: cx, cx_end
    logical :: errflag
    real(8) :: error, dx, ltemp
    real(8) :: Aend(NAMAX,NAMAX), A(NAMAX,NAMAX)
    real(8), allocatable :: work(:), eigval(:), eigval_end(:)
    real(8) :: val(NAMAX), val_end(NAMAX)
    real(8), allocatable:: dist(:,:),dsp(:,:),dist_end(:,:),dsp_end(:,:)
    character*2 :: label


    if (errflag) then

      error = BIG

    else

      ! Calulate graph error, accounting for permutational invariance.....
      !
      Call CompareGraphs(cx, cx_end, error)

    endif

    return
  end Subroutine EvaluateGraphError


  !************************************************************************
  !> CompareGraphs
  !!
  !! Compare two separate structures via their graphs. In particular,
  !! we compare the eigenvalues of the atomic-weighted graph Coulomb
  !! matrix, with the distances between atoms given by the shortest-path
  !! distances as calculated from the graph.
  !!
  !! cx1, cx2 - Two structures to compare. Note that the graphs MUST have
  !!            been already calculated.
  !! error - Calculate mean-square error between eigenvalues.
  !!
  !************************************************************************
  !
  Subroutine CompareGraphs(cx1, cx2, error)
    implicit none
    integer :: na, i, j, lwork, info,ic,nc,id1,id2, nel,vsum,ifound,k,iel1
    integer :: i1,j1,ina
    type(cxs) :: cx1, cx2
    logical :: errflag
    real(8) :: error, dx, ltemp
    real(8) :: A1(NAMAX,NAMAX), A2(NAMAX,NAMAX)
    real(8), allocatable :: work(:), eigval1(:), eigval2(:),ercomp(:)
    real(8) :: val1(NAMAX), val2(NAMAX)
    real(8), allocatable:: dist1(:,:),dsp1(:,:),dist2(:,:),dsp2(:,:)
    integer, allocatable :: vhist(:,:,:),vhist_targ(:,:,:)
    integer, allocatable :: nelid(:)
    character*2 :: label

    na = cx1%na
    error = 0.d0


    ! IGFUNC = 0: Use element-wise comparison of graphs (NOT PERMUTATIONALLY INVARIANT, AFFECTED
    ! BY ORDER OF ATOMS IN INPUT FILES):
    !
    if (igfunc == 0) then

      do i = 1, na-1
        do j = i+1,na
          error = error + (cx1%graph(i,j)-cx2%graph(i,j))**2
        enddo
      enddo

    else if (igfunc == 1 .or. igfunc == 3) then

      ! allocate( dist1(na,na), dsp1(na,na) )
      !
      ! do i = 1, na
      !    do j = i, na
      !      if ( cx1%graph(i,j) == 1) then
      !        dist1(i,j) = 1.d0
      !        dist1(j,i) = 1.d0
      !      else
      !        dist1(i,j) = BIG
      !        dist1(j,i) = BIG
      !      endif
      !    enddo
      !  enddo
      !  Call GetShortestPaths(na,dist1,dsp1)
      !
      ! allocate( dist2(na,na), dsp2(na,na) )
      !  do i = 1, na
      !    do j = i, na
      !      if ( cx2%graph(i,j) == 1) then
      !        dist2(i,j) = 1.d0
      !        dist2(j,i) = 1.d0
      !      else
      !        dist2(i,j) = BIG
      !        dist2(j,i) = BIG
      !      endif
      !    enddo
      !  enddo
      !  Call GetShortestPaths(na,dist2,dsp2)

      ! Calculate weighted path-distance matrices.
      !
      do i = 1, na
        label = cx1%atomlabel(i)
        id1 = LabelToNumber(label)
        do j = i, na

          label = cx1%atomlabel(j)
          id2 = LabelToNumber(label)

          if (igfunc == 1) then
            A1(i,j) = dble(id1 * id2) * cx1%graph(i,j)
            A2(i,j) = dble(id1 * id2) * cx2%graph(i,j)
          else if (igfunc == 3) then
            A1(i,j) = cx1%graph(i,j)
            A2(i,j) = cx2%graph(i,j)
          endif

          if (i /= j) then
            A1(j,i) = A1(i,j)
            A2(j,i) = A2(i,j)
          endif

        enddo

        if (igfunc == 3) then
          A1(i,i) = id1
          A2(i,i) = id1
        endif

      enddo

      ! Diagonalize.
      !
      lwork = 3*na
      allocate(work(lwork), eigval1(na), eigval2(na) )
      Call DSYEV( 'N', 'U', na, A1, NAMAX, eigval1, work,lwork,info)
      Call DSYEV( 'N', 'U', na, A2, NAMAX, eigval2, work,lwork,info)

      ! Calculate error in eigenvalues.
      !
      do i = 1, na
        error = error + (eigval1(i)-eigval2(i))**2
      enddo
      error = error / dble(na)
      !  error = error * 1d-3

      deallocate(work)
      ! deallocate(dsp1)
      ! deallocate(dist1)
      ! deallocate(dsp2)
      ! deallocate(dist2)
      deallocate(eigval1)
      deallocate(eigval2)


      ! Histogram-based method
    else if (igfunc == 2) then

      ! Count number of elements.
      !
      allocate(nelid(NELMAX))
      nel = 1
      label = cx2%atomlabel(1)
      id1 = LabelToNumber(label)
      nelid(1) = id1
      do i = 2, na
        label = cx2%atomlabel(i)
        id1 = LabelToNumber(label)
        ifound = 0
        do j = 1, nel
          if (id1 == nelid(j))ifound = 1
        enddo
        if (ifound == 0) then
          nel = nel + 1
          nelid(nel) = id1
        endif
      enddo

      ! Allocate space for valence histograms for each element.
      !
      allocate( vhist(nel,nel,0:nvalmax), vhist_targ(nel,nel,0:nvalmax) )
      vhist = 0
      vhist_targ = 0

      ! Calculate valence histograms.
      !
      do i = 1, na  ! For each atom...

        label = cx1%atomlabel(i)
        id1 = LabelToNumber(label)
        do j = 1, nel
          if (nelid(j) == id1) iel1 = j
        enddo
        do j = 1, nel      ! Loop over elements....
          vsum = 0         ! How many element j is i bonded to?
          do k = 1, na
            label = cx1%atomlabel(k)
            id1 = LabelToNumber(label)
            if (nelid(j) == id1) then
              vsum = vsum + cx1%graph(i,k)
            endif
          enddo
          vhist(iel1,j,vsum) = vhist(iel1,j,vsum) + 1
        enddo

        label = cx2%atomlabel(i)
        id1 = LabelToNumber(label)
        do j = 1, nel
          if (nelid(j) == id1) iel1 = j
        enddo
        do j = 1, nel      ! Loop over elements....
          vsum = 0         ! How many element j is i bonded to?
          do k = 1, na
            label = cx2%atomlabel(k)
            id1 = LabelToNumber(label)
            if (nelid(j) == id1) then
              vsum = vsum + cx2%graph(i,k)
            endif
          enddo
          vhist_targ(iel1,j,vsum) = vhist_targ(iel1,j,vsum) + 1
        enddo
      enddo


      error=0.d0
      do i = 1, nel
        do k = 1, nel
          do j = 0, NVALMAX
            error = error + dble(vhist(i,k,j) - vhist_targ(i,k,j))**2
          enddo
        enddo
      enddo

      deallocate(nelid)
      deallocate(vhist)
      deallocate(vhist_targ)

      ! This version is for when the final structure of a target molecule
      ! is given, rather than the full n-atom structure. This will enable
      ! search for paths leading to a defined molecule from a collection
      ! of molecules.
      !
    else if (igfunc == 4) then

      ! Calculate molecules for current graph.
      !
      Call Getmols(cx1)

      ! Loop over all molecules in cx1 and compare each to cx2
      ! (which is just a single molecule).
      !
      allocate( ercomp(cx1%nmol) )

      if (cx2%nmol /= 1)stop 'Target structure must have 1 molecule only for igfunc=4'


      do k = 1, cx1%nmol

        ! First add the contribution to the error from differences in the
        ! number of atoms.....
        !
        ercomp(k) = ( cx1%namol(k) - cx2%namol(1) )**2

        ! Now add permutationally invariant mass-weighted graph difference..
        !
        ! do i = 1, cx1%namol(k)
        !   i1 = cx1%molid(k,i)
        !   label = cx1%atomlabel(i1)
        !   id1 = LabelToNumber(label)
        !   do j = 1, cx1%namol(k)
        !     j1 = cx1%molid(k,j)
        !     label = cx1%atomlabel(j1)
        !     id2 = LabelToNumber(label)
        !     A1(i,j) = cx1%graph(i1,j1)
        !   enddo
        !   A1(i,i) = id1
        ! enddo
        !
        ! do i1 = 1, cx2%namol(1)
        !   label = cx2%atomlabel(i1)
        !   id1 = LabelToNumber(label)
        !   do j1 = 1, cx2%namol(1)
        !     label = cx2%atomlabel(j1)
        !     id2 = LabelToNumber(label)
        !     A2(i1,j1) = cx2%graph(i1,j1)
        !   enddo
        !   A2(i1,i1) = id1
        ! enddo


        allocate( dist1(cx1%namol(k),cx1%namol(k)), dsp1(cx1%namol(k),cx1%namol(k)) )

        do i = 1, cx1%namol(k)
          i1 = cx1%molid(k,i)
          do j = i, cx1%namol(k)
            j1 = cx1%molid(k,j)
            if ( cx1%graph(i1,j1) == 1) then
              dist1(i,j) = 1.d0
              dist1(j,i) = 1.d0
            else
              dist1(i,j) = BIG
              dist1(j,i) = BIG
            endif
          enddo
        enddo
        Call GetShortestPaths(cx1%namol(k),dist1,dsp1)
        do i = 1, cx1%namol(k)
          i1 = cx1%molid(k,i)
          label = cx1%atomlabel(i1)
          id1 = LabelToNumber(label)
          do j = 1, cx1%namol(k)
            A1(i,j) = dsp1(i,j)
          enddo
          A1(i,i) = id1
        enddo


        allocate( dist2(cx2%namol(1),cx2%namol(1)), dsp2(cx2%namol(1),cx2%namol(1)) )
        do i = 1, cx2%namol(1)
          i1 = cx2%molid(1,i)
          do j = i, cx2%namol(1)
            j1 = cx2%molid(1,j)
            if ( cx2%graph(i1,j1) == 1) then
              dist2(i,j) = 1.d0
              dist2(j,i) = 1.d0
            else
              dist2(i,j) = BIG
              dist2(j,i) = BIG
            endif
          enddo
        enddo
        Call GetShortestPaths(cx2%namol(1),dist2,dsp2)
        do i = 1, cx2%namol(1)
          i1 = cx2%molid(1,i)
          label = cx2%atomlabel(i1)
          id1 = LabelToNumber(label)
          do j = 1, cx2%namol(1)
            A2(i,j) = dsp2(i,j)
          enddo
          A2(i,i) = id1
        enddo
        deallocate(dist1)
        deallocate(dist2)
        deallocate(dsp1)
        deallocate(dsp2)

        ! Diagonalize each matrix...
        !
        lwork = 3*max(cx1%namol(k),cx2%namol(1))
        allocate(work(lwork), eigval1(cx1%namol(k)), eigval2(cx2%namol(1)) )
        Call DSYEV( 'N', 'U', cx1%namol(k), A1, NAMAX, eigval1, work,lwork,info)
        Call DSYEV( 'N', 'U', cx2%namol(1), A2, NAMAX, eigval2, work,lwork,info)

        ! Compare eigenvalues....
        ina = min( cx1%namol(k), cx2%namol(1) )

        do i = 1, ina
          ercomp(k) = ercomp(k) + (eigval1(i)-eigval2(i))**2
        enddo
        deallocate(work)
        deallocate(eigval1)
        deallocate(eigval2)
      enddo

      ! The final error is the minimum value of the error across all
      ! molecules.
      !
      error = 1d10
      do i = 1, cx1%nmol
        if (ercomp(i) < error) then
          error = ercomp(i)
          cx1%itargetmol = i
        endif
      enddo
      if (cx1%itargetmol == 0)stop '* Error: no closest match found in CompareGraphs'

      deallocate(ercomp)

    else
      stop '* ERROR: Unknown igfunc in CompareGraphs in pathfinder.f90'
    endif

    return
  end Subroutine CompareGraphs




  !************************************************************************
  !> PropagateGraphs
  !!
  !! Using the string of nrxn graph moves and their atom labels,
  !! this routine starts from the connectivity graph of cx_start
  !! and generates the sequence of nrxn graphs.
  !!
  !! cx_start - Chemical structure object, contains the starting graph.
  !! cx(1:nrxn) - Sequence of chemical structure objects generated
  !!              by graph moves.
  !! nrxn - Number of moves in the reaction string.
  !! movenum(:) - Integer move number for reaction i.
  !! moveatoms(:,:) - Indices of atoms involved in move i.
  !!
  !************************************************************************
  !
  Subroutine PropagateGraphs( cx_start, cx, nrxn, movenum, moveatoms, errflag, error )
    implicit none
    integer :: nrxn, i, j, k, irxn
    integer :: imove, sum,ii,jj,ifound
    type(cxs) :: cx_start, cx(nrxn)
    integer :: movenum(nrxn), moveatoms(nrxn,NAMOVEMAX)
    logical, allocatable :: atomchange(:), bondchange(:,:)
    integer :: nrx, na,rxindex(NAMAX), ix(NAMOVEMAX)
    logical :: errflag
    real(8) :: error


    ! Set up allowed change matrices.
    !
    na = cx_start%na
    allocate( atomchange(na) )
    allocate( bondchange(na,na) )


    ! Set up error flag indicating forbidden move.
    !
    errflag = .FALSE.


    ! Loop through the reactions, moving from one graph to the next...
    !
    outer: do irxn = 1, nrxn

      ! Get the allowed atoms for the current graph, and copy the graph from the
      ! previous step.
      !
      if (irxn == 1) then
        Call SetReactiveIndices( bondchange, atomchange, na, cx_start, rxindex, nrx )
        cx(irxn)%graph(:,:) = cx_start%graph(:,:)
      else
        Call SetReactiveIndices( bondchange, atomchange, na, cx(irxn-1), rxindex, nrx )
        cx(irxn)%graph(:,:) = cx(irxn-1)%graph(:,:)
      endif


      ! Identify the graph move and put the reactive atoms into the ix(:) array.
      !
      imove = movenum(irxn)

      ! If imove = 0 (it's a null move), then we've already copied the previous step
      ! graph so we can just move on to the next reaction.
      !
      if (imove == 0) then
        cycle outer
      endif


   !   write(6,*)'HERE 1',imove

      ! If we get here, we're actually doing a reaction...
      !
      do i = 1, namove(imove)
        ix(i) = moveatoms(irxn,i)
      enddo

      ! Check reactivity indices.
      !
      do i = 1, namove(imove)
  !      print*,'atomcheck:',i,ix(i),atomchange(ix(i))
        if (.not.atomchange(ix(i))) then
          errflag = .TRUE.
          exit outer
        endif
      enddo

!      print*,'HERE2'

!      do i = 1, namove(imove)-1
!        do j = i+1, namove(imove)
! !         print*,'bondcheck:',ix(i),ix(j),bondchange(ix(i),ix(j))
!          if (.not.bondchange(ix(i),ix(j))) then
!            errflag = .TRUE.
!            exit outer
!          endif
!        enddo
!      enddo


      ! Improved version
      do i = 1, namove(imove)-1
        do j = i+1, namove(imove)
          if (gmstart(irxn,i,j).eq.0.and.gmend(irxn,i,j).eq.0) then
            if (.not.bondchange(ix(i),ix(j))) then
              errflag = .TRUE.
              exit outer
            endif
          endif
        enddo
      enddo


   !   print*,'HERE3'


      ! Check that the selected atoms match the graph.
      !
      do i = 1, namove(imove)
        do j = 1, namove(imove)
          if (cx(irxn)%graph(ix(i),ix(j)) /= gmstart(imove,i,j))then
            errflag = .TRUE.
            exit outer
          endif
        enddo
      enddo


      ! Check that the labels match the target labels....
      !
      !movelabel(ngmove,i),i=1,namove(ngmove)
      !
      do i = 1, namove(imove)
        if (trim( movelabel(imove,i) ) /= '*') then
          if (trim(movelabel(imove,i)(1:2)) == '*!') then
            if (trim(movelabel(imove,i)(3:)) == trim(cx(irxn)%atomlabel(ix(i)))) then
              errflag = .TRUE.
              exit outer
            endif
          else
            if (trim(movelabel(imove,i)) /= trim(cx(irxn)%atomlabel(ix(i))) ) then
              errflag = .TRUE.
              exit outer
            endif
          endif
        endif
      enddo

      ! Check that the move involves atoms in essentialatomsinmols.....
      !
      if (nessentialatomsinmols > 0) then

        ! Check that any of the essentialatomsinmols are bonded (or equal to) one of the reactive atoms.
        !
        ifound = 0
        outer3: do i = 1, namove(imove)
          ii = ix(i)
          do j = 1, nessentialatomsinmols
            jj = essentialatomsinmols(j)
            if (irxn == 1) then
              if (cx_start%graph(ii,jj) == 1 .or. ii == jj) then
                ifound = 1
                exit outer3
              endif
            else
              if (cx(irxn-1)%graph(ii,jj) == 1 .or. ii == jj) then
                ifound = 1
                exit outer3
              endif
            endif
          enddo
        enddo outer3
        if (ifound == 0) then
          errflag = .TRUE.
          exit outer
        endif
      endif



      ! Perform the graph move by updating the graph.
      !
      do i = 1, namove(imove)
        do j = 1, namove(imove)
          cx(irxn)%graph(ix(i), ix(j)) = gmend(imove,i,j)
        enddo
      enddo


      ! Check that the final valences are sensible.
      !
      do i = 1, na
        sum = 0
        do j = 1, na
          if (i/=j) then
            sum = sum + cx(irxn)%graph(i,j)
          endif
        enddo

        do k = 1, nvalcon
          if (trim(cx(irxn)%atomlabel(i)) == trim(valatom(k))) then
            if (sum < valrange(k,1)) then
              !    error = error + atemp*(sum-valrange(k,1))**2
              errflag = .TRUE.
              exit outer
            else if ( sum > valrange(k,2))then

              errflag = .TRUE.
              exit outer
              !error = error + atemp*(sum-valrange(k,2))**2

            endif
          endif
        enddo
      enddo


      ! Check that any allowed bonding constraints are not violated.
      !
      do i = 1, na
        do k = 1, nallowbonds
          if (trim(cx(irxn)%atomlabel(i)) == trim( allowbondsatom(k,1) )) then
            sum = 0
            do j = 1, na
              if (i/=j .and. (trim(cx(irxn)%atomlabel(j)) == trim(allowbondsatom(k,2))) ) then
                sum = sum + cx(irxn)%graph(i,j)
              endif
            enddo
            if (sum > allowbondsmax(k)) then
              !error = error + atemp*(sum - allowbondsmax(k))**2
              errflag = .TRUE.
              exit outer
            endif
          endif
        enddo
      enddo

    enddo outer

    ! Clean up memory.
    !
    deallocate( atomchange )
    deallocate( bondchange )


    return
  end Subroutine PropagateGraphs



  !************************************************************************
  !> SelectMoveAtoms
  !!
  !! For move number imove, this subroutine selects new atoms for the
  !! move using the pre-defined list of reactive atom indices in
  !! rxindex.
  !!
  !! imove - A given movenumber. This corresponds to the index of
  !!         a move in the movefile.
  !! moveatoms(nrx,NAMOVEMAX) - Indices of atoms which are moved during
  !!                            a given graph-move.
  !! nrxn - Number of reactions in string.
  !! irxn - The reaction number which we're trying to change.
  !! rxindex(MAMAX) - Indices of reactive atoms in current graph. This
  !!                  array must be previously determined using
  !!                  SetReacctiveIndices() from module gds.
  !! nrx - Number of reactive atoms stored in rxindex(:).
  !! fail - logical flag indicating whether or not this subroutine has
  !!        successfully managed to find new atoms for the move imove.
  !! atomchange(na) - Logical array indicating which atoms can react.
  !! bondchange(na,na) - Logical array indicating which bonds can change.
  !! na - Number of atoms.
  !! cx(nrxn) - Chemical structure objects for each node in reaction
  !!            string.
  !! cx_start - Chemical structure object for start configuration.
  !!
  !************************************************************************
  !
  Subroutine SelectMoveAtoms(imove,moveatoms,nrxn,irx,rxindex,nrx,fail,atomchange, &
    bondchange,na,cx,cx_start)
    implicit none
    integer :: i, j, imove,na,ii,ij,k,imol
    integer :: itry,nrxn, ifound, jj
    integer :: rxindex(NAMAX), nrx, movenum_store, moveatoms_store(NAMOVEMAX)
    logical :: flag,atomchange(na), bondchange(na,na)
    integer :: moveatoms(nrxn,NAMOVEMAX), ix(NAMAX),irx
    type(cxs) :: cx(nrxn), cx_start
    logical :: fail
    real(8) :: ir

    ! Calculate molecules.
    !
    if (irx == 1) then
      Call getmols(cx_start)
    else
      Call getmols(cx(irx-1))
    endif


    ! Set fail flag.
    !
    fail = .FALSE.

    ! Loop over attempts to find atoms - cap this at 5000...
    !
    itry = 0
    inner2: do while (itry < 5000)  !!!????

      itry = itry + 1

      ! Select atoms for this move from reactive indices:
      !
      do i = 1, namove(imove)
      99  call random_number(ir)

        ii = int( dble(nrx) * ir) + 1
        ix(i) = rxindex(ii)
        do j = 1, i-1
          if (ix(i) == ix(j))goto 99
        enddo
      enddo


      ! Check reactivity flags:
      !
      do i = 1, namove(imove)
        if (.not.atomchange(ix(i)))cycle inner2
      enddo
      do i = 1, namove(imove)-1
        do j = i+1, namove(imove)
          if (.not.bondchange(ix(i),ix(j))) cycle inner2
        enddo
      enddo

      ! Check that the selected atoms match the graph.
      !
      do i = 1, namove(imove)
        do j = 1, namove(imove)
          if (irx == 1) then
            if (cx_start%graph(ix(i),ix(j)) /= gmstart(imove,i,j)) cycle inner2
          else
            if (cx(irx-1)%graph(ix(i),ix(j)) /= gmstart(imove,i,j)) cycle inner2
          endif
        enddo
      enddo


      ! Check that these moveatoms obey the move graph...
      !
      do i = 1, namove(imove)
        if (trim( movelabel(imove,i) ) /= '*') then
          if (irx == 1) then
            if (trim(movelabel(imove,i)(1:2)) == '*!') then
              if (trim(movelabel(imove,i)(3:)) == trim(cx_start%atomlabel(ix(i)))) cycle inner2
            else
              if (trim(movelabel(imove,i)) /= trim(cx_start%atomlabel(ix(i))) ) cycle inner2
            endif
          else
            if (trim(movelabel(imove,i)(1:2)) == '*!') then
              if (trim(movelabel(imove,i)(3:)) == trim(cx(irx-1)%atomlabel(ix(i)))) cycle inner2
            else
              if (trim(movelabel(imove,i)) /= trim(cx(irx-1)%atomlabel(ix(i))) ) cycle inner2
            endif
          endif
        endif
      enddo

      ! Check essentialmoveatoms.....
      !
      if (nessentialatoms > 0) then
        ifound = 0
        do i = 1, namove(imove)
          if (essentialmoves(ix(i)) .eqv. .TRUE.) then
            ifound = 1
          endif
        enddo
        if (ifound == 0) cycle inner2
      endif


      ! Are any of the atoms in a molecule involving an essential atom?
      !
      if (nessentialatomsinmols > 0) then

        ! Check that any of the essentialatomsinmols are bonded (or equal to) one of the reactive atoms.
        !
        ifound = 0
        outer3: do i = 1, namove(imove)
          ii = ix(i)
          do j = 1, nessentialatomsinmols
            jj = essentialatomsinmols(j)
            if (irx == 1) then
              if (cx_start%graph(ii,jj) == 1 .or. ii == jj) then
                ifound = 1
                exit outer3
              endif
            else
              if (cx(irx-1)%graph(ii,jj) == 1 .or. ii == jj) then
                ifound = 1
                exit outer3
              endif
            endif
          enddo
        enddo outer3
        if (ifound == 0) cycle inner2

        !
        ! ifound = 0
        ! outer3: do i = 1, namove(imove)
        !   ii = ix(i)
        !   imol = 0
        !   ! Identify molecule for move atom ii.
        !   !
        !   if (irx == 1) then
        !     molsearch: do j = 1, cx_start%nmol
        !       do k = 1, cx_start%namol(j)
        !         if (ii == cx_start%molid(j,k)) then
        !           imol = j
        !           exit molsearch
        !         endif
        !       enddo
        !     enddo molsearch
        !   else
        !     molsearch2: do j = 1, cx(irx-1)%nmol
        !       do k = 1, cx(irx-1)%namol(j)
        !         if (ii == cx(irx-1)%molid(j,k)) then
        !           imol = j
        !           exit molsearch2
        !         endif
        !       enddo
        !     enddo molsearch2
        !   endif
        !
        !   if (imol ==0)stop 'ERROR'
        !
        !   ! Now check that at least one of the other atoms in imol is one of the listed essentialatomsinmols.
        !   !
        !   if (irx == 1) then
        !     do j = 1, cx_start%namol(imol)
        !       ij = cx_start%molid(imol,j)
        !       do k = 1, nessentialatomsinmols
        !         if (ij == essentialatomsinmols(k)) then
        !           ifound = 1
        !           exit outer3
        !         endif
        !       enddo
        !     enddo
        !   else
        !     do j = 1, cx(irx-1)%namol(imol)
        !       ij = cx(irx-1)%molid(imol,j)
        !       do k = 1, nessentialatomsinmols
        !         if (ij == essentialatomsinmols(k)) then
        !           ifound = 1
        !           if (irx<=3)print*,'HERE222:',irx,i,ix(i),cx(irx-1)%molid(imol,1:cx(irx-1)%namol(imol))
        !           exit outer3
        !         endif
        !       enddo
        !     enddo
        !   endif
        ! enddo outer3
        ! if (ifound == 0) cycle inner2
      endif


      ! Once all checks are done, assign the new atoms to the move.
      !
      do i = 1, namove(imove)
        moveatoms(irx,i) = ix(i)
        !if (irx <= 3)print*,'MOVE: ',irx,i,moveatoms(irx,i)
      enddo

      exit inner2

    enddo inner2


    ! If we get here and itry >= 5000, then we've exceeded the number of allowed
    ! attempts to locate atoms....need to return a fail flag and restore
    !
    if (itry >= 5000) then
      fail = .TRUE.
    endif

    return
  end Subroutine SelectMoveAtoms



  !************************************************************************
  !> PathBondEnergy
  !!
  !! Using a simple additive model, calculate the bonding energy along the
  !! entire path, using typical bond energies....
  !!
  !! cx_start - Starting chemical structure.
  !! cx(:) - The set of structures along the path.
  !! nrxn - Number of reactions in string.
  !! vbe - Returned additive energy.
  !!
  !
  !************************************************************************
  !
  Subroutine PathBondEnergy( cx_start, cx, nrxn, vbe, movenum,errflag )
    implicit none
    type(cxs) :: cx_start, cx(nrxn)
    integer :: i,j,k,nrxn, n, movenum(nrxn)
    integer :: irxn,i1,j1,nbreak,nform
    integer, allocatable :: dgraph(:,:)
    real(8) :: vbe, bondstr(NTYPEMAX,NTYPEMAX)
    logical :: errflag

    ! Temporary - bond-strengths.
    !
    bondstr(:,:) = -1.d0
    bondstr(1,1) = 435.d0 / au_to_kjmol
    bondstr(6,6) = 400.d0 / au_to_kjmol
    bondstr(6,1) = 410.d0  / au_to_kjmol
    bondstr(1,6) = bondstr(6,1)
    bondstr(6,8) = 450.d0 / au_to_kjmol
    bondstr(8,6) = bondstr(6,8)
    bondstr(1,8) = 460.d0 / au_to_kjmol
    bondstr(8,1) = bondstr(1,8)
    bondstr(8,8) = 500.d0 / au_to_kjmol
    bondstr(78,6) = 580.d0 / au_to_kjmol
    bondstr(6,78) = bondstr(78,6)
    bondstr(78,8) = 420.d0 / au_to_kjmol
    bondstr(8,78) = bondstr(78,8)
    bondstr(78,1) = 330.d0 / au_to_kjmol
    bondstr(1,78) = bondstr(78,1)

    ! Cobalt = 27
    bondstr(27,6) = 180.d0 / au_to_kjmol
    bondstr(6,27) = bondstr(27,6)

    bondstr(27,1) = 220.d0 / au_to_kjmol
    bondstr(1,27) = bondstr(27,1)

    bondstr(27,8) = 180.d0 / au_to_kjmol
    bondstr(8,27) = bondstr(27,8)

    !
    !********


    ! Get number of atoms, and allocated workspace for difference graph.
    !
    n = cx_start%na
    allocate (dgraph(n,n))


    ! Loop over reaction string, accumulate error
    !
    vbe = 0.d0
    do irxn = 1, nrxn

      ! vbe = vbe + alphavbe * namove(movenum(irxn))**2

      ! Calculate difference graph.
      !
      if ( irxn == 1) then
        do i = 1, n
          do j = 1, n
            dgraph(i,j) = cx(irxn)%graph(i,j) - cx_start%graph(i,j)
          enddo
        enddo
      else
        do i = 1, n
          do j = 1, n
            dgraph(i,j) = cx(irxn)%graph(i,j) - cx(irxn-1)%graph(i,j)
          enddo
        enddo
      endif


      nbreak = 0
      outer2: do i = 1, n
        do j = 1, n
          if (abs(dgraph(i,j))>0)then
            nbreak = nbreak + 1
            !     print*,'REACTIVE: ',irxn,i
            cycle outer2
          endif
        enddo
      enddo outer2
      ! print*,'ACTIVE ATOMS = ',irxn,nbreak
      vbe = vbe + alphavbe * nbreak**2


      ! Calculate error contribution.
      !
      ! do i = 1, n-1
      !    do j = i+1, n
      !
      !       if (dgraph(i,j) /= 0) then
      !
      !          i1 = LabelToNumber(cx_start%atomlabel(i))
      !          j1 = LabelToNumber(cx_start%atomlabel(j))
      !          if (dgraph(i,j)<0.d0) then
      !            vbe = vbe - dble(dgraph(i,j)) * bondstr(i1,j1)
      !          endif
      !          if (bondstr(i1,j1) < 0.d0)then
      !             stop '*** ERROR in PathBondEnergy in pathfinder.f90'
      !          endif
      !
      !       endif
      !    enddo
      ! enddo

      ! New bit...penalize null reactions.
      !
      ! if (movenum(irxn) ==0)vbe = vbe + 0.5d0

    enddo
    return
  end Subroutine PathBondEnergy


  !************************************************************************
  !> GetPathFitness
  !!
  !! Calculate the fitness value associated with a given reaction mechanism.
  !!
  !! This TotalError is defined as the sum of a GraphError, describing how far
  !! away a propagated graph is from a target graph, as well as a vbe contribution
  !! which includes energy contributions to the path fitness.
  !!
  !! Note that the energy contribution is only calculated when alphavbe > 1d-3.
  !!
  !! cx_start - The chemical structure object for the reactants.
  !! cx_end - The chemical structure object for the products.
  !! cx(:) - Chemical structure objects for intermediates along path.
  !! nrxn - Number of reactions in path.
  !! novenum(nrxn) - Move numbers for each reaction in the path.
  !! moveatomos(nrxn,NAMOVEMAX) - Atoms which are acted upon in each reaction.
  !! errflag - Logical flag indicating whether or not a chemically-sensible reaction
  !!           path has been generated.
  !! GraphError - Error term arising due to differences between propagated graph
  !!              and target graph.
  !! TotalError - GraphError + alphavbe * vbe.
  !! vbe - Calculated bonding energy.
  !!
  !************************************************************************
  !
  Subroutine GetPathFitness(cx_start, cx_end, cx, nrxn, movenum, moveatoms, errflag, &
    GraphError, TotalError,vbe )
    implicit none
    integer :: nrxn
    type(cxs) :: cx_start, cx_end, cx(nrxn)
    integer :: moveatoms(nrxn,NAMOVEMAX), movenum(nrxn)
    real(8) :: GraphError, TotalError, vbe, error
    logical :: errflag

    ! Propagate the graph from the start to generate the final graph:
    !
    GraphError = 0.d0
    Call PropagateGraphs( cx_start, cx, nrxn, movenum, moveatoms, errflag, GraphError )

  !  print*,'HERE: ',errflag
  !  stop

    ! Evaluate the distance from the target graph:
    !
    Call EvaluateGraphError( cx(nrxn), cx_end, grapherror, errflag )

    ! If alphavbe is greater than zero, then calculate the effective energy along the reaction path.
    !
    if ( alphavbe > 1d-3) then
      Call PathBondEnergy( cx_start, cx, nrxn, vbe, movenum,errflag )
      TotalError = GraphError + alphavbe * vbe
    else
      TotalError = GraphError
    endif

    return
  end Subroutine GetPathFitness


  !************************************************************************
  !> TrimPath
  !!
  !! When the target molecule is detected in a reaction, the rest of the
  !! movenums are set to zero. In other words, if the target is found at
  !! some reaction step before the end, we then set the rest of the
  !! reactions to zero because we don't care what they are any more..
  !!
  !! cx_start - The chemical structure object for the reactants.
  !! cx_end - The chemical structure object for the products.
  !! cx(:) - Chemical structure objects for intermediates along path.
  !! nrxn - Number of reactions in path.
  !! novenum(nrxn) - Move numbers for each reaction in the path.
  !! moveatoms(nrxn,NAMOVEMAX) - Atoms which are acted upon in each reaction.
  !!
  !************************************************************************
  !
  Subroutine TrimPath(cx_start, cx_end, cx, nrxn, movenum, moveatoms, iend )
    implicit none
    integer :: nrxn, i, j, iend
    type(cxs) :: cx_start, cx_end, cx(nrxn)
    integer :: moveatoms(nrxn,NAMOVEMAX), movenum(nrxn)
    real(8) :: grapherror
    logical :: errflag

    write(logfile,'(/"* Checking when target molecule is created...trimming path..."/)')

    ! Propagate the graph from the start to generate the final graph:
    !
    Call PropagateGraphs( cx_start, cx, nrxn, movenum, moveatoms, errflag, GraphError )

    ! Evaluate the distance from the target graph:
    !
    outer1: do i = 1, nrxn
      Call EvaluateGraphError( cx(i), cx_end, grapherror, errflag )
      write(logfile,'("- Reaction: ",1x,i3,": GraphError = ",1x,f16.8)')i,grapherror

      if (grapherror <= 1d-5) then
        write(logfile,'(/"* TARGET FOUND AFTER REACTION ",i3)')i
        write(logfile,'("* Setting remainder of reactions to null moves...")')
        do j = i+1, nrxn
          movenum(j) = 0
        enddo
        iend = i
        exit outer1
      endif
    enddo outer1

    ! Now we need to set the target graphs of the remaining reactions to be the
    ! same as the graph of the final frame when the products was formed.
    !
    do j = iend+1, nrxn
      Call CopyCXS( cx(iend), cx(i) )
    enddo

    write(logfile,'(/"* Successfully finished path trimming...."/)')

    return
  end Subroutine TrimPath


  !************************************************************************
  !> PrintMolsAlongPath
  !!
  !! Outputs (to logfile) the molecular structure present at each step
  !! along a reaction-mechanism,so that the user can get a feel for what is going on in the path.
  !!
  !! nrxn - Number of reactions in path.
  !! cx_start - The chemical structure object for the reactants.
  !! cx_end - The chemical structure object for the products.
  !! cx(:) - Chemical structure objects for intermediates along path.
  !! movenum(nrxn) - Move numbers for each reaction in the path.
  !!
  !************************************************************************
  !
  Subroutine PrintMolsAlongPath( nrxn, cx_start, cx_end, cx, movenum,chargemove,changecharges,moveatoms )
    implicit none
    integer :: nrxn, i, j, movenum(nrxn),chargemove(nrxn,NMOLMAX),moveatoms(nrxn,NAMOVEMAX)
    type(cxs) :: cx_start, cx_end, cx(nrxn)
    logical :: changecharges

    write(logfile,'("*********  Molecular analysis of final reaction-string  ********"/)')

    ! Print information for start-point:
    !
    write(logfile,'("======== Start-point molecules ========"/)')
    Call GetMols(cx_start)
    do j = 1, cx_start%nmol
      write(logfile,'("- Molecule number:",3x,i4)')j
      write(logfile,'("- Chemical formula:",3x,a/)')MolecularFormula( cx_start, j )
    enddo

    do i = 1, nrxn
      write(logfile,'("=== Reaction number:",3x,i4/)')i
      write(logfile,'("- Selected move number:",3x,i4/)')movenum(i)
      write(logfile,*)"- Atom numbers: ",moveatoms(i,1:namove(movenum(i)))
      Call GetMols(cx(i))
      do j = 1, cx(i)%nmol
        write(logfile,'("- Molecule number:",3x,i4)')j
        write(logfile,'("- Chemical formula:",3x,a/)')MolecularFormula( cx(i), j )
      enddo
      if (changecharges) then
        write(logfile,'("- Charge Moves:")')
        do j = 1, cx(i)%nmol
          write(logfile,*)j,chargemove(i,j)
        enddo
      endif
    enddo

    write(logfile,'("======== End-point molecules ========"/)')
    Call GetMols(cx_end)
    do j = 1, cx_end%nmol
      write(logfile,'("- Molecule number:",3x,i4)')j
      write(logfile,'("- Chemical formula:",3x,a/)')MolecularFormula( cx_end, j )
    enddo
    write(logfile,'("********  Molecular analysis of reaction-string complete *******"/)')

    return
  end Subroutine PrintMolsAlongPath


  !**********************************************************************************************
  !> GraphsToCoords
  !!
  !! Converts a sequence of connectivity graphs to xyz coordinates by optimization
  !! under the action of the graph-restraining potential.
  !!
  !! The coordinates of the reactant structure are contained in cx_start; these are
  !! used as the starting point for a sequence of geometry optimizations. The resulting
  !! structures obey the target bonding graph at each step.
  !!
  !! If optaftermove = .TRUE. in the input file, then geometry optimization of the structures
  !! is also performed AFTER graph-potential optimization.
  !!
  !! cx_start - Chemical structure object for reactants.
  !! cx(nrxn) - Chemical structures along the path.
  !! nrxn - Number of reactions along the path.
  !! printflag - Logical flag indicating whether to print out xyz files or not.
  !! printfile - If printflag = .TRUE., printfile indicates the file to print to.
  !!
  !
  !**********************************************************************************************
  !
  Subroutine GraphsToCoords(cx_start, cx, nrxn, printflag, printfile)
    implicit none
    integer :: nrxn
    type (cxs) :: cx_start, cx(nrxn)
    logical :: printflag, success
    character(len=*) :: printfile
    integer :: j,k,na,irxn,isum,i
    integer, allocatable :: grstore(:,:)
    real(8) :: x,y,z

    ! Set a useful local parameter - number of atoms.
    !
    na = cx_start%na
    allocate( grstore(na,na))


    ! If requested by user, perform geometry optimization for the start point.
    !
    if (optaftermove) then

      ! Store the original graph....
      !
      grstore(:,:) = cx_start%graph(:,:)
      Call AbInitio( cx_start, 'optg', success )
      Call SetCXSconstraints( cx_start ,NDOFconstr, FixedDOF, Natomconstr, FixedAtom)

      ! Check that the new graph matches the original - otherwise atoms have moved during optimization.
      ! If this has happened, we set the energy to some large value.
      !
      Call GetGraph(cx_start)
      isum = 0
      do i = 1, na
        do j = 1, na
          if (cx_start%graph(i,j) /= grstore(i,j))isum=isum+1
        enddo
      enddo
      if (isum /= 0) cx_start%vcalc = 1d4

    endif


    ! First, print out the starting structure.
    !
    if (printflag) then
      open(93,file = trim(printfile), status='unknown')
      write(93,'(i5)')na
      write(93,'("Starting structure.",1x,"Energy in au =",1x,f14.8)') &
      cx_start%vcalc
      do j = 1, na
        x = cx_start%r(1,j) * bohr_to_ang
        y = cx_start%r(2,j) * bohr_to_ang
        z = cx_start%r(3,j) * bohr_to_ang
        write(93,'(a2,2x,3(f14.8,2x))')cx_start%atomlabel(j),x,y,z
      enddo
      call flush(93)
    endif

    ! Loop over reactions, generating coordinates based on changes to previous step.
    !
    do irxn = 1, nrxn

      ! Assign starting coordinates as the coordinates of the prevous step.
      !
      if (irxn == 1) then
        do j = 1, na
          do k = 1, 3
            cx(irxn)%r(k,j) = cx_start%r(k,j)
          enddo
          cx(irxn)%atomlabel(j) = cx_start%atomlabel(j)
        enddo
      else
        do j = 1, na
          do k = 1, 3
            cx(irxn)%r(k,j) = cx(irxn-1)%r(k,j)
          enddo
          cx(irxn)%atomlabel(j) = cx(irxn-1)%atomlabel(j)
        enddo
      endif


      ! Set constraints.
      !
      Call SetCXSconstraints( cx(irxn) ,NDOFconstr, FixedDOF, Natomconstr, FixedAtom)


      ! Optimize coordinates under action of current graph.
      !
      Call GetMols( cx(irxn) )

      !    Call PrintCXSToFile(cx(irxn),'temp1.xyz',1.d0)

      ! Optimize changed graph ic under action of graph restraining potential.
      !
      !  cx(irxn)%dvdr(:,:) = 0.D0
      !  Call OptimizeGRP(cx(irxn),success, gdsrestspring, nbstrength, nbrange, &
      !  kradius, ngdsrelax, gdsdtrelax )

      ! NEW - optimize under double-ended GRP.
      !
      if (irxn == 1) then
        Call OptimizeGRP_DoubleEnded(cx(irxn), cx_start, success, gdsrestspring, nbstrength, nbrange, &
        kradius, ngdsrelax, gdsdtrelax )
      else
        Call OptimizeGRP_DoubleEnded(cx(irxn), cx(irxn-1), success, gdsrestspring, nbstrength, nbrange, &
        kradius, ngdsrelax, gdsdtrelax )
      endif

      !  Call PrintCXSToFile(cx(irxn),'temp2.xyz',1.d0)

      ! Call GetMols( cx(irxn) ) need this?

      ! print*,'NMOL AFTER OptGRP = ',cx(irxn)%nmol
      !  stop
      if (.not.success)then
        !    print*,'WTF'
        !    stop
      endif

      ! If requested by user, perform geometry optimization after move.
      !
      if (optaftermove) then
        grstore(:,:) = cx(irxn)%graph(:,:)

        Call AbInitio( cx(irxn), 'optg', success )
        Call SetCXSconstraints( cx(irxn) ,NDOFconstr, FixedDOF, Natomconstr, FixedAtom)

        Call GetGraph(cx(irxn))
        isum = 0
        do i = 1, na
          do j = 1, na
            if (cx(irxn)%graph(i,j) /= grstore(i,j))then
              !          print*,'PROBLEM:',i,j,cx(irxn)%atomlabel(i),cx(irxn)%atomlabel(j),cx(irxn)%graph(i,j),grstore(i,j)
              !    stop
              isum=isum+1
            endif
          enddo
        enddo
        !  if (isum /= 0) cx(irxn)%vcalc = 1d4

      endif


      ! Output coordinates to final_path.xyz
      !
      if (printflag) then
        write(93,'(i5)')na
        write(93,'("Product of reaction",1x,i3,1x,"Energy in au =",1x,f14.8)') &
        irxn,cx(irxn)%vcalc
        do j = 1, na
          x = cx(irxn)%r(1,j) * bohr_to_ang
          y = cx(irxn)%r(2,j) * bohr_to_ang
          z = cx(irxn)%r(3,j) * bohr_to_ang
          write(93,'(a2,2x,3(f14.8,2x))')cx(irxn)%atomlabel(j),x,y,z
        enddo
        call flush(93)
      endif

    enddo

    if (printflag)close(93)

    deallocate(grstore)

    return
  end Subroutine GraphsToCoords


  !**********************************************************************************************
  !> PrintMechanismPaths
  !!
  !! Prints out nimage intermediate snapshots for each of the generated reactions
  !! in the final reaction-path.
  !!
  !! Note that the end-point coordinates for each reaction are assumed to be stored in
  !! cx_start and cx(nrxn).
  !!
  !! We use linear interpolation between these end-points to generate internal images.
  !!
  !! The results are stored to the files final_path_rx_*.xyz.
  !!
  !! nrxn - Number of reactions in mechanism.
  !! cx_start - Chemical structure object for reactants.
  !! cx(nrxn) - Chemical structure objects for products of each reaction. These come from the
  !!            subroutine GraphsToCoords().
  !! file_root - Root name of printing file.
  !!
  !
  !**********************************************************************************************
  !
  Subroutine PrintMechanismPaths(nrxn, cx_start, cx, file_root,maxbarrier,bsum,movenum)

    implicit none

    integer :: nrxn, na, irxn, j, movenum(NRXNMAX), idum
    type(rxp) :: rp
    type(cxs) :: cx_start, cx(nrxn)
    character (len=7) :: x1
    character (len=8) :: fmt
    character (len=20) :: file_root
    logical :: success,idppguess
    real(8) :: maxbarrier,bsum, brxn

    maxbarrier = -1d6

    ! Set number of atoms.
    !
    na = cx_start%na

    ! Create workspace for reaction-path generation.
    !
    Call NewPath(rp, .FALSE., startfile, endfile, pathfile, nimage, &
    pathinit, .TRUE., na)


    ! Loop over each reaction.
    !
    open(95,file=trim(file_root)//'_energy.dat')
    bsum = 0.d0
    outer3: do irxn = 1, nrxn

      ! Set the CX for the end-point of this reaction.
      !
      if (irxn == 1) then
        rp%cx(1) = cx_start
        do j = 2, nimage-1
          rp%cx(j) = cx_start
        enddo
        rp%cx(nimage) = cx(irxn)
      else
        rp%cx(1) = cx(irxn-1)
        do j = 2, nimage-1
          rp%cx(j) = cx(irxn-1)
        enddo
        rp%cx(nimage) = cx(irxn)
      endif
      Call SetPathConstraints(rp, NDOFconstr, FixedDOF, Natomconstr, Fixedatom)

      ! Set as linear path - this needs to be done to give sensible initial
      ! coordinates to the internal beads before IDPP (if required).
      !
      rp%coeff(:,:,:) = 0.0
      Call FourierToPath( rp )

      ! Use the Fourier coefficients to calculate the initial path.
      !
      if (idpppath) then
        Call FindIDPPPath( rp, NEBIter*250, NEBConv*0.1d0, NEBstep, NEBspring)
      endif


      ! Could put NEB here.....


      ! Print path to file.
      !
      if (irxn <= 9) then
        fmt = '(I1.1)'
        write (x1,fmt)irxn
      else
        fmt = '(I2.2)'
        write (x1,fmt)irxn
      endif
      Call PrintPathToFile(rp,trim(file_root)//'_rx_'//trim(x1)//'.xyz')

    enddo outer3
    close(95)

    ! Clean up memory.
    !
    Call DeletePath(rp)

    return
  end Subroutine PrintMechanismPaths



  !**********************************************************************************************
  !> UpdateMechanism
  !!
  !! Performs a trial move by modifying a current mechanism string by either changing
  !! a move and its atoms, or just changing the atoms of an exisiting move.
  !!
  !! Note that the end-point coordinates for each reaction are assumed to be stored in
  !! cx_start and cx(nrxn).
  !!
  !! nrxn - Number of reactions in mechanism.
  !! cx_start - Chemical structure object for reactants.
  !! cx(nrxn) - Chemical structure objects for products of each reaction. These come from the
  !!            subroutine GraphsToCoords().
  !!
  !
  !**********************************************************************************************
  !
  Subroutine UpdateMechanism(nrxn,movenum,moveatoms,bondchange,atomchange,na,cx_start,cx,rxindex,cyc, &
    movenum_store,moveatoms_store)

    implicit none
    integer :: na, irx, itype, imove, nmove,im,i,j
    real(8) :: rx, ir
    type(cxs) :: cx_start, cx(nrxn)
    integer :: nrxn, movenum(nrxn), moveatoms(nrxn,NAMOVEMAX)
    logical :: bondchange(na,na), atomchange(na), fail, cyc
    integer :: movenum_store(nrxn), moveatoms_store(nrxn,NAMOVEMAX)
    integer :: rxindex(NAMAX), nrx


    ! Set the loop cycle flag for return
    !
    cyc = .FALSE.

    ! Store the previous move.
    !
    movenum_store(:) = movenum(:)
    moveatoms_store(:,:) = moveatoms(:,:)


    ! Decide how many moves to attempt - the maximum is set by nmechmove in the input file.
    !
    if (nmechmove == 1) then
      nmove = 1
    else
      call random_number(ir)
      nmove = int(dble(nmechmove) * ir) + 1
    endif

    do im = 1, nmove

      ! Decide which reaction to change....
      !
      call random_number(ir)

      irx = int( dble(nrxn) * ir ) + 1


      ! Decide which type of move to perform - change atoms or change both move and atoms.
      ! However, if the move-number of reaction irx is zero (i.e. a null reaction) then we must
      ! insist on itype = 2.
      !
      if (movenum(irx) == 0) then
        itype = 2
      else
        call random_number(rx)

        if (rx < 0.5d0) then
          itype = 1
        else
          itype = 2
        endif
      endif



      ! Get the allowed reactivity indices for the graph preceding this reaction...
      !
      if (irx == 1) then
        Call SetReactiveIndices( bondchange, atomchange, na, cx_start, rxindex, nrx )
      else
        Call SetReactiveIndices( bondchange, atomchange, na, cx(irx-1), rxindex, nrx )
      endif


      ! Now change reaction irx to some new reaction.....
      !
      select case(itype)

        !************************************************************
        ! Change atom labels of existing move.
        !************************************************************
        !
      case (1)

        ! Select new reactive atoms from rxindex.
        !
        imove = movenum(irx)
        Call SelectMoveAtoms(imove,moveatoms,nrxn,irx,rxindex,nrx,fail,atomchange, &
        bondchange,na,cx,cx_start)


        ! If we've failed to find atoms, restore the previous and cycle...
        !
        if (fail) then
          movenum(:) = movenum_store(:)
          moveatoms(:,:) = moveatoms_store(:,:)
          cyc = .TRUE.
          return
        endif


        !************************************************************
        ! Change both atom labels and move type.
        !************************************************************
        !
      case (2)

        ! Choose a new move-type - note that this CAN include a null move, so we
        ! don't add "1" at the end....
        !
      51  call random_number(ir)

        movenum(irx) = int( dble(ngmove+1) * ir )         !!!!!!!!!!
        !  51    movenum(irx) = 1 + int( dble(ngmove) * ran2(irun,0.d0,1.d0))

        if (movenum(irx) == movenum_store(irx))goto 51

        ! Select new reactive atoms from rxindex.
        !
        imove = movenum(irx)

        if (imove == 0) then
          moveatoms(irx,1:NAMOVEMAX) = 0
        else

          ! Select new reactive atoms from rxindex.
          !
          Call SelectMoveAtoms(imove,moveatoms,nrxn,irx,rxindex,nrx,fail,atomchange, &
          bondchange,na,cx,cx_start)

          !
          ! if (irx==1)then
          !   if (movenum(irx)==1)then
          !     moveatoms(irx,1) = 1
          !     moveatoms(irx,2) = 5
          !     print*,'HERE:',irx,movenum(irx),movenum_store(irx),moveatoms(irx,1:namove(movenum(irx)))
          !   endif
          ! endif
          !

          ! If we've failed to find atoms, restore the previous and cycle...
          !
          if (fail) then
            movenum(:) = movenum_store(:)
            moveatoms(:,:) = moveatoms_store(:,:)
            cyc = .TRUE.
            return
          endif

        endif

      end select

    enddo


    return
  end Subroutine UpdateMechanism


  !**********************************************************************************************
  !> ReactionErrorCorrection
  !!
  !! Performs error correction using geometry-optimized structures, making sure that the reactions
  !! listed in movenum and moveatoms are actually correct for a given reaction string.
  !!
  !! Note that the end-point coordinates for each reaction are assumed to be stored in
  !! cx_start and cx(nrxn).
  !!
  !!
  !
  !**********************************************************************************************
  !
  Subroutine ReactionErrorCorrection(cx_start, cx, nrxn, movenum, moveatoms, errx)
    implicit none
    integer :: nrxn, irxn, i, j, isum, na, ip,nratom
    integer :: ita(4), it1,it2,it3,i1,j1,k,l,ifound,imove
    integer :: movenum(nrxn), moveatoms(nrxn,NAMOVEMAX)
    integer, allocatable :: movenum_temp(:), moveatoms_temp(:,:), irx(:)
    type(cxs) :: cx_start, cx(nrxn), cxinit
    type(cxs), allocatable :: cxtemp(:)
    logical :: errflag, errx
    real(8) :: error

    ! Print final....
    !
    ! do irxn = 1, nrxn
    !   print*,'**** INITIAL REACTION = ',irxn
    !   print*,'MOVENUM = ',movenum(irxn)
    !   do j = 1, namove(movenum(irxn))
    !     print*,'ATOMS: ',moveatoms(irxn,j)
    !   enddo
    !   print*
    ! enddo

    ! Initialize error flag.
    !
    errx = .FALSE.

    ! Setup up workspace.
    !
    allocate(cxtemp(1))
    Call CopytoNewCXS( cx_start, cxtemp(1) )
    Call SetMass(cxtemp(1))
    Call CopytoNewCXS( cx_start, cxinit )
    Call SetMass(cxinit)

    allocate( movenum_temp(1), moveatoms_temp(1,NAMOVEMAX), irx(NAMOVEMAX) )
    na = cx_start%na

    ! Make sure we have the correct graphs for each structure.
    !
    Call GetGraph( cx_start )
    do irxn = 1, nrxn
      Call GetGraph(cx(irxn))
    enddo


    ! Loop over reactions, checking movenum and moaveatoms data is still correct....
    !
    outer1: do irxn = 1, nrxn        !!!!!!!!!!!!!!!!!

      ! print*,'REACTION: ',irxn
      ! if (movenum(irxn) == 0) then
      !   print*,'ORGINAL = 0 0'
      ! else
      !   print*,'ORGINAL = ',movenum(irxn),moveatoms(irxn,1:namove(movenum(irxn)))
      ! endif

      ! Move the reaction details into the temporary array.
      !
      movenum_temp(1) = movenum(irxn)
      moveatoms_temp(1,:) = moveatoms(irxn,:)

      ! Make a CXS object containing the start-point..
      !
      if (irxn==1)then
        Call CopyCXS(cx_start,cxinit)
      else
        Call CopyCXS(cx(irxn-1),cxinit)
      endif
      Call CopyCXS(cx(irxn),cxtemp(1))       !!!!!! REMOVE
      !
      ! ! Trick the propagategraphs routine to run this reaction..the result will be in cxtemp.
      ! !
      ! Call PropagateGraphs( cxinit, cxtemp, 1, movenum_temp, moveatoms_temp, errflag, error )
      !
      ! print*,'ERRFLAG = ',errflag
      !
      ! ! Now check whether cxtemp is the same graph as the expected cx(irxn)...
      ! !
      ! isum = 0
      ! do i = 1, na-1
      !   do j = i+1, na
      !     isum = isum + abs(cxtemp(1)%graph(i,j) - cx(irxn)%graph(i,j))
      !   enddo
      ! enddo
      !
      !
      ! ! If isum is zero, everything is fine and there is nothing else to do...
      ! !
      ! if (isum ==0) cycle outer1


      ! If isum is not zero, we need to work out which reaction did happen....
      !
      ! Count number of reactive atoms....
      !
      nratom = 0
      irx(:) = 0
      outer2: do i = 1, na
        do j = 1, na
          !      print*,'WTF:',irxn,i,j,cxinit%graph(i,j),cxtemp(1)%graph(i,j)
          if (abs(cxinit%graph(i,j)-cxtemp(1)%graph(i,j)) > 0 ) then
            nratom = nratom + 1
            irx(nratom) = i
            cycle outer2
          endif
        enddo
      enddo outer2
      !  print*,'NRATOM = ',nratom


      ! If nratom is zero, we've found a null move.
      !
      if (nratom == 0) then
        ifound = 1
        movenum(irxn) = 0
        moveatoms(irxn,:) = 0
      else if (nratom==2) then

        ! Check 2-atom moves...
        !
        ifound = 0
        outer4: do imove = 1, ngmove

          if (namove(imove) == 2) then

            if (cxinit%graph(irx(1),irx(2)) == gmstart(imove,1,2)) then
              if (cxtemp(1)%graph(irx(1),irx(2)) == gmend(imove,1,2)) then
                movenum(irxn) = imove
                moveatoms(irxn,1) = irx(1)
                moveatoms(irxn,2) = irx(2)
                ifound = 1
                exit outer4
              endif
            endif
          endif

        enddo outer4

      else if (namove(imove) == 3) then

        ifound = 0
        outer5: do imove = 1, ngmove

          if (namove(imove) == 3) then

            ! Check permutations...
            !
            it1 = irx(1)
            it2 = irx(2)
            it3 = irx(3)
            perm1: do ip = 1, 6
              if (ip == 1) then
                irx(1) = it1
                irx(2) = it2
                irx(3) = it3
              else if (ip == 2) then
                irx(1) = it2
                irx(2) = it1
                irx(3) = it3
              else if (ip == 3) then
                irx(1) = it2
                irx(2) = it3
                irx(3) = it1
              else if (ip == 4) then
                irx(1) = it1
                irx(2) = it3
                irx(3) = it2
              else if (ip == 5) then
                irx(1) = it3
                irx(2) = it2
                irx(3) = it1
              else if (ip == 6) then
                irx(1) = it3
                irx(2) = it1
                irx(3) = it2
              endif

              ! Compare graphs....
              !
              do i = 1, 3
                do j = 1, 3
                  if (cxinit%graph(irx(i),irx(j)) /= gmstart(imove,i,j)) then
                    cycle perm1
                  endif
                enddo
              enddo

              do i = 1, 3
                do j = 1, 3
                  if (cxtemp(1)%graph(irx(i),irx(j)) /= gmend(imove,i,j)) then
                    cycle perm1
                  endif
                enddo
              enddo

              ! If we get here, we've identified the move...
              !
              ifound = 1
              movenum(irxn) = imove
              moveatoms(irxn,1:3) = irx(1:3)
              exit outer5
            enddo perm1
          endif
        enddo outer5

      else if (namove(imove) == 4) then

        ifound = 0
        outer6: do imove = 1, ngmove

          if (namove(imove) == 4) then

            ! Check permutations...
            !
            ita(1) = irx(1)
            ita(2) = irx(2)
            ita(3) = irx(3)
            ita(4) = irx(4)
            loopi: do i = 1, 4
              loopj: do j = 1, 4
                if (i == j)cycle loopj
                loopk: do k = 1, 4
                  if (k == i .or. k ==j)cycle loopk
                  loopl: do l = 1, 4
                    if (l == k .or. l == i .or. l == j)cycle loopl
                    irx(1) = ita(i)
                    irx(2) = ita(j)
                    irx(3) = ita(k)
                    irx(4) = ita(l)

                    ! Compare graphs....
                    !
                    do i1 = 1, 4
                      do j1 = 1, 4
                        if (cxinit%graph(irx(i1),irx(j1)) /= gmstart(imove,i1,j1)) then
                          cycle loopl
                        endif
                      enddo
                    enddo

                    do i1 = 1, 4
                      do j1 = 1, 4
                        if (cxtemp(1)%graph(irx(i1),irx(j1)) /= gmend(imove,i1,j1)) then
                          cycle loopl
                        endif
                      enddo
                    enddo

                    ! If we get here, we've found the move...
                    !
                    ifound = 1
                    movenum(irxn) = imove
                    moveatoms(irxn,1:4) = irx(1:4)
                    exit outer6
                  enddo loopl
                enddo loopk
              enddo loopj
            enddo loopi
          endif
        enddo outer6

      else
        ifound = 0
      endif

      ! If we get here and ifound is zero, we haven't been able to ID the move - need to return
      ! with an error flag.
      !
      if (ifound == 0) then
        print*,'FAILED TO ID REACTION:',irxn
        stop
        errx = .TRUE.
        return
      endif

    enddo outer1

    ! Deallocate workspace
    !
    Call DeleteCXS (cxtemp(1))
    Call DeleteCXS(cxinit)
    deallocate(cxtemp)
    deallocate( movenum_temp)
    deallocate(moveatoms_temp)
    deallocate( irx )

    ! Print final....
    !
    ! do irxn = 1, nrxn
    !   print*,'**** FINAL REACTION = ',irxn
    !   print*,'MOVENUM = ',movenum(irxn)
    !   do j = 1, namove(movenum(irxn))
    !     print*,'ATOMS: ',moveatoms(irxn,j)
    !   enddo
    !   print*
    ! enddo

    return
  end Subroutine ReactionErrorCorrection


  !**********************************************************************************************
  !> UpdateCharges
  !!
  !! Update the charges on molecules along the reaction path. These moves are performed in
  !! addition to reaction changes.
  !!
  !
  !**********************************************************************************************
  !
  Subroutine UpdateCharges(nrxn,cx,chargemove, chargemove_store, errflag)
    implicit none
    integer :: nrxn, irx, imol,i,j,k,nn,ncharge,TotalCharge,ia
    type(cxs) :: cx(nrxn)
    real(8) :: ir,r1
    integer :: chargemove(nrxn,NMOLMAX), chargemove_store(nrxn,NMOLMAX)
    logical :: errflag

    ! Store the current charge moves.
    !
    chargemove_store = chargemove

    ! Decide which reaction to change....
    !
    call random_number(ir)
    irx = int( dble(nrxn) * ir ) + 1

    ! Get molecules
    !
    Call GetMols(cx(irx))

    ! Choose a random molecule.
    !
    call random_number(ir)
    imol = 1 + int( dble(cx(irx)%nmol) * ir )

    ! Choose to increase or decrease charge on molecule imol with equal probability.
    !
    call random_number(r1)

    if (r1 < 0.5d0) then
      cx(irx)%molcharge(imol) = cx(irx)%molcharge(imol) - 1
      chargemove(irx,imol) = -1
    else
      cx(irx)%molcharge(imol) = cx(irx)%molcharge(imol) + 1
      chargemove(irx,imol) = +1
    endif

    ! Check whether the charge on any molecule exceeds min/max allowed, and undo the move if it does.
    !
    if ( cx(irx)%molcharge(imol) > maxmolcharge ) then
      cx(irx)%molcharge(imol) = cx(irx)%molcharge(imol) - 1
      chargemove(irx,imol) = 0
    else if (cx(irx)%molcharge(imol) < minmolcharge ) then
      cx(irx)%molcharge(imol) = cx(irx)%molcharge(imol) + 1
      chargemove(irx,imol) = 0
    endif

    ! Run some checks on the total charges along the path....
    !
    nn = 0
    do i = 1, nrxn

      ia = 0
      ncharge = 0
      do k = 1, NMOLMAX
        if (chargemove(i,k) /= 0)then
          ia = 1
          ncharge = ncharge + 1
        endif
      enddo
      if (ia /= 0)nn = nn + 1

      TotalCharge = 0
      do j = 1, cx(i)%nmol
        TotalCharge = totalCharge + cx(i)%molcharge(j)
      enddo

      ! Number of charged molecules exceeded...
      if (ncharge > Nchargemol) then
        errflag = .true.
      endif

      ! Total maximum charge exceeded...
      if (abs(TotalCharge) > MaxTotalCharge) then
        errflag = .true.
      endif
    enddo

    ! Total number of charge-changing steps exceeded.
    !
    if (nn > maxstepcharge) then
      errflag = .TRUE.
    endif


    return
  end Subroutine UpdateCharges



  !************************************************************************
  !> RunNetGrow
  !!
  !! Runs automatic growth of a reaction network, based on allowed
  !! graph-moves.
  !!
  !************************************************************************
  !
  Subroutine RunNetGrow( )
    implicit none
    real(8) :: error, error_old, rx, vbest, ir
    real(8) :: x,y,z, gsum, Qvalue, Qvalue2,dE, Qmax,error_temp,bsum
    integer :: ii, imove, irx, itype,k,it,ncount,imcount,nsuccess,ia,nn,ncharge,TotalCharge
    integer :: i, j, imc, ix(NAMOVEMAX), na,isum,isum2,i1,i2,im1(NRXNMAX),im2(NRXNMAX)
    integer :: iend,counter
    integer, allocatable ::  movenum(:), moveatoms(:,:)
    integer, allocatable ::  msmovenum(:,:), msmoveatoms(:,:,:)
    integer, allocatable :: chargemove(:,:), chargemove_store(:,:)
    type (cxs) :: cx_start, cx_end
    type(cxs), allocatable :: cx(:), cxtemp(:)
    logical :: errflag, fail, success
    logical :: printflag,cyc, errx,ChangeCharges
    logical, allocatable :: atomchange(:), bondchange(:,:)
    integer :: rxindex(NAMAX), nrx
    integer, allocatable :: movenum_store(:), moveatoms_store(:,:)
    integer :: itry, idof, irxn
    real(8) :: beta, fac, mctstore, vbe, TotalError, grapherror,maxbarrier
    type(rxp) :: rp
    character (len=7) :: x1
    character (len=8) :: fmt
    character (len=20) :: printfile, file_root
    character (len=15) :: fout


    ! Whatever is in the input file, set igfunc to zero here...it's irrelevant for the
    ! purposes of single-ended network generation:
    !
    igfunc = 0


    ! Open output file mcopt.dat to print graph-error function.
    !
    write(logfile,'("* Running network-generation calculation...."/)')
    call flush(logfile)


    ! Assign the start-point structure from the input files, startfile.
    ! Note that startfile is read from the main input file.
    !
    write(logfile,'(/"* Reading reactant structure...")')
    call flush(logfile)
    Call ReadCXS( cx_start, startfile )
    Call SetMass(cx_start)
    Call GetGraph( cx_start )
    Call Getmols(cx_start)
    Call PrintCXSGraphInfo(cx_start,logfile,"Reactant structure")

    ! Make a copy in cx_end so that the evaluation of graph-error is OK.
    !
    Call CopytoNewCXS(cx_start,cx_end)


    ! Allocate space for atomchange and bondchange arrays - these will indicate at each MC search
    ! step which atoms and bonds are allowed to change.
    !
    na = cx_start%na
    allocate( atomchange(na) )
    allocate( bondchange(na,na) )


    ! Read the graphmoves.
    !
    Call ReadGraphMoves( movefile )


    ! Based on the movefile, decide whether or not we're also going to have to
    ! consider changes in charge states too.
    !
    ChangeCharges = .FALSE.
    do i = 1, ngmove
      if (namove(i) == 0) then
        ChangeCharges = .TRUE.
      endif
    enddo


    ! Create space for the TOTAL reaction-string. This is defined by
    ! nrxn in the input file. Here, nrxn defines the (maximum) number of reactions
    ! which can occur. Here, cx(1) is the structure formed after application of reaction
    ! 1 to the cx_start, cx(2) is the structure formed after application of reaction 2 to
    ! cx(1), and so on....
    !
    allocate( cx(nrxn) )


    ! Copy all aspects of the end-point structures to the intermediates for consistency:
    !
    do i = 1, nrxn
      Call CopytoNewCXS( cx_start, cx(i) )
      Call SetMass(cx(i))
    enddo


    ! Allocate space for the moves:
    !
    allocate( movenum(nrxn) )
    allocate( moveatoms(nrxn,NAMOVEMAX) )
    allocate( movenum_store(nrxn) )
    allocate( moveatoms_store(nrxn,NAMOVEMAX) )
    allocate( chargemove(nrxn,NMOLMAX) )
    allocate( chargemove_store(nrxn,NMOLMAX))


    ! Set all the initial moves to be null moves.
    !
    write(logfile,'("* Setting all initial moves to null...")')
    call flush(logfile)
    do i = 1, nrxn
      movenum(i) = 0
      moveatoms(i,1:NAMOVEMAX) = 0
    enddo


    ! Zero the charges for all molecules along the starting path...
    !
    cx_start%molcharge(:) = 0
    do i = 1, nrxn
      cx(i)%molcharge(:) = 0
      chargemove(i,:) = 0
    enddo


    ! Evaluate the initial error, defined as the graph-distance between the
    ! graph resulting from applying nrxn reactions and the graph for the
    ! end-point:
    !

    Call GetPathFitness( cx_start, cx_end, cx, nrxn, movenum, moveatoms, errflag, &
    GraphError, TotalError, vbe )


    ! If an error flag is returned here, something is wrong with the initial path....
    !
    if (errflag)Stop '* Initial errflag for path is TRUE - something weird going &
    on with initial path....go check it out....'


    ! Output initial errors to logfile.
    !
    write(logfile,'(/"* Initial Graph Error = ",2x,f10.5/)')GraphError
    write(logfile,'("* Initial Total Error = ",2x,f10.5/)')TotalError
    call flush(logfile)


    ! Loop over MC moves, trying to anneal down to zero error:
    !
    write(logfile,'("*** Starting reaction-path search ***"/)'); call flush(logfile)
    mctstore = mcrxntemp
    counter = 0
    outer: do imc = 1, nmcrxn

      ! Make a change to the current reaction-mechanism...
      !
      Call UpdateMechanism(nrxn,movenum,moveatoms,bondchange,atomchange,na,cx_start,cx,rxindex,cyc,&
      movenum_store, moveatoms_store)
      if (cyc) cycle outer


      ! Propagate the graph and evaluate the new error
      !
      Call GetPathFitness( cx_start, cx_end, cx, nrxn, movenum, moveatoms, errflag, &
      GraphError, TotalError, vbe )


      ! Do we change the charges on the molecules along the path?
      !
      if (ChangeCharges) Call UpdateCharges(nrxn,cx,chargemove, chargemove_store, errflag)


      ! If errors detected, replace with last sensible path.
      !
      if (errflag) then
        movenum(:) = movenum_store(:)
        moveatoms(:,:) = moveatoms_store(:,:)
        Totalerror = error_old
        if (changecharges) then
          chargemove = chargemove_store
        endif
      endif


      ! Output the structure along the reaction path - if
      ! optaftermove = .TRUE., we also do geometry optimization.
      !
      fmt = '(I7.7)'
      counter = counter + 1
      write (x1,fmt) counter
      fout = 'rxn_'//trim(x1)//'.xyz'
      Call GraphsToCoords(cx_start, cx, nrxn, .TRUE., fout )


    enddo outer

    stop




    ! At this point, we find the first instance of the target molecule(s)
    ! being formed then set the remainder of movenum to be zero....
    !
    Call TrimPath(cx_start, cx_end, cx, nrxn, movenum, moveatoms, iend )


    ! Print molecules generated along reaction:
    !
    Call PrintMolsAlongPath( nrxn, cx_start, cx_end, cx, movenum,chargemove, &
    changecharges,moveatoms )


    ! At this point, we need to turn graphs into coordinates....
    !
    write(logfile,'("* Generating structures from connectivity matrices..."/)')
    write(logfile,'("* Reaction end-points printed to final_path.xyz"/)')
    Call GraphsToCoords(cx_start, cx, nrxn, .TRUE., 'final_path.xyz' )
    write(logfile,'("* Finished structure generation."/)')


    ! Now we need to do error correction - if the end-point structures have changed due to geometry
    ! optimization, such that the structures no longer represent the zero-error reaction path, we need to
    ! work out which reactions "did" happen, then modify the movenum(:) and moveatoms(:,:). We also need to
    ! check that the reaction-string is still valid.
    !
    ! write(6,*)'OPT = ',optaftermove
    ! if (optaftermove) then
    !    Call ReactionErrorCorrection(cx_start, cx, nrxn, movenum, moveatoms,errx)
    !    if (errx) then
    !      write(logfile,'(/"*** ERROR CORRECTION ATTEMPT FAILED ***"/)')
    !    else
    !      write(logfile,'(/"*** ERROR CORRECTION ATTEMPT SUCCESSFUL ***"/)')
    !      write(logfile,'("\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/")')
    !      write(logfile,'("*** CORRECTED MECHANISM ***")')
    !      Call PrintMolsAlongPath( nrxn, cx_start, cx_end, cx, movenum,chargemove,changecharges,moveatoms )
    !
    !     ! Recalcualte energy.
    !      write(6,*)'INTO AbInitio here...',cx_start%nmol;call flush(6)
    !     Call AbInitio( cx_start, 'ener', success )
    !     do i = 1, nrxn
    !       print*,'INTO AbInitio here...',cx(i)%nmol,cx(i)%molcharge(1:cx(i)%nmol)
    !       Call AbInitio( cx(i), 'ener', success)
    !     enddo
    !
    !   endif
    ! endif

    ! Calculate and output the energy information, including the calculated Q-value:
    !
    if (optaftermove) then
      open(94,file = 'final_path_energy.dat', status='unknown')
      write(94,'("# Reaction index | Absolute energy in au | Relative energy of reaction in kJ/mol")')
      write(94,'("# Energy relative to reactants in kJ/mol")')
      write(logfile,'("===============        ENERGY INFORMATION         ==============="/)')
      write(logfile,'("* Initial energy =",3x,f14.8,2x,"Eh")')cx_start%vcalc
      call flush(logfile)
      write(94,*)'0  ',cx_start%vcalc, '  0.00 ','   0.00 '
      Qvalue = 0.d0
      Qvalue2 = 0.d0
      Qmax = -1d6
      ncount = 0
      do irxn = 1, nrxn
        if (irxn == 1) then
          dE = (cx(irxn)%vcalc - cx_start%vcalc) * au_to_kjmol
        else
          dE = (cx(irxn)%vcalc - cx(irxn-1)%vcalc) * au_to_kjmol
        endif
        write(94,*)irxn,cx(irxn)%vcalc, dE,(cx(irxn)%vcalc-cx_start%vcalc)*au_to_kjmol
        if (abs(de) > Qmax)Qmax = abs(de)
        Qvalue2 = Qvalue2 + abs(dE)
        if (movenum(irxn) /= 0) then
          ncount = ncount + 1
          Qvalue = Qvalue + abs(dE)
        endif
        write(logfile,'("* Reaction #",1x,i4,1x,":",2x,"dE = ",2x,f14.8)')irxn,dE
      enddo
      Qvalue = Qvalue / dble(ncount)
      write(logfile,'(/"* Final Q-value (NOT averaged over non-zero rxns) =",2x,f16.8,2x,"kJ/mol"/)')Qvalue2
      write(logfile,'(/"* Number of non-zero reactions =",2x,i4,2x,"kJ/mol"/)')ncount
      write(logfile,'(/"* Final Q-value (averaged over non-zero rxns) =",2x,f16.8,2x,"kJ/mol"/)')Qvalue
      write(logfile,'(/"* Maximum dE =",2x,f16.8,2x,"kJ/mol"/)')Qmax
      write(logfile,'("================================================================"/)')
      close(94)
    endif

    ! For each reaction along the path, print nimage intermediate snapshots to file.
    ! These snapshots can be used for subsequenct NEB calculations.
    !
    write(logfile,'("* Reaction paths printed to final_path_rx_*.xyz")'); call flush(logfile)
    file_root = 'final_path'
    Call PrintMechanismPaths(nrxn, cx_start, cx, file_root, maxbarrier,bsum,movenum )
    write(logfile,*)"* Maximum barrier = ",maxbarrier * au_to_kjmol,"kJ/mol";call flush(logfile)
    write(logfile,*)"* Barrier sum = ",bsum,"au";call flush(logfile)


    ! If we're running a calculation with igfunc=4, so we're trying to target formation
    ! of a specific molecule (or set of molecules), then we need to only focus on
    ! outputting mechanism and energy data for the set of molecules leading to the target.
    !
    if (igfunc == 4)Call AdjustPaths(cx_start, cx_end, cx, nrxn, movenum, moveatoms,iend )


    ! Clean up memory.
    !
    deallocate( movenum )
    deallocate( moveatoms )
    deallocate( movenum_store )
    deallocate( moveatoms_store )

    41 Format(i7,3x,f16.7,3x,f16.7,3x,f16.7,3x,f16.7)

    return
  end Subroutine RunNetGrow



end Module pathfinder
