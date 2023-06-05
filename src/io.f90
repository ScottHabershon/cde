!
!***************************************************************************************
!
!> @brief Main input/output functions for CDE.
!!
!***************************************************************************************
!

Module IO

  use constants
  use globaldata
  use functions
  implicit none

contains

  !*************************************************************************
  !
  !> ReadInput
  !!
  !! Reads the main CDE input file, and fills in the data in globaldata.f90.
  !!
  !! The input filename is read as the first argument from the command-line.
  !!
  !************************************************************************
  !!
  Subroutine ReadInput()
    implicit none
    integer :: narg, line, pos, ios, i,j,k, imin,imax,id
    integer :: i1,i2,ic
    character (len=2) :: ctemp, a1,a2
    character (len=8) :: c_date
    character (len=10):: c_time
    character(len=100) :: buffer, label, comment
    character (len=40) :: buffer2, buffer3, buffer4, templatefile,label2
    logical :: there
    real(8) :: r1, sum

    ! Get the current date and time.
    !
    call date_and_time(date=c_date, time=c_time)


    ! Read the input filename from the command line.
    !
    narg = iargc()
    if ( narg .ne. 1 ) then
      write(6,'(/"* Program usage: cde.x [input file]"/)')
      stop
    endif
    call getarg (1, inputfile )


    ! Open the input file and read the input parameters.
    !
    inquire( file = inputfile, EXIST=THERE )
    if (.not.there)stop '* Specified input file does not exist'
    Open(10, file = inputfile, status = 'unknown')


    ! Open a *.log file with the same root as the input file.
    !
    Open(logfile, file = trim( inputfile )//'.log', status = 'unknown')
    write(logfile,'(/"************************ CDE simulation ************************"/)')
    write(logfile,'("- Start date: ",3x,a2,"/",a2,"/",a4)')c_date(7:8),c_date(5:6),c_date(1:4)
    write(logfile,'("- Start time: ",3x,a2,":",a2,":",a2/)')c_time(1:2),c_time(3:4),c_time(5:6)
    write(logfile,'("- The input filename is:",1x,a/)')inputfile
    write(logfile,'("================================================================")')
    write(logfile,'("                INPUT FILE INTERPRETATION FOLLOWS               ")')
    write(logfile,'("================================================================"/)')
    write(logfile,'("- The input filename is:",1x,a/)')inputfile


    ! Set default values for everything.
    !
    Call SetIODefaults()


    ! Loop over lines in the file, reading keywords and arguments.
    !
    line = 0
    ios = 0
    do while (ios == 0)
      read(10, '(A)', iostat=ios) buffer
      if (ios == 0) then
        line = line + 1

        ! Find the first space
        pos = scan(buffer, ' ')
        label = buffer(1:pos)
        buffer = buffer(pos+1:)


        ! If pos == 1, there is a space in the first position and
        ! we assume a blank line.
        if (pos == 1)then
          write(logfile,*)
          CYCLE
        endif

        ! Read input based on label.
        !
        select case (label)

        case ('# ')
          read(buffer,*,iostat=ios) comment
          write(logfile,'("# COMMENT:",1x,a)')buffer

          !***********************************************************************
          ! GENERAL INPUT PARAMETERS.
          !***********************************************************************
          !
        case ('nimage')       ! Number of images along paths.
          read(buffer,*,iostat=ios) nimage
          write(logfile,'("- Number of images per path: ",1x,i4)')nimage
          if (nimage <=0 )stop '* ERROR: nimage <= 0'


        case ('calctype')
          ! 'optpath' for path optimization,
          ! 'molopt' for molecular property optimization (EXPERIMENTAL),
          ! 'pathfind' for path-finding calculation,
          ! 'evb' for EVB-based initial MEP generation (EXPERIMENTAL).
          read(buffer, *, iostat=ios) calctype
          if (calctype == 'optpath') then
            write(logfile,'("- Calculation type:  PATH REFINEMENT")')
          else if (calctype == 'interp') then
            write(logfile, '("- Calculation type: PATH INTERPOLATION")')
          else if (calctype == 'molopt') then
            write(logfile,'("- Calculation type:  MOLECULE PROPERTY OPTIMIZATION")')
          else if (calctype == 'pathfind') then
            write(logfile,'("- Calculation type:  REACTION-PATH FINDER")')
          else if (calctype == 'netgrow' .or. calctype == 'netgrow2') then
            write(logfile,'("- Calculation type:  NETWORK GENERATION")')
          else if (calctype == 'breakdown') then
            write(logfile, '("Calculation type: BREAKDOWN NETWORK GENERATION")')
          else
            stop 'ERROR: Unknown calctype in input file (options: optpath, interp, molopt, &
                & pathfind, netgrow, breakdown)'
          endif


        case ('pathinit')  ! Path initialization method - generally has to be 'linear', not
          ! really used other than to give internal images position values
          read(buffer, *, iostat=ios) pathinit
          if (pathinit == 'linear') then
            write(logfile,'("- Initial path initialization type:  LINEAR")')
          else
            stop 'ERROR: Unknown PathInit in input file (options: linear)'
          endif


        case ('startfile')  ! For some calculations, reactants and product structures are required.
          ! startfile is the xyz file containing the reactant structure....
          read(buffer, *, iostat=ios) startfile
          write(logfile,'("- Start-point structure file: ",1x,a)')startfile

        case ('endfile')    ! ..and endfile is the xyz file containing the product structure.
          read(buffer, *, iostat=ios) endfile
          write(logfile,'("- End-point structure file: ",1x,a)')endfile

        case ('pathfile') ! xyz file containing an entire reaction path - number of images here
          ! must be same as nimage in input file.
          read(buffer, *, iostat=ios) pathfile
          write(logfile,'("- Starting path file: ",1x,a)')pathfile

        case ('startfrompath')
          read(buffer, *, iostat = ios) startfrompath
          write(logfile,'("- Start from path?: ",1x,L)')startfrompath

        case ('ranseed')
          read(buffer,*,iostat=ios)irun
          write(logfile,'("- Random number seed: ",1x,i4)')irun


          !***********************************************************************
          ! EVB CALCULATION INPUT PARAMETERS.
          !***********************************************************************
          !
        case ('evbtype')   ! EVB-based MEP generation type (integer)
          read(buffer,*,iostat=ios)evbtype
          if (evbtype == 1) then
            write(logfile,'("- EVB type 1: HARMONIC EVB METHOD")')
          else if (evbtype == 2) then
            write(logfile,'("- EVB type 2: DISTANCE-BASED EVB METHOD")')
          else if (evbtype == 3) then
            write(logfile,'("- EVB type 3: PES-BASED EVB METHOD")')
          else
            stop 'ERROR: Unknown evbtype in input file '
          endif

        case ('evbiter')  ! Number of iterations for optimization during EVB.
          read(buffer,*,iostat=ios)evbiter
          if (evbiter <= 0)stop 'ERROR: evbiter not acceptable'
          write(logfile,'("- EVB iterations: ",1x,i7)')evbiter

        case ('nevbqm')  ! Number of iterations for quickmin optimization after EVB.
          read(buffer,*,iostat=ios)nevbqm
          if (evbiter < 0)stop 'ERROR: nevbqm not acceptable'
          write(logfile,'("- EVB quick-min iterations: ",1x,i7)')nevbqm

        case ('nevbl')  ! Number of intermediate structures in EVB.
          read(buffer,*,iostat=ios)nevbl
          if (nevbl <= 0)stop 'ERROR: nevbl not acceptable'
          write(logfile,'("- Number of EVB intermediates: ",1x,i7)')nevbl

        case ('nevbdist')  ! Number of distribution tests in EVB.
          read(buffer,*,iostat=ios)nevbdist
          if (nevbdist <= 0)stop 'ERROR: nevbdist not acceptable'
          write(logfile,'("- Number of EVB intermediate tests: ",1x,i7)')nevbdist

        case ('evbstep') ! EVB steepest-descent step length.
          read(buffer,*,iostat=ios)evbstep
          if (evbstep < 0.d0)stop 'ERROR: evbstep < 0'
          write(logfile,'("- EVB gradient step: ",1x,f14.8)')evbstep

        case ('evbmaxdl') ! Maximum target distance between adjacent
          ! images in EVB path
          read(buffer,*,iostat=ios)evbmaxdl
          if (evbmaxdl <= 0.d0 )stop 'ERROR: evbmaxdl <= 0'
          write(logfile,'("- EVB maximum image-distance: ",1x,f14.8)')evbmaxdl

        case ('evbalpha1')  ! \alpha_1 parameter for EVB MEP
          read(buffer,*,iostat=ios)evbalpha1
          if (evbalpha1 < 0.d0) stop 'ERROR: evbalpha1 < 0'
          write(logfile,'("- EVB alpha1: ",1x,f14.8)')evbalpha1

        case ('evbalpha2') ! \alpha_2 parameter for EVB MEP
          read(buffer,*,iostat=ios)evbalpha2
          if (evbalpha1 < 0.d0) stop 'ERROR: evbalpha2 < 0'
          write(logfile,'("- EVB alpha2: ",1x,f14.8)')evbalpha2

        case ('evbvrep')       ! Logical flag indicating whether to use repulsive
          ! potential during EVB calculation.
          read(buffer,*,iostat=ios)evbvrep
          if (evbvrep) then
            write(logfile,'("- Use Vrep in EVB: TRUE")')
          else
            write(logfile,'("- Use Vrep in EVB: FALSE")')
          endif

        case ('alphavbe')   ! Strength parameter for adding-in approximate additive
          ! bond-strength changes during double-ended pathfinding.
          read(buffer,*,iostat=ios)alphavbe
          write(logfile,'("- Strength parameter for additive bonds in &
          double-ended GDS: ",1x,f14.8)')alphavbe


          !***********************************************************************
          ! MEP CALCULATION INPUT PARAMETERS.
          !***********************************************************************
          !
        case ('pathoptmethod')  ! MEP refinement method - 'cineb' for (climbing image) nudged
          ! elastic band, or 'zts' for zero-temperature string (experimental)
          read(buffer, *, iostat=ios) pathoptmethod
          if (pathoptmethod == 'cineb') then
            write(logfile,'("- Path refinement method:  CLIMBING-IMAGE NUDGED ELASTIC BAND")')
          else if (pathoptmethod == 'zts') then
            write(logfile,'("- Path refinement method:  ZERO-TEMPERATURE STRING")')
          else
            stop 'ERROR: Unknown pathoptmethod in input file (options: cineb or zts)'
          endif

        case ('nebmethod')  ! Method use for optimization in NEB calculations.
          ! pathquad - quadratic interpolation.
          ! steepest - steepest descents method
          ! quickmin - The Quick-min method (recommended)
          ! fire - Fast inertial relaxation engine.
          ! dfp - BFGS-type scheme, currently unstable!
          ! bmk - Runge-Kutta method
          ! hball - Heavy rolling ball method.
          read(buffer, *, iostat=ios) NEBmethod
          if (NEBmethod == 'pathquad') then
            write(logfile,'("- (CI)NEB optimization method:  QUADRATIC INTERPOLATION")')
          else if (NEBmethod == 'steepest') then
            write(logfile,'("- (CI)NEB optimization method:  STEEPEST DESCENT")')
          else if (NEBmethod == 'quickmin') then
            write(logfile,'("- (CI)NEB optimization method:  QUICKMIN")')
          else if (NEBmethod == 'fire') then
            write(logfile,'("- (CI)NEB optimization method:  FIRE")')
          else if (NEBmethod == 'dfp') then
            write(logfile,'("- (CI)NEB optimization method:  DFP (BFGS)")')
          else if (NEBmethod == 'bmk') then
            write(logfile,'("- (CI)NEB optimization method:  BMk (Runge-Kutta)")')
          else if (NEBmethod == 'hball') then
            write(logfile,'("- (CI)NEB optimization method:  HEAVY-BALL ")')
          else
            stop 'ERROR: UNKNown NEBmethod in input file (options: pathquad, steepest, &
            & quickmin, fire, dfp, bmk, hball)'
          endif

        case ('nebconv')        ! NEB convergence criteria on RMS forces, in au.
          read(buffer,*,iostat=ios) NEBconv
          if (NEBconv <= 0.d0) then
            stop '* ERROR: NEBconv <= 0.0'
          else
            write(logfile,'("- NEB NORM force-convergence criterion (au): ",1x,f14.8)')NEBconv
          endif

        case ('nebmaxconv')   ! Convergence criterion on maximum force during NEB, in au.
          read(buffer,*,iostat=ios) NEBMAXconv
          if (NEBmaxconv <= 0.d0) then
            stop '* ERROR: NEBmaxconv <= 0.0'
          else
            write(logfile,'("- NEB MAXIMUM force-convergence criterion (au): ",1x,f14.8)')NEBmaxconv
          endif

        case ('cithresh')   ! Threshold on RMS forces at which climbing-image variant is engaged.
          read(buffer,*,iostat=ios) CIthresh
          if (CIthresh <= 0.d0) then
            stop '* ERROR: CIthresh <= 0.0'
          else
            write(logfile,'("- Climbing-image start criterion (au): ",1x,f14.8)')CIthresh
          endif

        case ('stripinactive') ! Remove inactive molecules from path.
          read(buffer, *, iostat = ios) stripinactive
          write(logfile,'("- Strip inactive molecules before path optimization?: ",1x,L)')stripinactive

        case ('alignedatoms')
          read(10, *, iostat = ios) atomidx(1),atomidx(2),atomidx(3)
          write(logfile,'("- Aligned atoms for NEB refinements: ",1x,3i5)')atomidx(1:3)

        case ('nebiter')
          read(buffer,*,iostat=ios)NEBiter
          write(logfile,'("- Maximum number of NEB iteractions: ",1x,i6)')NEBiter

        case ('neboutfreq')
          read(buffer,*,iostat=ios)NEBoutfreq
          write(logfile,'("- NEB output frequency: ",1x,i6)')NEBoutfreq

        case ('nebspring')
          read(buffer,*,iostat=ios)NEBspring
          write(logfile,'("- NEB spring constant: ",1x,f14.6)')NEBspring

        case ('nebstep')
          read(buffer,*,iostat=ios)NEBstep
          write(logfile,'("- NEB optimization step-size: ",1x,f14.6)')NEBstep

        case ('optendsbefore')
          read(buffer,*,iostat = ios)optendsbefore
          write(logfile,'("- Optimize end-points before NEB run: ",1x,L)')optendsbefore

        case ('reconnect')
          read(buffer,*,iostat=ios)reconnect
          write(logfile,'("- Reconnect end-points after optimization: ",1x,L)')reconnect

        case ('optendsduring')
          read(buffer,*,iostat = ios)optendsduring
          write(logfile,'("- Optimize end-points during NEB run: ",1x,L)')optendsduring

        case ('idppguess')
          read(buffer, *, iostat=ios) idpppath
          write(logfile,'("- IDPP interpolation : ",1x,L)')idpppath


        case ('projforcetype')
          read(buffer, *, iostat=ios) projforcetype
          write(logfile,'("- Using NEB projected forces of type: ",1x,I2)')projforcetype

        case ('nebrestrend')
          read(buffer,*,iostat = ios)nebrestrend
          write(logfile,'("- Use graph-restrain potential during NEB: ",1x,L)')nebrestrend

        case ('vsthresh')
          read(buffer, *, iostat = ios) VSthresh
          write(logfile,'("- Make variable spring strengths with threshold?: ",1x,f8.6)')VSthresh


          !***********************************************************************
          ! GDS PARAMETERS.
          !***********************************************************************
          !
        case('reactiveatomtypes{')
          nreactivetypes = 0
          do i = 1, NTYPEMAX
            read(10,*,iostat=ios) ctemp
            if (ctemp /= '}') then
              nreactivetypes = nreactivetypes + 1
              reactivetype(nreactivetypes) = ctemp
            else if (ctemp == '}') then
              EXIT
            endif
          enddo
          write(logfile,'("- Number of reactive atom types: ",1x,i4)')nreactivetypes
          do i = 1, nreactivetypes
            write(logfile,'("- Reactive atom labels: ",1x,A)')reactivetype(i)
          enddo

        case('valencerange{')
          nvalcon = 0
          do i = 1, NTYPEMAX
            read(10,'(A)',iostat=ios)buffer2
            if (buffer2(1:1) /= '}') then
              a2 = ''
              if (index(buffer2,'fz').ne. 0 ) &
              read(buffer2,*)a1,imin,imax, a2
              if (index(buffer2,'fz').eq. 0 ) &
              read(buffer2,*)a1,imin,imax
              nvalcon = nvalcon + 1
              valatom(nvalcon) = a1
              valrange(nvalcon,1) = imin
              valrange(nvalcon,2) = imax
              if (trim(a2) .eq. 'fz' ) valfz(nvalcon) = .true.
            else if (buffer2(1:1) == '}') then
              EXIT
            endif
          enddo
          write(logfile,'("- Number of valence range constraints: ",1x,i4)')nvalcon
          do i = 1, nvalcon
            write(logfile,'("- Valence range constraint: ",1x,i3,1x,A,1x,i4,1x,i4)')i,valatom(i),valrange(i,1),valrange(i,2)
          enddo

        case('reactivevalence{')
          nrxval = 0
          do i = 1, NTYPEMAX
            read(10,'(A)',iostat=ios)buffer2
            if (buffer2(1:1) /= '}') then
              read(buffer2,*)a1,a2,imin,imax
              nrxval = nrxval + 1
              rxvalatom(nrxval,1) = a1
              rxvalatom(nrxval,2) = a2
              rxvalrange(nrxval,1) = imin
              rxvalrange(nrxval,2) = imax
            else if (buffer2(1:1) == '}') then
              EXIT
            endif
          enddo
          write(logfile,'("- Number of reactive valence constraints: ",1x,i4)')nrxval
          do i = 1, nrxval
            write(logfile,'("- Reactive valence range constraint: ",1x,i3,1x,A,1x,A,1x,i4,1x,i4)') &
            i,rxvalatom(i,1),rxvalatom(i,2),rxvalrange(i,1),rxvalrange(i,2)
          enddo

        case('reactiveatoms{')
          reactive(:) = .FALSE.
          do i = 1, NAMAX
            read(10,'(A)',iostat=ios) buffer2
            pos = scan(buffer2, ' ')
            label2 = buffer2(1:pos)
            buffer2 = buffer2(pos+1:)
            select case(label2)
            case('range')
              read(buffer2,*,iostat=ios)imin,imax
              write(logfile,'("- Reactive atom range: ",1x,i4," - ",i4)')imin,imax
              do j = imin, imax
                reactive(j) = .TRUE.
              enddo
            case('id')
              read(buffer2,*,iostat=ios)id
              reactive(id) = .TRUE.
            case('all')
              if (.not. startfrompath) &
              open(12, file = startfile, status = 'unknown')
              if ( startfrompath) &
              open(12, file = pathfile, status = 'unknown')
              read(12, *, iostat=ios) imax
              do j = 1, imax
                reactive(j) = .TRUE.
              enddo
              close(12)
            case('}')
              EXIT
            case default
              print*,'* Error in reactiveatoms{} block'
              stop
            end select
          enddo

        case('essentialmoveatoms{')
          essentialmoves(:) = .FALSE.
          do i = 1, NAMAX
            read(10,'(A)',iostat=ios) buffer2
            pos = scan(buffer2, ' ')
            label2 = buffer2(1:pos)
            buffer2 = buffer2(pos+1:)
            select case(label2)
            case('range')
              read(buffer2,*,iostat=ios)imin,imax
              write(logfile,'("- Essential atom range: ",1x,i4," - ",i4)')imin,imax
              do j = imin, imax
                essentialmoves(j) = .TRUE.
              enddo
            case('id')
              read(buffer2,*,iostat=ios)id
              essentialmoves(id) = .TRUE.
            case('all')
              open(12, file = startfile, status = 'unknown')
              read(12, *, iostat=ios) imax
              do j = 1, imax
                essentialmoves(j) = .TRUE.
              enddo
              close(12)
            case('}')
              EXIT
            case default
              print*,'* Error in essentialatoms{} block'
              stop
            end select
          enddo
          ! Count number of essential moves...
          !
          nessentialatoms = 0
          do i = 1, NAMAX
            if (essentialmoves(i))nessentialatoms = nessentialatoms + 1
          enddo

        case('essentialatomsinmols{')
          essentialatomsinmols(:) = 0
          ic = 0
          do i = 1, NAMAX
            read(10,'(A)',iostat=ios) buffer2
            pos = scan(buffer2, ' ')
            label2 = buffer2(1:pos)
            buffer2 = buffer2(pos+1:)
            select case(label2)
            case('range')
              read(buffer2,*,iostat=ios)imin,imax
              write(logfile,'("- Essential atoms in mols range: ",1x,i4," - ",i4)')imin,imax
              do j = imin, imax
                ic = ic + 1
                essentialatomsinmols(ic) = j
              enddo
            case('id')
              read(buffer2,*,iostat=ios)id
              ic = ic + 1
              essentialatomsinmols(ic) = id
            case('}')
              EXIT
            case default
              print*,'* Error in essentialatoms{} block'
              stop
            end select
          enddo

          ! Count number of essential moves...
          !
          nessentialatomsinmols = ic

        case('essentialatoms{')
          essential(:) = .FALSE.
          do i = 1, NAMAX
            read(10,'(A)',iostat=ios) buffer2
            pos = scan(buffer2, ' ')
            label2 = buffer2(1:pos)
            buffer2 = buffer2(pos+1:)
            select case(label2)
            case('range')
              read(buffer2,*,iostat=ios)imin,imax
              write(logfile,'("- Essential atom range: ",1x,i4," - ",i4)')imin,imax
              do j = imin, imax
                essential(j) = .TRUE.
              enddo
            case('id')
              read(buffer2,*,iostat=ios)id
              essential(id) = .TRUE.
            case('all')
              open(12, file = startfile, status = 'unknown')
              read(12, *, iostat=ios) imax
              do j = 1, imax
                essential(j) = .TRUE.
              enddo
              close(12)
            case('}')
              EXIT
            case default
              print*,'* Error in essentialatoms{} block'
              stop
            end select
          enddo

        case('allowedbonds{')
          nallowbonds = 0
          do i = 1, NTYPEMAX
            read(10,'(A)',iostat=ios)buffer2
            if (buffer2(1:1) /= '}') then
              read(buffer2,*)a1,a2,imin
              nallowbonds = nallowbonds + 1
              allowbondsatom(nallowbonds,1) = a1
              allowbondsatom(nallowbonds,2) = a2
              allowbondsmax(nallowbonds) = imin
            else if (buffer2(1:1) == '}') then
              EXIT
            endif
          enddo
          write(logfile,'("- Number of allowed bonding constraints: ",1x,i4)')nallowbonds
          do i = 1, nallowbonds
            write(logfile,'("- Allowed bonding constraints: ",1x,i3,1x,A,1x,A,1x,i4,1x,i4)') &
            i,allowbondsatom(i,1),allowbondsatom(i,2),allowbondsmax(i)
          enddo

        case('fixedbonds{')
          write(logfile,'("- Fixed bonds for GDS:")')
          fixedbonds(:,:) = .FALSE.
          nfixtype = 0
          do i = 1, NAMAX
            read(10,'(A)',iostat=ios)buffer2
            if (buffer2(1:1) == '}') then
              EXIT
            else if ( .not.is_numeric(buffer2(1:1)) ) then
              read(buffer2,*,iostat=ios)a1,a2
              write(logfile,'("- Bonds are fixed for atom types: ",1x,a2," - ",a2)')a1,a2
              nfixtype = nfixtype + 1
              fixedbondtype(nfixtype,1) = a1
              fixedbondtype(nfixtype,2) = a2
            else
              read(buffer2,*,iostat=ios)i1,i2
              fixedbonds(i1,i2) = .TRUE.
              fixedbonds(i2,i1) = .TRUE.
              write(logfile,'("- Bonds are fixed for atom numbers: ",1x,i4," - ",i4)')i1,i2
            endif
          enddo

        case('gdsspring')
          read(buffer,*,iostat=ios)gdsspring
          write(logfile,'("- Spring constant for inter-bead springs: ",1x,f14.8)')gdsspring

        case('gdsrestspring')
          read(buffer,*,iostat=ios)gdsrestspring
          write(logfile,'("- Spring constant for bonded-graph restraints: ",1x,f14.8)')gdsrestspring

        case('nbstrength')
          read(buffer,*,iostat=ios)nbstrength
          write(logfile,'("- Strength of non-bonded exponential graph restraints: ",1x,f14.8)')nbstrength

        case('nbrange')
          read(buffer,*,iostat=ios)nbrange
          write(logfile,'("- Range of non-bonded exponential graph restraints: ",1x,f14.8)')nbrange

        case('kradius')
          read(buffer,*,iostat=ios)kradius
          write(logfile,'("- Spring constant for intermolecular repulsion restraints: ",1x,f14.8)')kradius


        case ('movefile')
          read(buffer, *, iostat=ios) movefile
          write(logfile,'("- Graph-moves file: ",1x,a)')movefile

        case ('forbidgraphs')
          read(buffer, *, iostat=ios) forbidgraphs
          write(logfile,'("- Monitor forbidden graph-moves file: ",1x,L)')forbidgraphs

        case ('forbidfile')
          read(buffer, *, iostat=ios) forbidfile
          write(logfile,'("- Forbidden-moves file: ",1x,a)')forbidfile

        case ('gdsoutfreq')
          read(buffer,*,iostat=ios) gdsoutfreq
          write(logfile,'("- GDS output frequency: ",1x,i7)') gdsoutfreq

        case ('gdsthresh')
          read(buffer,*,iostat=ios) gdsthresh
          write(logfile,'("- GDS attempt probability: ",1x,f14.8)') gdsthresh

        case ('ngdsrelax')
          read(buffer,*,iostat=ios) ngdsrelax
          write(logfile,'("- Number of SD steps after GDS graph move: ",1x,i7)') ngdsrelax

        case ('gdsdtrelax')
          read(buffer,*,iostat=ios) gdsdtrelax
          write(logfile,'("- GDS SD relaxation step-size: ",1x,f14.8)') gdsdtrelax

        case ('gds_intra_cutoff')
          read(buffer,*,iostat=ios) intra_cutoff
          write(logfile,'("- GDS intramolecular atom move cutoff distance: ",1x,f14.8)') intra_cutoff

        case ('optaftermove')
          read(buffer, *, iostat = ios) optaftermove
          if (.not. optaftermove .and. lmoldata) then
            skiprepeats = .true.
            write(logfile,'("- Moldata has been set to true, so optaftermove is overridden ")')
          endif
          write(logfile,'("- PES geometry optimization after graph move?: ",1x,L)')optaftermove

        case ('ignoreinvalidgraphopt')
          read(buffer, *, iostat = ios) ignoreinvalidgraphopt
          write(logfile, '("- Ignore geometry optimisations that invalidate molecular graphs? ", &
              & 1x, L)') ignoreinvalidgraphopt

        case ('doinitialopt')
          read(buffer,*,iostat=ios) doInitialOpt
          write(logfile,'("- Perform initial geometry optimisation on starting system: ", 1x, L)') &
              & doInitialOpt


          !***********************************************************************
          ! CONSTRAINT PARAMETERS.
          !***********************************************************************
          !
        case ('dofconstraints')
          read(buffer,*,iostat=ios) NDOFconstr
          write(logfile,'("- Number of fixed DOFs: ",1x,i4)')NDOFconstr
          if (NDOFconstr > 0) then
            read(10,*,iostat=ios) (FixedDOF(i),i=1,NDOFconstr)
            do i = 1, NDOFconstr
              write(logfile,'("- DOF",1x,i4," is fixed! ")')FixedDOF(i)
            enddo
          endif

        case ('atomconstraints')
          read(buffer,*,iostat=ios) Natomconstr
          write(logfile,'("- Number of fixed atoms: ",1x,i4)')Natomconstr
          if (Natomconstr > 0) then
            read(10,*,iostat=ios) (FixedAtom(i), i=1,Natomconstr)
            do i = 1, Natomconstr
              write(logfile,'("- Atom",1x,i4," is fixed! ")')FixedAtom(i)
            enddo
          endif


          !***********************************************************************
          ! PES CALCULATION INPUT PARAMETERS.
          !***********************************************************************
          !
        case('pesfull')
          read(buffer,*,iostat=ios)PESfull
          if (PESfull) then
            write(logfile,'("- PES evaluations for FULL system being used...")')
          else
            write(logfile,'("- PES evaluations for INDEPENDENT MOLECULES being used...")')
          endif

        case('pestype')
          read(buffer,*,iostat=ios)PEStype
          write(logfile,'("- PES type: ",1x,A)')trim(PEStype)

        case('pesfile')
          read(buffer,*,iostat=ios)PESfile
          write(logfile,'("- PES calculation template file: ",1x,A)')trim(PESfile)

        case('pesexecutable')
          read(buffer,'(A100)',iostat=ios)PESexecutable
          if (index(PESexecutable,'#') .ne. 0 ) &
          PESexecutable = PESexecutable(:index(PESexecutable,'#')-1)
          write(logfile,'("- PES calculation executable: ",1x,A)')trim(PESexecutable)

        case('pesoptexecutable')
          read(buffer,'(A100)',iostat=ios)PESOptexecutable
          if (index(PESOptexecutable,'#') .ne. 0 ) &
          PESOptexecutable = PESOptexecutable(:index(PESOptexecutable,'#')-1)
          write(logfile,'("- PES geometry optimization executable: ",1x,A)')trim(PESOptexecutable)

        case('pesopttype')
          read(buffer,*,iostat=ios)PESOpttype
          write(logfile,'("- PES type for geometry optimization: ",1x,A)')trim(PESOpttype)

        case('pesoptfile')
          read(buffer,*,iostat=ios)PESOptfile
          write(logfile,'("- Geometry optimisation template file: ",1x,A)')trim(PESOptfile)


          !***********************************************************************
          ! PATHFINDING PARAMETERS.
          !***********************************************************************
          !
        case('nrxn')
          read(buffer,*,iostat=ios)nrxn
          write(logfile,'("- Number of reactions in sampled strings for pathfinder: ",1x,i4)')nrxn

        case('nmcrxn')
          read(buffer,*,iostat=ios)nmcrxn
          write(logfile,'("- Number of MC steps in pathfinder calculation:: ",1x,i7)')nmcrxn


        case('nmechmove')
          read(buffer,*,iostat=ios)nmechmove
          write(logfile,'("- Maximum number of mechanism moves in pathfinder calculation:: ",1x,i7)')nmechmove

        case('mcrxntemp')
          read(buffer,*,iostat=ios)mcrxntemp
          write(logfile,'("- Initial temperature for pathfinder (in K): ",1x,e14.4)')mcrxntemp

        case('graphfunctype')
          read(buffer,*,iostat=ios)igfunc
          if (igfunc == 0) then
            write(logfile,'("- Pathfinding function type: Standard element-wise comparison")')
          else if (igfunc == 1) then
            write(logfile,'("- Pathfinding function type: Eigenvalue comparison 1")')
          else if (igfunc == 2) then
            write(logfile,'("- Pathfinding function type: Histogram-based comparison")')
          else if (igfunc == 3) then
            write(logfile,'("- Pathfinding function type: Eigenvalue comparison 2")')
          else if (igfunc == 4) then
            write(logfile,'("- Pathfinding function type: Eigenvalue comparison for single molecule")')
          else
            stop '* ERROR: Unknown graphfunctype in input file. Must be either 0 (element-wise) or 1 (permutational)'
          endif


          ! min total charge on any one molecule.
        case('minmolcharge')
          read(buffer,*,iostat=ios)minmolcharge
          write(logfile,'("- Minimum allowed molecular charge: ",1x,i4)')minmolcharge

          ! Max total charge on any molecule
        case('maxmolcharge')
          read(buffer,*,iostat=ios)maxmolcharge
          write(logfile,'("- Maximum allowed molecular charge: ",1x,i4)')maxmolcharge

          ! Max total charge of the system at any reaction step.
          !
        case('nchargemol')
          read(buffer,*,iostat=ios)nchargemol
          write(logfile,'("- Maximum number of charged molecules allowed at any step: ",1x,i4)')nchargemol

          ! Max number of reation steps which can involve charge changes.
        case('maxstepcharge')
          read(buffer,*,iostat=ios)maxmolcharge
          write(logfile,'("- Maximum allowed molecular charge: ",1x,i4)')maxstepcharge

          ! Max total charge in a given reaction step.
        case('maxtotalcharge')
          read(buffer,*,iostat=ios)maxtotalcharge
          write(logfile,'("- Maximum total system charge: ",1x,i4)')maxtotalcharge

        case default
          print*
          print*, '* Error at line: ', line
          print*, '* Unknown input: ', label
          print*
          stop
        end select
      endif
    enddo
    close(10)
    write(logfile,'(/"================================================================")')
    write(logfile,'("                 END OF INPUT FILE INTERPRETATION                   ")')
    write(logfile,'("================================================================")')

    ! Wrap it up....
    !
    write(logfile,'(/"*************** Finished reading input parameters **************"/)')
    call flush(logfile)

    return
  end Subroutine ReadInput


  !
  !******************************************************************
  !> ReadGraphMoves
  !!
  !! Reads the set of graph moves used during GDS simulation from
  !! the specified file 'movefile'.
  !!
  !! - movefile: The file containing the graph-move information.
  !!
  !******************************************************************
  !
  Subroutine ReadGraphMoves( movefile )
    implicit none
    character(len=50) :: movefile
    logical :: there
    integer :: line, ios, pos, i, pos2
    character(len=100) :: buffer, label, comment, buffer2, label2, gwritefmt
    real(8) :: sum

    write(logfile,'(/"* Reading graph moves file...")')
    write(logfile,'("* Graph moves input file:",1x,A/)')trim(movefile)

    ! Check and open movefile.
    !
    inquire( file = movefile, EXIST=THERE )
    if (.not.there)stop '* Specified move file does not exist'
    Open(20, file = movefile, status = 'unknown')
    write(logfile,'("================================================================")')
    write(logfile,'("              MOVE FILE INTERPRETATION FOLLOWS...               ")')
    write(logfile,'("================================================================"/)')


    ! Loop over lines in the file, reading keywords and arguments.
    !
    ngmove = 0
    line = 0
    ios = 0
    do while (ios == 0)
      read(20, '(A)', iostat=ios) buffer
      if (ios == 0) then
        line = line + 1

        ! Find the first space
        pos = scan(buffer, ' ')
        label = buffer(1:pos)
        buffer = buffer(pos+1:)

        ! If pos == 1, there is a space in the first position and we assume a blank line.
        if (pos == 1) then
          ! write(logfile,*)label,buffer
          CYCLE
        endif

        ! Read input based on label.
        !
        select case (label)

        ! Comment line
        case ('# ')
          read(buffer, *, iostat=ios) comment
          write(logfile, '("Comment line: ")')
          write(logfile, '(A, " ", A, /)') trim(label), trim(buffer)

          ! Graph move
        case ('move')
          ngmove = ngmove + 1
          write(logfile, '("*** Graph move number:", 1x, i4, " ***"/)') ngmove

          ! The first entry is the number of atoms.
          read(20, '(A)', iostat=ios) buffer2
          pos2 = scan(buffer2 , ' ')
          label2 = buffer2(1:pos2)
          buffer2 = buffer2(pos2+1:)
          read(buffer2, *, iostat=ios) namove(ngmove)
          write(logfile, '("- Number of atoms:", 1x, i4)') namove(ngmove)

          ! Is this a bond-changing move....?
          if (namove(ngmove) > 0) then
            write(logfile, '("- Bond-changing move detected ")')
            write(gwritefmt, '("(", I1, "I3)")') namove(ngmove)

            ! Next, read the starting and ending graph - the first line is a '-', followed
            ! by the namove x namove matrix on the next namove lines.
            write(logfile, '("- Start graph: ")')
            read(20, '(A)', iostat=ios) buffer2
            do i = 1, namove(ngmove)
              read(20, *, iostat=ios) gmstart(ngmove, i, 1:namove(ngmove))
              write(logfile, gwritefmt) gmstart(ngmove, i, 1:namove(ngmove))
            enddo
            read(20, '(A)', iostat=ios) buffer2
            write(logfile, '("- End graph: ")')
            do i = 1, namove(ngmove)
              read(20, *, iostat=ios) gmend(ngmove, i, 1:namove(ngmove))
              write(logfile, gwritefmt) gmend(ngmove, i, 1:namove(ngmove))
            enddo
            read(20, '(A)', iostat=ios) buffer2

            ! Now read the labels of the atoms, if any.
            read(20, '(A)', iostat=ios) buffer2
            pos2 = scan(buffer2, ' ')
            label2 = buffer2(1:pos2)
            buffer2 = buffer2(pos2+1:)
            read(buffer2, *, iostat=ios) (movelabel(ngmove, i), i=1, namove(ngmove))
            write(logfile, '("- Allowed atom labels: ")')
            write(gwritefmt, '("(", I1, "A2)")') namove(ngmove)
            write(logfile, gwritefmt) movelabel(ngmove, 1:namove(ngmove))

            ! Now read the move probability.
            read(20, '(A)', iostat=ios) buffer2
            pos2 = scan(buffer2, ' ')
            label2 = buffer2(1:pos2)
            buffer2 = buffer2(pos2+1:)
            read(buffer2, *, iostat = ios) moveprob(ngmove)
            write(logfile, '("Move probability: ", F5.3/)') moveprob(ngmove)
            call flush(logfile)

          ! ..or is it a molecular charge-changing move?
          else if (namove(ngmove) == 0) then

            write(logfile, '("- Molecular charge-changing move detected ")')
            read(20, '(A)', iostat=ios) buffer2   ! Reading the '-'...
            read(20, '(A)', iostat=ios) buffer2
            pos2 = scan(buffer2, ' ')
            label2 = buffer2(1:pos2)
            movetype(ngmove) = label2
            read(20, '(A)', iostat=ios) buffer2   ! Reading the '-'...

            ! Now read the move probability.
            read(20, '(A)', iostat=ios) buffer2
            pos2 = scan(buffer2, ' ')
            label2 = buffer2(1:pos2)
            buffer2 = buffer2(pos2+1:)
            read(buffer2, *, iostat = ios) moveprob(ngmove)

          endif

        end select

      endif
    enddo

    ! Make sure to normalize probabilities of graph moves.
    sum = 0.d0
    do i = 1, ngmove
      sum = sum + moveprob(i)
    enddo
    do i = 1, ngmove
      moveprob(i) = moveprob(i) / sum
    enddo
    write(logfile,'("================================================================")')
    write(logfile,'("                MOVE FILE INTERPRETATION COMPLETE               ")')
    write(logfile,'("================================================================")')

    write(logfile,'(/"* Finished reading graph moves file..."/)')

    return
  end Subroutine ReadGraphMoves


  !
  !******************************************************************
  !> ReadForbiddenGraphs
  !!
  !! Reads the set of graphs which are not allowed to exist.
  !!
  !! - forbidfile: The file containing the forbidden
  !!               graph information.
  !!
  !******************************************************************
  !
  Subroutine ReadForbiddenGraphs( forbidfile )
    implicit none
    character(len=50) :: forbidfile
    logical :: there
    integer :: line, ios, pos, i, pos2
    character(len=100) :: buffer, label, comment, buffer2, label2

    write(logfile,'(/"- Reading forbidden graph file..."/)')
    write(logfile,'("- Forbidden graphs input file:",1x,A/)')trim(forbidfile)

    ! Check and open movefile.
    !
    inquire( file = forbidfile, EXIST=THERE )
    if (.not.there)stop '* Specified forbidfile does not exist'
    Open(20, file = forbidfile, status = 'unknown')


    ! Loop over lines in the file, reading keywords and arguments.
    !
    nforbid = 0
    line = 0
    ios = 0
    do while (ios == 0)

      read(20, '(A)', iostat=ios) buffer
      if (ios == 0) then
        line = line + 1

        ! Find the first space
        pos = scan(buffer, ' ')
        label = buffer(1:pos)
        buffer = buffer(pos+1:)


        ! If pos == 1, there is a space in the first position and we assume a blank line.
        if (pos == 1) then
          write(logfile,*)label,buffer
          CYCLE
        endif


        ! Read input based on label.
        !
        select case (label)


          ! Comment line
        case ('# ')
          read(buffer,*,iostat=ios) comment
          write(logfile,*)trim(label)," ",trim(buffer)


          ! Graph move
        case ('forbid')

          nforbid = nforbid + 1

          write(logfile,'("*** Graph move number:",1x,i4," *** ")')nforbid


          ! The first entry is the number of atoms.
          !
          read(20,'(A)',iostat=ios)buffer2
          pos2 = scan(buffer2, ' ')
          label2 = buffer2(1:pos2)
          buffer2 = buffer2(pos2+1:)
          read(buffer2,*,iostat=ios)naforbid(nforbid)
          write(logfile,'("- Number of atoms:",1x,i4)')naforbid(nforbid)


          ! Next, read the forbidden graph - the first line is a '-', followed
          ! by the namove x namove matrix on the next namove lines.
          !
          write(logfile,'("- Forbidden graph: ")')
          read(20,'(A)',iostat=ios)buffer2
          do i = 1, naforbid(nforbid)
            read(20,*,iostat=ios)gforbid(nforbid,i,1:naforbid(nforbid))
            write(logfile,*)gforbid(nforbid,i,1:naforbid(nforbid))
          enddo
          read(20,'(A)',iostat=ios)buffer2


          ! Now read the labels of the atoms, if any.
          !
          read(20,'(A)',iostat=ios)buffer2
          pos2 = scan(buffer2, ' ')
          label2 = buffer2(1:pos2)
          buffer2 = buffer2(pos2+1:)
          read(buffer2,*,iostat=ios)(forbidlabel(nforbid,i),i=1,naforbid(nforbid))
          write(logfile,'("- Allowed atom labels: ",A)')(forbidlabel(nforbid,i),i=1,naforbid(nforbid))

        end select

      endif
    enddo

    write(logfile,'(/"- Finished reading forbidden graphs file..."/)')


    return
  end Subroutine ReadForbiddenGraphs


  !
  !******************************************************************
  !> SetIODefaults
  !!
  !! Sets default input parameters to "sensible" values.
  !!
  !******************************************************************
  !
  Subroutine SetIODefaults()
    implicit none

    nimage = 10
    NEBiter = 0
    NEBoutfreq = 10
    NEBconv = 1d-4
    CIthresh = 1d-3
    temperature = 10.d0
    nebspring = 0.05d0
    NEBstep = 20.d0
    calctype = 'pathfind'
    pathinit = 'linear'
    pathoptmethod = 'cineb'
    NEBmethod = 'quickmin'
    startfrompath = .FALSE.
    stripinactive = .FALSE.
    optendsbefore = .FALSE.
    optendsduring = .TRUE.
    NEBrestrend = .FALSE.
    NDOFconstr = 0
    Natomconstr = 0
    PEStype = 'null'
    PESopttype = 'null'
    PESfull = .TRUE.
    nreactivetypes = 0
    gdsspring = 0.025d0
    gdsrestspring = 0.05d0
    nbstrength = 0.03d0
    nbrange = 2.3d0
    kradius = 0.02
    gdstemperature = 100.d0
    gdscoefftemp = 100.d0
    gdscoeffmass = 50000.d0
    timestep = 0.5
    ngmove = 0
    nessentialatoms = 0
    igfunc = 0
    gdsoutfreq = 10
    gdsthresh = 0.05d0
    ngdsrelax = 2500
    gdsdtrelax = 0.1d0
    optaftermove = .FALSE.
    skiprepeats = .FALSE.
    nvalcon = 0
    nrxval = 0
    ngen = 0
    npop = 0
    ncross = 0
    nmut = 0
    nheteromax = 30
    nheteromin = 5
    nelprob = 0
    nmccxs = 2500
    mctemperature = 1000.d0
    mcbondprob = 0.05d0
    nmolpenalty = 5.d0
    ringpenalty3 = 5.d0
    ringpenalty4 = 5.d0
    nheterolimit = 60
    essentialmoves(:) = .FALSE.
    essential(:) = .FALSE.
    forbidgraphs = .FALSE.
    shimmybeads = .FALSE.
    VSthresh = 0.0d0
    idpppath =  .FALSE.
    fraginterpol = .FALSE.
    optfragorient = .FALSE.
    atomidx = (/0,0,0/)
    anebb = 0
    gatherreactivemol = .FALSE.
    lmoldata = .FALSE.
    projforcetype = 2
    readcore = .FALSE.
    simpleopt = .false.
    evbvrep = .false.
    evbtype = 1
    evbiter = 100
    evbalpha1 = 0.5
    evbalpha2 = 2.5
    evbstep = 1d-3
    evbmaxdl = 2.0
    moleculemax = 100
    intra_cutoff = -1.0d0

    return
  end Subroutine SetIODefaults

end Module IO
