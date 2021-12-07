!
!***************************************************************************************
!
!> @brief Potential energy surface calculations.
!
!> Defines subroutines for evaluating energies and forces using various
!! external codes.
!
!***************************************************************************************
!

Module pes

  use constants
  use chemstr
  use globaldata
  implicit none

  integer, parameter :: NLINEMAX = 200
  character (len=6) :: vtype
  character (len=25) :: vfile
  character (len=6) :: vopttype
  character (len=25) :: voptfile
  character (len=100) :: PESExec, PESOPTEXEC
  character (len=100), dimension(NLINEMAX) :: PESlines, PESoptlines
  integer :: nline, nlineopt
  integer :: coordsline, optcoordsline
  save

contains


  !
  !************************************************************************
  !> SetupEnergyCalc
  !!
  !! Initializes the PES calculation.
  !!
  !! - PEStype: Identifies the type of PES calculation, 'ORCA'...
  !! - PESfile: Identifies the file containing the template input file
  !!            for PES evaluation. This should have the characters 'XXX'
  !!            in the line where the coordinates will be input.
  !!
  !************************************************************************
  !
  Subroutine SetupEnergyCalc(PEStype, PESfile, PESExecutable)
    character (len=6) :: PEStype
    character (len=25) :: PESfile
    character (len=100) :: PESExecutable
    logical :: there
    integer :: ios
    character (len=100) :: buffer

    vtype = Pestype
    vfile = PESfile
    PESEXEC = PESExecutable

    ! Read the PES file.
    !
    if (trim(pestype) == 'orca' .or. trim(pestype) == 'dftb'&
        .or. trim(pestype) == 'lammps' .or. trim(pestype) == 'psi4' .or. &
        trim(pestype) == 'molpro') then

      inquire(file = vfile, EXIST=THERE)
      if (.not. there) then
        print*,'* Input file does not exist for SetupEnergyCalc: ', vfile
        stop
      endif
      open(18, file = vfile, status = 'unknown')
      ios = 0
      nline = 0
      do while (ios == 0)
        read(18, '(A)', iostat=ios)buffer
        if (ios == 0) then
          if (buffer(1:3) == 'XXX') then
            nline = nline + 1
            coordsline = nline
          else
            nline = nline + 1
            PESlines(nline) = buffer
          endif
        endif
      enddo
      close(18)

    endif

    return
  end Subroutine SetupEnergyCalc


  !
  !************************************************************************
  !> SetupGeomOpt
  !!
  !! Initializes the geometry optimisation calculation.
  !!
  !! - PESopttype: Identifies the PES for geometry optimisation, 'ORCA'...
  !! - PESoptfile: Identifies the file containing the template input file
  !!               for PES evaluation. This should have the characters 'XXX'
  !!               in the line where the coordinates will be input.
  !!
  !************************************************************************
  !
  Subroutine SetupGeomOpt(PESOpttype, PESOptfile, PESOptExecutable)
    implicit none
    character (len=6) :: PESopttype
    character (len=25) :: PESoptfile
    character (len=100) :: PESoptExecutable
    logical :: there
    integer :: ios
    character (len=100) :: buffer

    vopttype = Pesopttype
    voptfile = PESoptfile
    PESOPTEXEC = PESOptExecutable

    ! Read the PES file.
    !
    if (trim(pestype) == 'orca' .or. trim(pestype) == 'dftb'&
        .or. trim(pestype) == 'lammps' .or. trim(pestype) == 'psi4' .or. &
        trim(pestype) == 'molpro') then

      inquire(file=voptfile, EXIST=THERE)
      if (.not. there) then
        print*,'* Input file does not exist for SetupGeomOpt: ', voptfile
        stop
      endif
      open(18, file = voptfile, status = 'unknown')
      ios = 0
      nlineopt = 0
      do while (ios == 0)
        read(18, '(A)', iostat=ios) buffer
        if (ios == 0) then
          if (buffer(1:3) == 'XXX') then
            nlineopt = nlineopt + 1
            optcoordsline = nlineopt
          else
            nlineopt = nlineopt + 1
            PESoptlines(nlineopt) = buffer
          endif
        endif
      enddo
      close(18)

    else if (trim(pesopttype) == 'uff') then
      ! If pesopttype is UFF, there isn't actually anything to do!.
      !
    endif

    return
  end Subroutine SetupGeomOpt


  !
  !************************************************************************
  !> AbInitio
  !!
  !! Calculates Ab-initio properties for the chemical structure in cx.
  !!
  !! - cx: Chemical structure object.
  !! - abtype: 'ener' = energy calculation ; 'grad' = energy+gradient
  !!           'hess' = energy+gradient+hessian ; 'optg' = optimization
  !! - success: A logical flag indicating whether or not the PES
  !!            evaluation was successful. For example, might not
  !!            have converged the SCF cycles in DFT/HF, etc.
  !!
  !************************************************************************
  !
  Subroutine AbInitio(cx, abtypein, success)
    implicit none
    logical :: minimize, success
    type(cxs) :: cx
    type(cxs), dimension(:), allocatable :: cxtemp
    character (len=6) :: ptype
    real(8) :: Vtot
    integer :: ii, natom, nmol, i, j
    real(8), dimension(:), allocatable :: xtemp, ytemp, ztemp
    character(len=2), dimension(:), allocatable :: labeltemp
    character(len=4)              :: abtype
    character(len=4),intent(in)   :: abtypein

    abtype = abtypein

    ! If we're running optg or hess calculations, turn on minimization.
    !
    minimize = .false.
    if (index(abtype, 'optg') .ne. 0) minimize = .true.
    if (index(abtype, 'hess') .ne. 0) minimize = .true.

    ! If there is only 1 atom, turn off minimization!
    !
    if (index(abtype, 'optg') .ne. 0) then
      if (cx%na .eq. 1) abtype = 'ener'
    endif

    ! Decide which calculation type we need to run, based on the calculation
    ! type requested for either single-point or geometry optimization.
    !
    if (.not. minimize) then
      ptype = vtype
    else if (minimize) then
      ptype = vopttype
    endif

    ! CHECK - HESSIAN ONLY IMPLEMENTED FOR ORCA AT MOMENT (SH)
    !
    if ((index(abtype, 'hess') .ne. 0) .and. (ptype .ne. 'orca')) then
      stop 'ERROR: HESSIAN IS ONLY IMPLEMENTED FOR ORCA CALCULATION'
    endif

    ! Are we running PES calculation for the full system? (PESfull = .TRUE.)
    !
    if (PESfull) then

      ! Run the calculation.
      select case (ptype)
        case('orca')
          call ORCAcalc(cx, abtype, success)

        case('dftb')
          call DFTBcalc(cx, minimize, success)

        case('lammps')
          call LAMMPScalc(cx, abtype, success)

        case('psi4')
          call PSI4calc(cx, minimize, success)

        case('molpro')
          call MOLPROcalc(cx, abtype, success)

        case('uff')
          call UFFcalc(cx, minimize, success)

        case('null')
          cx%vcalc = 0.d0
          cx%dvdr(:, :) = 0.d0

        case default
          stop '* Unknown calculation type in PEScalc'

      end select

      ! ...or are we running PES evaluations for each molecule in cx?
      !
    else if (.not. PESfull) then

      ! Make individual CXS objects for each molecule.
      !
      nmol = cx%nmol
      allocate(cxtemp(nmol))
      allocate(labeltemp(NAMAX))
      allocate(xtemp(NAMAX), ytemp(NAMAX), ztemp(NAMAX))
      do i = 1, nmol
        natom = cx%namol(i)
        do j = 1, natom
          ii = cx%molid(i, j)
          labeltemp(j) = cx%AtomLabel(ii)
          xtemp(j) = cx%r(1, ii)
          ytemp(j) = cx%r(2, ii)
          ztemp(j) = cx%r(3, ii)
        enddo
        call CreateCXS(cxtemp(i), natom, labeltemp, xtemp, ytemp, ztemp)
        cxtemp(i)%molcharge(1) = cx%molcharge(i)
        cxtemp(i)%nmol = 1
      enddo

      ! For each CXS object, calculate the energy.
      !
      do i = 1, nmol

        if (cxtemp(i)%na .eq. 1) then
          abtype='ener'
        endif

        ! Run the calculation.
        select case (ptype)

          case('orca')
            call ORCAcalc(cxtemp(i), abtype, success)

          case('dftb')
            call DFTBcalc(cxtemp(i), minimize, success)

          case('lammps')
            call LAMMPScalc(cxtemp(i), abtype, success)

          case('psi4')
            call PSI4calc(cxtemp(i), minimize, success)

          case('molpro')
            call MOLPROcalc(cxtemp(i), abtype, success)

          case('uff')
            call UFFcalc(cxtemp(i), minimize, success)

          case('null')
            cxtemp(i)%vcalc = 0.d0
            cxtemp(i)%dvdr(:, :) = 0.d0

          case default
            Stop '* Unknown calculation type in PEScalc'

        end select
      enddo

      ! Now recombine the results for each molecule.
      !
      Vtot = 0.d0
      do i = 1, nmol
        Vtot = Vtot + cxtemp(i)%vcalc
        natom = cx%namol(i)
        cx%molen(i) = cxtemp(i)%vcalc
        print *, 'MOLEN VALUES: ', i, cx%molen(i)
        do j = 1, cxtemp(i)%na
          ii = cx%molid(i, j)
          cx%r(1:3, ii) = cxtemp(i)%r(1:3, j)
          cx%dvdr(1:3, ii) = cxtemp(i)%dvdr(1:3, j)
        enddo
        call DeleteCXS(cxtemp(i))
      enddo
      cx%vcalc = Vtot

      deallocate(cxtemp)
      deallocate(labeltemp)
      deallocate(xtemp)
      deallocate(ytemp)
      deallocate(ztemp)

    endif

    return
  end Subroutine AbInitio

  !
  !************************************************************************
  !> MOLPROcalc
  !!
  !! Performs a single-point energy and force calculation using MOLPRO.
  !!
  !! - cx: Chemical structure object.
  !! - abtype: indicates what time of calculation we shold perform
  !!
  !************************************************************************
  !
  Subroutine MOLPROcalc(cx,abtype,success)
    implicit none
    type(cxs) :: cx
    real(8):: xx,yy,zz, rdum
    logical :: minimize, flag, there,success
    character (len=100) :: buffer(NLINEMAX), cmsg, exec, cline
    character (len=125) :: str
    integer :: n, iline, i, j, estat, cstat, ios, gat, sz, TotalCharge, TotalSpin, wfi
    integer :: nel, opti
    character (len=20) :: cdum, cnum, gstr1, gstr2, gstr3, gstr4, string, nelec, spin
    character (len=50) :: string1, string3, ia1, ia2
    character (len=4) :: abtype
    success = .true.

    minimize = .false.
    if (abtype .eq. 'optg') minimize = .true.
    ! Create the file template based on the calculation type.
    !
    call CreateFileTemplate(minimize, exec, n, buffer, iline)

    call execute_command_line("rm -f f.out e.out temp.*",exitstat = estat, cmdstat = cstat, &
                              cmdmsg = cmsg )
    ! Open temporary input file called temp.in.
    !
    open(21,file='temp.in',status='unknown')
    do i = 1, iline-1
      write(21,*) trim(buffer(i))
    enddo
    !print*, 'MINIMIZE', minimize
    do j = 1, cx%na
      if (minimize) then
        write(21,'(A)',advance='no') trim(cx%atomlabel(j))//',  ,'
        do i = 1, 3
          if (cx%fixedatom(j) .or. cx%fixeddof((j-1)*3+i)) then
            write(gstr1,*) j ; write(gstr2,*) i
            write(21,'(A)',advance='no') trim(cx%atomlabel(j))//trim(adjustl(gstr1))//trim(adjustl(gstr2))//' '
          else
            write(21,'(f14.8)',advance='no') cx%r(i,j) * bohr_to_ang
          endif
        enddo
        write(21,*)
      else
        xx = cx%r(1,j) * bohr_to_ang
        yy = cx%r(2,j) * bohr_to_ang
        zz = cx%r(3,j) * bohr_to_ang
        write(21,'(a2,",,",1x,3(f14.8,", "))')cx%atomlabel(j),xx,yy,zz
      endif
    enddo
    write(21,*) '}'
    if (minimize) then
      do j = 1, cx%na
        do i = 1, 3
          if (cx%fixedatom(j) .or. cx%fixeddof((j-1)*3+i)) then
            write(gstr1,*) j ; write(gstr2,*) i ; write(gstr3,'(f14.8)') cx%r(i,j) * bohr_to_ang
            write(21,*) trim(cx%atomlabel(j))//trim(adjustl(gstr1))//trim(adjustl(gstr2))//'='//trim(adjustl(gstr3))
          endif
        enddo
        write(21,*)
      enddo
    endif

    nel = 0
    do j = 1, cx%na
      nel = nel + LabelToNumber(cx%atomlabel(j))
    enddo
    do j = 1, cx%nmol
      nel = nel - cx%molcharge(j)
    enddo
    write(nelec,*) nel
    nelec = trim(adjustl(nelec))
    if (mod(nel,2).eq.0) then
      spin = '0'
    else
      spin = '1'
    endif
    wfi = 0
    do i = 1, n
      if ( index(buffer(i),'wf') .ne. 0 ) wfi = i
    enddo
    do i = iline+2, wfi-1
      write(21,*) trim(buffer(i))
    enddo
    write(21,*) 'wf,'//trim(nelec)//',1,'//trim(spin)
    if (minimize .and. cx%na .gt. 1) then
      do i = wfi+1, n
        write(21,*) trim(buffer(i))
        if (index(buffer(i),'optg') .ne. 0 ) then
          opti = i
          exit
        endif
      enddo
      if (any(cx%fixedatom) .or. any(cx%fixeddof)) then
        str = ''
        write(21,'(A)',advance='no')  'inactive'
        do j = 1, cx%na
          do i = 1, 3
            if (cx%fixedatom(j) .or. cx%fixeddof((j-1)*3+i)) then
              write(gstr1,*) j ; write(gstr2,*) i
              write(21,'(A)',advance='no') ','//trim(cx%atomlabel(j))//trim(adjustl(gstr1))//trim(adjustl(gstr2))
            endif
          enddo
        enddo
        write(21,*)
      endif
      do i = opti+1,n
        write(21,*) trim(buffer(i))
      enddo
    else
      do i = wfi+1, n
        write(21,*) trim(buffer(i))
      enddo
    endif
    if (cx%na .gt. 1) then
      if (abtype .ne. 'ener') write(21,*) 'force'
      if ( minimize ) write(21,*) 'put,xyz,molpro_opt.xyz,new'
      if (abtype .ne. 'ener') write(21,*) 'put,xyzgrad,molpro_grad.xyz,new'
    endif
    close(21)
    ! Molpro
    str = trim(Exec) // " temp.in > pes_molpro.out"
    Call ExecuteCalculation(str, success)
    if (.not.success)stop
    open(21,file='molpro_grad.xyz',status='old')
    read(21,*) ; read(21,*)
    read(21,'(A400)',iostat=ios) cline
    do i = 1, cx%na
      read(cline,'(A3,3F17.10,F6.2,3F14.7)') cdum,rdum,rdum,rdum,rdum,cx%dvdr(1:3,i)
      cx%dvdr(1:3,i) = -cx%dvdr(1:3,i)/(ang_to_bohr*au_to_ev)
    enddo
    close(21)
    if (minimize) then
      open(21,file='molpro_opt.xyz',status='old')
      read(21,*) ; read(21,*)

      do i = 1, cx%na
        read(21,'(A4,3F20.10)') cdum,cx%r(1:3,i)
        cx%r(1:3,i) = cx%r(1:3,i)*ang_to_bohr
      enddo
      close(21)
    endif
    open(21,file='temp.out',status='old',position='append')
    backspace(21)
    there = .false.
    do while(.true.)
      backspace(21) ; backspace(21)
      read(21,'(A)') cline
      cline = trim(adjustl(cline))
      j = scan(cline,' ') ; if (j .eq. 0) j = len(trim(cline))
      if (is_numeric(cline(:j)) ) then
        read(cline(:j),'(F16.8)') cx%vcalc
        there = .true.
        exit
      endif
    enddo
    if (.not. there) success = .false.
    return
  end Subroutine MOLPROcalc


  !
  !************************************************************************
  !> ORCAcalc
  !!
  !! Performs a single-point energy and force calculation using ORCA.
  !!
  !! - cx: Chemical structure object.
  !! - abtype: indicates what type of calculation we shold perform
  !! - success: return flag indicating non-convergence of calculation
  !!            if success = .FALSE.
  !!
  !************************************************************************
  !
  Subroutine ORCAcalc(cx,abtype,success)
    implicit none
    type(cxs) :: cx
    real(8):: xx,yy,zz
    logical :: minimize, flag, there, success
    character(len=4)    :: abtype
    character (len=100) :: buffer(NLINEMAX), cmsg, exec, buffer_store(NLINEMAX)
    character (len=125) :: str
    integer :: n, iline, i, j, estat, cstat, ios, gat, sz, TotalCharge,TotalSpin
    integer :: nel
    character (len=20) :: cdum, cnum, gstr1, gstr2, gstr3, gstr4, string
    character (len=50) :: string1, string3, ia1, ia2
    success = .true.

    minimize = .false.
    if (abtype .eq. 'optg' .or. abtype .eq. 'hess') minimize = .true.


    ! Create the file template based on the calculation type.
    !
    Call CreateFileTemplate(minimize, exec, n, buffer, iline)

  !  if (abtype .eq. 'hess') minimize = .false.  ! DODGY HACK !

    ! Make sure that the output files f.out, e.out and temp.* are all removed.
    ! Otherwise, ORCA might try to restart....
    !
    call execute_command_line("rm -f f.out e.out temp.* temp_*",exitstat = estat, &
                              cmdstat = cstat, &
                              cmdmsg = cmsg )


    ! Open temporary input file called temp.in.
    !
    open(21,file='temp.in',status='unknown')


    ! Write the file. Note: our coordinates are in atomic units (Bohr), and
    ! convert them to Angstroms for use in ORCA.
    !
    buffer_store(:) = buffer(:)
    do i = 1, n

      if ( index(buffer(i),'!') .ne. 0 ) then

        ! If single atom, remove kdiis...
        !   if ( cx%na .eq. 1 .and. ( index((buffer(i)),'KDIIS') .ne. 0 )) then
        !    write(21,*) Replace_Text((buffer(i)),'KDIIS',' ') ; cycle
        !   Endif



       ! add gradient keyword if required
    !   if ( abtype .eq. 'grad' .and. ( index((buffer(i)),'ENGRAD') .eq. 0 )) then
    !    write(21,*) trim(buffer(i))//' ENGRAD ' ; cycle
    !   endif

       ! opt does not require engrad kwy
    !   if ( abtype .eq. 'optg' .and. ( index((buffer(i)),'ENGRAD') .ne. 0 )) then
    !    write(21,*) Replace_Text((buffer(i)),'ENGRAD',' ') ; cycle
    !   Endif

       ! opt does not require FREQ kwy
    !   if ( abtype .eq. 'optg' .and. ( index((buffer(i)),'FREQ') .ne. 0 )) then
    !    write(21,*) Replace_Text((buffer(i)),'FREQ',' ') ; cycle
    !   Endif

       ! hess does not require engrad kwy
    !   if ( abtype .eq. 'hess' .and. ( index((buffer(i)),'ENGRAD') .ne. 0 )) then
    !    write(21,*) Replace_Text((buffer(i)),'ENGRAD',' ') ; cycle
    !   Endif

       ! remove engrad if its only an energy calculation..
    !   if ( abtype .eq. 'ener' .and. ( index((buffer(i)),'ENGRAD') .ne. 0 )) then
    !    write(21,*) Replace_Text((buffer(i)),'ENGRAD',' ') ; cycle
    !   endif

       ! add FREQ keyword if required for hessian
       if ( abtype .eq. 'hessian' .and. ( index((buffer(i)),'FREQ') .eq. 0 )) then
        write(21,*) trim(buffer(i))//' FREQ ' ; cycle
       endif

      endif


      ! NEW - account for charges and spin. Here, we replace the '* xyz ...' line read from the
      ! input file with the relevant charge and spin states calculated for the molecule.
      !
      if (adjustl(buffer(i)(1:4)) == '*xyz' .or. adjustl(buffer(i)(1:5)) == '* xyz') then

        ! Count electrons...
        !
        nel = 0
        do j = 1, cx%na
          nel = nel + LabelToNumber(cx%atomlabel(j))
        enddo
        TotalCharge = 0
        do j = 1, cx%nmol
          nel = nel - cx%molcharge(j)
          TotalCharge = TotalCharge + cx%molcharge(j)
        enddo
        if (mod(nel,2) == 0) then
          TotalSpin = 1
        else
          TotalSpin = 2
        endif
        write(ia1,'(i2)')TotalCharge
        write(ia2,'(i2)')TotalSpin
        buffer(i) = '* xyz ' // trim(ia1) // ' ' // trim(ia2)

      endif

      ! Are we at the line for coordinate input?
      !
      if ( i == iline ) then

        do j = 1, cx%na
          xx = cx%r(1,j) * bohr_to_ang
          yy = cx%r(2,j) * bohr_to_ang
          zz = cx%r(3,j) * bohr_to_ang
          write(21,'(a2,1x,3f14.8)')cx%atomlabel(j),xx,yy,zz
        enddo

      !   if (pesfull) then  ! Constraints implemented for full system only....
      !
      !     if (minimize .and. (any(cx%fixedatom) .or. any(cx%fixeddof))) then
      !       do j = 1, cx%na
      !         print*,'J = ',j,cx%fixeddof((j-1)*3+1),cx%fixeddof((j-1)*3+2),cx%fixeddof((j-1)*3+3)
      !         if ( cx%fixedatom(j) ) cycle
      !         if (cx%fixeddof((j-1)*3+1)) then
      !           write(gstr1,'(F14.8)') bohr_to_ang*cx%r(1,j)
      !           write(gstr2,'(F14.8)') bohr_to_ang*cx%r(1,j)+10.0d0
      !           write(gstr3,'(F14.8)') bohr_to_ang*cx%r(2,j)+15.0d0
      !           write(gstr4,'(F14.8)') bohr_to_ang*cx%r(3,j)
      !           write(21,'("DA     '//gstr1//' '//gstr3//' '//gstr4//'")')
      !           write(21,'("DA     '//gstr2//' '//gstr3//' '//gstr4//'")')
      !         endif
      !         if (cx%fixeddof((j-1)*3+2)) then
      !           write(gstr1,'(F14.8)') bohr_to_ang*cx%r(2,j)
      !           write(gstr2,'(F14.8)') bohr_to_ang*cx%r(2,j)+10.0d0
      !           write(gstr3,'(F14.8)') bohr_to_ang*cx%r(1,j)+15.0d0
      !           write(gstr4,'(F14.8)') bohr_to_ang*cx%r(3,j)
      !           write(21,'("DA     '//gstr3//' '//gstr1//' '//gstr4//'")')
      !           write(21,'("DA     '//gstr3//' '//gstr2//' '//gstr4//'")')
      !         endif
      !         if (cx%fixeddof((j-1)*3+3)) then
      !           write(gstr1,'(F14.8)') bohr_to_ang*cx%r(3,j)
      !           write(gstr2,'(F14.8)') bohr_to_ang*cx%r(3,j)+10.0d0
      !           write(gstr3,'(F14.8)') bohr_to_ang*cx%r(1,j)+15.0d0
      !           write(gstr4,'(F14.8)') bohr_to_ang*cx%r(2,j)
      !           write(21,'("DA     '//gstr3//' '//gstr4//' '//gstr1//'")')
      !           write(21,'("DA     '//gstr3//' '//gstr4//' '//gstr2//'")')
      !         endif
      !       enddo
      !     endif
      !
      !   else
      !      if (index(buffer(i),'xyz') .ne. 1 .and. minimize .and.(any(cx%fixedatom) .or. any(cx%fixeddof))) then
      !       gat = 1
      !       write(21,'("%geom Constraints")')
      !       do j = 1, cx%na
      !         if ( cx%fixedatom(j) ) then
      !           write(string,*) j-1
      !           write(21,'("        {C '//trim(adjustl(string))//' C}")')
      !         else
      !           if (cx%fixeddof((j-1)*3+1)) then
      !             write(string,*) j-1
      !             write(gstr1,*) cx%na+gat-1
      !             write(gstr2,*) cx%na+gat
      !             gat = gat + 2
      !             write(21,'("        {A '//trim(string)//trim(gstr1)//' '//trim(gstr2)//'  C} ")')
      !             write(21,'("        {C '//trim(gstr1)//' C} ")')
      !             write(21,'("        {C '//trim(gstr2)//' C} ")')
      !           endif
      !           if (cx%fixeddof((j-1)*3+2)) then
      !             write(string,*) j-1
      !             write(gstr1,*) cx%na+gat-1
      !             write(gstr2,*) cx%na+gat
      !             gat = gat + 2
      !             write(21,'("        {A '//trim(string)//trim(gstr1)//' '//trim(gstr2)//'  C} ")')
      !             write(21,'("        {C '//trim(gstr1)//' C} ")')
      !             write(21,'("        {C '//trim(gstr2)//' C} ")')
      !           endif
      !           if (cx%fixeddof((j-1)*3+3)) then
      !             write(string,*) j-1
      !             write(gstr1,*) cx%na+gat-1
      !             write(gstr2,*) cx%na+gat
      !             gat = gat + 2
      !             write(21,'("        {A '//trim(string)//trim(gstr1)//' '//trim(gstr2)//'  C} ")')
      !             write(21,'("        {C '//trim(gstr1)//' C} ")')
      !             write(21,'("        {C '//trim(gstr2)//' C} ")')
      !           endif
      !         endif
      !       enddo
      !       write(21,'("     end")')
      !       write(21,'("end")')
      !       write(21,*)buffer(i)
      !     endif
      !   endif
       else
        write(21,*)buffer(i)
      endif

    enddo

    ! Replace the buffer lines with the original (???)
    buffer(:) = buffer_store(:)

    ! Close ORCA file.
    !
    close(unit=21)


    ! Run the ORCA calculation using the relevant executable...
    !
    str = trim(Exec) // " temp.in > pes_orca.out"
    Call ExecuteCalculation(str, success)

    !if (.not.success)stop ????????

  !  if (abtype .ne. 'ener') then
      ! Extract the forces from the 'temp.engrad' file.
      !
      string1 = "grep 'The current gradient' temp.engrad -A"
      write(cnum,'(I5)')3*cx%na+2
      string3 = " > f.out"
      str = string1 // adjustl(trim(cnum))//string3
      Call ExtractData(str)
      inquire(file='f.out',SIZE=sz)

      if (sz <= 0) then
        success = .FALSE.
        return
      endif

      ! Read the forces - the output is already in atomic units.
      !
      Call ReadForcesIndividually(cx, 2)

      ! Extract the energy from the 'temp.engrad' file.
      !
      str = "grep 'The current total energy' temp.engrad -A 2 > e.out"
      Call ExtractData(str)

      ! Read in the energy from file (in atomic units).
      !
      Call ReadEnergy(cx, 2)
  !  else
      !str = "grep 'Value' temp_property.txt  | tail -1 | awk '{print $2}' > e.out"
      ! DIFFERENT VERSION OF ORCA:
  !    str = "grep 'SCF Energy:' temp_property.txt  | tail -1 | awk '{print $3}'> e.out "
  !    Call ExtractData(str)
      ! Read in the energy from file (in atomic units).
      !
  !    Call ReadEnergy(cx, 0)
  !  endif


    ! Read the optimized coordinates from temp.xyz
    !
    if (minimize) then

      ! TEMP
      if (abtype .ne. 'hess')then
       Call ReadOptimizedCoordinates(cx, 'temp.xyz')
      endif

       ! Temporary fix - if requiring Hessian, use FREQ in orca.min.
       !
       if (abtype .eq. 'hess') then
           Call ReadHessian(cx, 'temp.hess')
       endif

    endif

    return
  end Subroutine ORCAcalc



  !
  !************************************************************************
  !> DFTBcalc
  !!
  !! Performs a single-point energy and force calculation using DFTB.
  !! Note that the appropriate Slater-Koster files should be present in the
  !! run directory.
  !!
  !! - cx: Chemical structure object.
  !! - minimize: logical flag indicating whether or not to run a
  !!             geometry optimization calculation.
  !! - success: Flag indicating success of calculation.
  !!
  !************************************************************************
  !
  Subroutine DFTBcalc(cx, minimize, success)
    implicit none
    type(cxs) :: cx
    real(8):: xx,yy,zz
    logical :: minimize, success, there
    character (len=100) :: cmsg, exec, str
    character (len=100), dimension(NLINEMAX) :: buffer
    integer :: n, iline, i, j, estat, cstat, ios, ntype
    integer :: ifound, idum, sz
    character (len=2) :: atype(NAMAX), c2
    character (len=20) :: cdum
    character (len=50) :: string1, string3
    character (len=3) :: cnum

    ! Create the file template based on the calculation type.
    !
    call CreateFileTemplate(minimize, exec, n, buffer, iline)

    ! Set flag indicating whether or not the calculation exectued correctly.
    !
    success = .true.

    ! Clear out the current dftb_in.hsd
    !
    call execute_command_line("rm -f dftb_in.hsd detailed.out f.out", exitstat = estat, cmdstat = cstat, &
        cmdmsg = cmsg, WAIT = .TRUE. )

    ! Write the dftb_in.hsd file.
    !
    open(21, file='dftb_in.hsd', status='unknown')
    write(cnum, '(I3)') cx%na
    write(21, *) buffer(1)

    ! Set number of atoms.
    !
    write(21, *) adjustl(trim(cnum)),' C' !< This sets number of atoms and defines a "cluster" calculation (no PBC)

    ! Identify the different atom types.
    !
    atype(1) = cx%atomlabel(1)
    ntype = 1
    do i = 2, cx%na
      ifound = 0
      do j = 1, ntype
        if (cx%atomlabel(i) == atype(j)) ifound = 1
      enddo
      if (ifound == 0) then
        ntype = ntype + 1
        atype(ntype) = cx%atomlabel(i)
      endif
    enddo

    ! Catenate atom labels together and write to file
    !
    str = atype(1)
    do i = 2, ntype
      str = trim(str) // ' ' // trim(atype(i))
    enddo
    write(21, *) str

    ! Now print coordinates...
    !
    do i = 1, cx%na
      do j = 1, ntype
        if (cx%atomlabel(i) == atype(j)) then
          xx = cx%r(1, i) * bohr_to_ang
          yy = cx%r(2, i) * bohr_to_ang
          zz = cx%r(3, i) * bohr_to_ang
          write(21, '(i3,1x,i2,1x,3f14.8)') i, j, xx, yy, zz
        endif
      enddo
    enddo

    do i = iline+1, n
      write(21, *) buffer(i)
    enddo

    ! Close DFTB file.
    !
    close(unit=21)

    ! Run the calculation - note that DFTB+ just assumes that the input file is called
    ! dftb_in.hsd, so doesn't need any arguments.
    !
    str = trim(Exec) // " > dftb.output "
    call ExecuteCalculation(str, success)

    ! Check that f.out and e.out exists....
    !

    ! Extract the forces from the 'detailed.out' file.
    !
    str = "grep 'Total Forces' detailed.out -A "//adjustl(trim(cnum))//" > f.out"
    call ExtractData(str)
    inquire(file='f.out', SIZE=sz)

    if (sz <= 0) then
      success = .FALSE.
      return
    endif

    ! Read the forces - the output is already in atomic units.
    !
    call ReadForcesTogether(cx, 1)

    ! Change forces to derivatives.
    !
    cx%dvdr(:, :) = -cx%dvdr(:, :)

    ! Extract the energy from the 'detailed.out' file.
    !
    str = "grep 'Total energy' detailed.out | awk '{print $3}'> e.out"
    call ExtractData(str)
    inquire(file='e.out', SIZE=sz)
    if (sz <= 0) then
      success = .FALSE.
      return
    endif

    ! Read in the energy from file (in atomic units).
    !
    call ReadEnergy(cx, 0)

    ! If we're running optimization, we need to read in the optimized coordinates.
    !
    if (minimize) then
      call ReadOptimizedCoordinates(cx, 'geo_end.xyz')
    endif

    return
  end Subroutine DFTBcalc


  !
  !************************************************************************
  !> UFFcalc
  !!
  !! Performs a single-point energy and force calculation using UFF, as
  !! implemented in openbabel (obabel).
  !!
  !! - cx: Chemical structure object.
  !! - minimize: logical flag indicating whether or not to run a
  !!             geometry optimization calculation.
  !!
  !!             For openbabel, we must have minimize = .TRUE.
  !!
  !! - success: Flag indicating success of calculation.
  !!
  !! DEPENDENCIES: Requires babel and obminimize.
  !!
  !************************************************************************
  !
  Subroutine UFFcalc(cx,minimize,success)
    implicit none
    type(cxs) :: cx
    real(8):: xx, yy, zz, del, rr(3,cx%na), vm, vp
    logical :: minimize, success, there
    character (len=100) :: buffer(NLINEMAX), cmsg, exec
    character (len=100) :: str
    integer :: n, iline, i, j, k, estat, cstat, ios, ntype, l
    integer :: ifound, idum
    character (len=2) :: atype(NAMAX), c2
    character (len=50) :: cdum
    character (len=50) :: string1, string3
    character (len=3) :: cnum

    if (.not. minimize) then
         Call PrintCXStofile(cx,'temporary.xyz',0.d0)
         call system("babel -ixyz temporary.xyz -opdb tempi.pdb &>/dev/null")
         call system("obenergy -ff uff tempi.pdb | grep 'TOTAL ENE' > f.out")
         open(24,file='f.out',status='old') ; read(24,'(A50)') string1 ; close(24)
         string1 = string1(index(string1,'=')+1:index(string1,'k')-1)
         read(string1,*) cx%vcalc
         cx%vcalc = cx%vcalc / au_to_kjmol
       del = 0.025
       cx%dvdr = 0.0
       ! calculates the gradient with uff - slow
       !
       open(24,file='temporary.xyz',status='unknown')
       do i = 1, cx%na
         do k = 1, 3
           cx%r(k,i) = cx%r(k,i)-del
           write(24,'(I5)') cx%na
           write(24,'(A5)') ' 0.0000000000'
           do j =1, cx%na
             write(24,'(a2,2x,3(f14.8,2x))') cx%atomlabel(j), (cx%r(l,j)*bohr_to_ang, l = 1,3)
           enddo
           cx%r(k,i) = cx%r(k,i)+2.0d0*del
           write(24,'(I5)') cx%na
           write(24,'(A5)') ' 0.0000000000'
           do j =1, cx%na
             write(24,'(a2,2x,3(f14.8,2x))') cx%atomlabel(j), (cx%r(l,j)*bohr_to_ang, l = 1,3)
           enddo
           cx%r(k,i) = cx%r(k,i)-del
         enddo
       enddo
       close(24)
       call system("babel -ixyz temporary.xyz -opdb tempi.pdb &>/dev/null")
       call system("obenergy -ff uff tempi.pdb 2>/dev/null | grep 'TOTAL ENE' | awk '{print $4}'> f.out")
       open(24,file='f.out',status='old')
       do i = 1, cx%na
         do k = 1, 3
           read(24,*) vm
           vm = vm / au_to_kjmol
           read(24,*) vp
           vp = vp / au_to_kjmol
           cx%dvdr(k,i) = (vp-vm)/(2.0d0*del)
         enddo
       enddo
       close(24)
       success = .TRUE.
    else
       Call PrintCXStofile(cx,'temporary.xyz',0.d0)
       call system("babel -ixyz temporary.xyz -opdb temp.pdb &>/dev/null")
       call system("obminimize -n 8000 -ff UFF temp.pdb > temp2.pdb 2>/dev/null")
!       call system("obminimize -sd -n 8000 -h -ff UFF temp.pdb > temp2.pdb 2>/dev/null")
       call system("babel -ipdb temp2.pdb -oxyz temp2.xyz &>/dev/null")

!       stop

       ! Read the cxs from file.
       !
       Call ReadOptimizedCoordinates(cx, 'temp2.xyz')
       Call GetGraph(cx)
       Call GetMols(cx)
       success = .TRUE.
    endif
    return
  end Subroutine UFFcalc



  !
  !************************************************************************
  !> LAMMPScalc
  !!
  !! Performs a single-point energy and force calculation using LAMMPS.
  !!
  !! - cx: Chemical structure object.
  !! - abtype: indicates what time of calculation we shold perform
  !!
  !************************************************************************
  !
  Subroutine LAMMPScalc(cx,abtype,success)
    implicit none
    type(cxs) :: cx
    real(8):: xx,yy,zz
    logical :: minimize, there, success
    character (len=100) :: buffer(NLINEMAX), cmsg, exec, string
    character (len=100) :: str
    integer :: n, iline, i, j, estat, cstat, ios, loc, q, ntype, ifound
    character (len=20) :: cdum,cnum,pairstyle
    character (len=50) :: string1,string3,gstr1,gstr2,gstr3,gstr4
    character (len=2)  :: atype(NAMAX)
    character (len=4) :: abtype


    success = .true.
    minimize = .false.
    IF (abtype .eq. 'optg') minimize = .true.

    ! Create the file template based on the calculation type.
    !
    Call CreateFileTemplate(minimize, exec, n, buffer, iline)


    ! Make sure that the output files f.out, e.out and temp.* are all removed.
    !
!    call execute_command_line("rm -f f.out e.out temp.*",exitstat=estat, cmdstat = cstat, &
!         cmdmsg = cmsg )
    call execute_command_line("rm -f f.out e.out temp.* *.lmp pes_lmp.log log.*")

    ! Open temporary input file called temp.lmp.
    !

    open(21,file='temp.lmp',status='unknown')

    ! Identify the different atom types.
    !
    atype(1) = cx%atomlabel(1)
    ntype = 1
    do i = 2, cx%na
      ifound = 0
      do j = 1, ntype
        if ( cx%atomlabel(i) == atype(j) ) ifound=1
      enddo
      if (ifound == 0) then
        ntype = ntype + 1
        atype(ntype) = cx%atomlabel(i)
      endif
    enddo


    ! Fills str with atom types
    !
    str = ""
    do q = 1, ntype
       str = trim(str) // " " // atype(q)
    enddo

    ! Write the file. Note: we assume that coordinates are in atomic units (Bohr), and
    ! convert them to Angstroms for use in LAMMPS.
    !
    do i = 1, n

       ! Are we at the line for coordinate input?
       !
       if ( i == iline ) then

          ! Finds maximum values to be accommodated by the box
          !
          xx = maxval(abs(cx%r(:,:))) * bohr_to_ang + 800.0d0

          ! Creates box using the maximum values
          !
          write(21, '(1x,a6,1x,a6,1x,a5,1x,6(2X,f14.8))')&
              'region','thebox','block',-xx,xx,-xx,xx,-xx,xx
          write(21, '(1x,a10,1x,i1,1x,a6)')'create_box',ntype,'thebox'

          ! Writes the atom types
          !
          do q = 1, ntype
             write(21,'(1x,a4,1x,i2,1x,f14.8)')&
                  'mass', q, MASS(LabelToNumber(atype(q)))*me_to_amu
          enddo

          ! Writes the creation of atoms
          !
          do j = 1, cx%na
             xx = cx%r(1,j) * bohr_to_ang
             yy = cx%r(2,j) * bohr_to_ang
             zz = cx%r(3,j) * bohr_to_ang

             ! Assigns current atom type
             !
             do q = 1, ntype
                if ( atype(q) == cx%atomlabel(j) ) then
                   loc = q
                endif
             enddo
             write(21,'(1x,a12,1x,i2,1x,a6,1x,3(f14.8,1x))')&
                  'create_atoms',loc,'single',xx,yy,zz
          enddo
       ! Assigns pairstyle variable
       !
       else if ( index(buffer(i)(1:10),'pair_style') .ne. 0 ) then
          pairstyle = trim(buffer(i)(11:))
          write(21,*)buffer(i)

       ! Adds atom types to pair_coeff line if using reax/c
       !
       else if ( (index(buffer(i)(1:10),'pair_coeff') .ne. 0)&
           .and. (index(pairstyle,'reax/c') /= 0) ) then

          write(21,*)trim(buffer(i)),str
          if (minimize) then
            gstr1 = '' ; gstr2 = '' ; gstr3 = '' ; gstr4 = ''
            do j = 1, cx%na
              if ( cx%fixedatom(j) ) then
                if ( len(trim(gstr1)) .eq. 0 ) gstr1 =  ' group 11 id '
                write(string,*) j
                gstr1 = trim(gstr1)//' '//trim(adjustl(string))
              else
                if (cx%fixeddof((j-1)*3+1)) then
                  if (len(trim(gstr2)) .eq. 0 ) gstr2 = ' group 12 id '
                  write(string,*) j
                  gstr2 = trim(gstr2)//' '//trim(adjustl(string))
                endif
                if (cx%fixeddof((j-1)*3+2)) then
                  if (len(trim(gstr3)) .eq. 0 ) gstr3 = ' group 13 id '
                   write(string,*) j
                  gstr3 = trim(gstr3)//' '//trim(adjustl(string))
                endif
                if (cx%fixeddof((j-1)*3+3)) then
                  if (len(trim(gstr4)) .eq. 0 ) gstr4 = ' group 14 id '
                  write(string,*) j
                  gstr4 = trim(gstr4)//' '//trim(adjustl(string))
                endif
              endif
            enddo
            if (len(trim(gstr1)) .ne. 0 ) then
               write(21,'("'//trim(gstr1)//'")')
               write(21,'(" fix 11 11 setforce 0.0 0.0 0.0 ")')
            endif
            if (len(trim(gstr2)) .ne. 0 ) then
               write(21,'("'//trim(gstr2)//'")')
               write(21,'(" fix 12 12 setforce 0.0 NULL NULL ")')
            endif
            if (len(trim(gstr3)) .ne. 0 ) then
               write(21,'("'//trim(gstr3)//'")')
               write(21,'(" fix 13 13 setforce NULL 0.0 NULL ")')
            endif
            if (len(trim(gstr4)) .ne. 0 ) then
               write(21,'("'//trim(gstr4)//'")')
               write(21,'(" fix 14 14 setforce NULL 0.0 0.0 ")')
            endif
          endif
       else
          write(21,*)buffer(i)
       endif

    enddo
    ! Writes final dump command after minimization
    !
    if ( minimize ) then
      write(21,'(" dump 1 all xyz 1 temp.xyz ")')
      write(21,*) 'dump_modify 1 format line "%s %14.8f %14.8f %14.8f"'
      write(21,'(" dump_modify 1 element'//trim(str)//' ")')
      Write(21,'(" run 0 ")')
    endif

    close(unit=21)
    ! Run the LAMMPS calculation using the relevant executable.
    !
    str = trim(Exec) // " < temp.lmp 2>/dev/null > pes_lmp.log"
    Call ExecuteCalculation(str, success)
    !if (minimize) print*, 'FINISHED EXEC RUN'
    !tempi = tempi + 1
    !write(string1,*) tempi
    !call system('cp temp.lmp  see_daa.lmp_'//adjustl(trim(string1)))

    !write(6,'("HERE1")')

    ! Extract the forces from the 'temp.force' file.
    !
    open(31,file='temp.force',status='unknown')
    ios = 0
    do while (ios .eq. 0)
       read(31,'(A100)',iostat=ios) string
       if (ios/=0) stop '* Error 2 reading f.out in pes.f90/ReadForcesTogether'
       if ( index(string,'fx fy fz').ne.0)exit
    enddo
    !if (minimize) print*, 'EXITED LOOP 1', ios
    do i = 1, cx%na
       read(31,*,iostat=ios)cx%dvdr(1,i),cx%dvdr(2,i),cx%dvdr(3,i)
       if (ios/=0) stop '* Error 2 reading f.out in pes.f90/ReadForcesTogether'
    enddo
    !if (minimize) print*, 'EXITED LOOP 2', ios
    close(31)
    ! Convert force to atomic units.
    !
    ! I think this is right -
    if (abtype .ne. 'ener') &
    cx%dvdr(:,:) = -kcalmol_to_au * cx%dvdr(:,:) / ang_to_bohr
    !write(6,'("HERE2")')
    ! Extract the energy from the 'pes_lmp.log' file.
    !
    open(31,file='pes_lmp.log',status='old')
    ios = 0
    read(31,'(A100)',iostat=ios) string
    cx%vcalc = 0.0d0
    do while (ios .eq. 0)
       if (ios/=0) stop '* Error 2 reading f.out in pes.f90/ReadForcesTogether'
       if (cx%vcalc == 11111.0d0) read(string,*) cx%vcalc
       if ( index(string,'PotEng').ne.0) cx%vcalc = 11111.0d0
       read(31,'(A100)',iostat=ios) string
    enddo
    !if (minimize) print*, 'EXITED LOOP 3', ios
    close(31)
    !write(6,'("HERE2")')

    ! Convert energy to atomic units.
    !
    cx%vcalc = cx%vcalc * kcalmol_to_au

    ! Read the optimized coordinates from temp.xyz
    !
    if (minimize) then
       Call ReadOptimizedCoordinates(cx, 'temp.xyz')
    endif

    return
  end Subroutine LAMMPScalc

  !
  !************************************************************************
  !> PSI4calc
  !!
  !! Performs a single-point energy and force calculation using psi4.
  !!
  !! - cx: Chemical structure object.
  !! - minimize: logical flag indicating whether or not to run a
  !!             geometry optimization calculation.
  !!
  !************************************************************************
  !
  Subroutine PSI4calc(cx,minimize,success)
    implicit none
    type(cxs) :: cx
    real(8):: xx,yy,zz
    logical :: minimize, flag, there,success
    character (len=100) :: buffer(NLINEMAX), cmsg, exec
    character (len=100) :: str
    integer :: n, iline, i, j, estat, cstat, ios
    character (len=20) :: cdum,cnum
    character (len=50) :: string1,string3

    success = .true.

    ! Create the file template based on the calculation type.
    !
    Call CreateFileTemplate(minimize, exec, n, buffer, iline)


    ! Make sure that the output files f.out, e.out and temp.* are all removed.
    !
    call execute_command_line("rm -f f.out e.out temp.*",exitstat = estat, cmdstat = cstat, &
                              cmdmsg = cmsg )


    ! Open temporary input file called temp.in.
    !
    open(21,file='temp.in',status='unknown')

    ! Write the file. Note: we assume that coordinates are in atomic units (Bohr), and
    ! convert them to Angstroms for use in psi4.
    !
    do i = 1, n

      ! Are we at the line for coordinate input?
      !
      if ( i == iline ) then

        do j = 1, cx%na
          xx = cx%r(1,j) * bohr_to_ang
          yy = cx%r(2,j) * bohr_to_ang
          zz = cx%r(3,j) * bohr_to_ang
          write(21,'(a2,1x,3f14.8)')cx%atomlabel(j),xx,yy,zz
        enddo

      else
        write(21,'(A)')trim(buffer(i))
      endif
    enddo

    ! Close psi4 file.
    !
    close(unit=21)


    ! Run the psi4 calculation using the relevant executable.
    !
    str = trim(Exec) // " temp.in"
    Call ExecuteCalculation(str, success)


    ! Read the forces - the output is already in atomic units.
    !
    Call ReadForcesIndividually(cx, 0)


    ! Read in the energy from file (in atomic units).
    !
    Call ReadEnergy(cx, 0)


    ! Read the optimized coordinates from temp.xyz
    !
    if (minimize) then
       Call ReadOptimizedCoordinates(cx, 'temp.xyz')
    endif


    return
  end Subroutine PSI4calc


  Subroutine ReadForcesIndividually(cx, dummyLines)
    implicit none
    type(cxs) :: cx
    logical :: fileExists
    integer :: i, j, dummyLines, ios
    character (len=20) :: cdum

    inquire(file = 'f.out', exist=fileExists)
    if (.not. fileExists) then
      stop '* f.out does not exist in pes.f90/ReadForcesIndividually'
    endif
    open (21, file = 'f.out', status = 'unknown')

    do i = 1, dummyLines
      read(21, *, iostat=ios) cdum
      if (ios/=0) then
        stop '* Error 1 reading f.out in pes.f90/ReadForcesIndividually'
      endif
    enddo

    do i = 1, cx%na
      do j = 1, 3
        read(21, *, iostat=ios) cx%dvdr(j, i)
        if (ios/=0) then
          stop '* Error 2 reading f.out in pes.f90/ReadForcesIndividually'
        endif
      enddo
    enddo
    close(unit = 21)

    return
  end Subroutine


  Subroutine ReadForcesTogether(cx, dummyLines)
     implicit none
     type(cxs) :: cx
     logical :: fileExists
     integer :: i, j, dummyLines, ios
     character (len=20) :: cdum

    inquire(file = 'f.out', exist=fileExists)
    if (.not. fileExists) then
      stop '* f.out does not exist in pes.f90/ReadForcesTogether'
    endif
    open(21, file = 'f.out', status = 'unknown')

    do i = 1, dummyLines
      read(21, *, iostat=ios) cdum
      if (ios/=0) then
        stop '* Error 1 reading f.out in pes.f90/ReadForcesTogether'
      endif
    enddo

    do i = 1, cx%na
      read(21, *, iostat=ios) cx%dvdr(1, i), cx%dvdr(2, i), cx%dvdr(3, i)
      if (ios/=0) then
        stop '* Error 2 reading f.out in pes.f90/ReadForcesTogether'
      endif
!      if (cx%dvdr(1,i) .ne. cx%dvdr(1,i) .or. cx%dvdr(2,i) .ne. cx%dvdr(2,i) .or. &
!          cx%dvdr(3,i) .ne. cx%dvdr(3,i) ) stop 'Error 3 reading f.out in pes.f90/ReadForcesTogether'
    enddo
    close(unit = 21)

    return
  end Subroutine


  Subroutine ReadEnergy(cx, dummyLines)
    implicit none
    type(cxs) :: cx
    logical :: fileExists
    integer :: i, j, dummyLines, ios
    character (len=20) :: cdum

    inquire(file = 'e.out', exist=fileExists)
    if (.not. fileExists) stop '* e.out does not exist in pes.f90/ReadEnergy'
    open (21, file = 'e.out', status = 'unknown')

    do i = 1, dummyLines
      read(21, *, iostat=ios) cdum
      if (ios/=0) stop '* Error 1 reading e.out in pes.f90/ReadEnergy'
    enddo

    read(21, *, iostat=ios) cx%vcalc
    if (ios/=0) stop '* Error 2 reading e.out in pes.f90/ReadEnergy'
    close(unit = 21)

    return
  end Subroutine


  Subroutine ReadHessian(cx, fileName)
    implicit none
    type(cxs) :: cx
    logical :: fileExists
    integer :: i, k, ioff, ii, irow, j
    real(8), dimension(6) :: xx
    character (len=*) :: fileName
    character (len=20) :: cdum
    character (len=3) :: cnum
    character (len=80) :: string1, string3

  !  print*,'NOW HERE';call flush(6)

    inquire(file = fileName, exist=fileExists)
    if (.not. fileExists) then
      return
    endif

   !  print*,'NOW HERE2';call flush(6)

    na = cx%na
    cx%Hessian(:, :) = 0.d0
    string1 = adjustl(trim("grep '$hessian' temp.hess -A "))
    write(6, *) 'NUMBER OF LINES: ', 1 + (3*na+1)*(1 + int(3*dble(na)/5))
    write(cnum, '(I3)') 1 + (3*na+1)*(1 + int(3*dble(na)/5))
    string3 = trim(" > hess.out")
    print *, 'COMMAND: ', string1//adjustl(trim(cnum))//string3
    call system(string1//adjustl(trim(cnum))//string3)
    open(21, file='hess.out', status='unknown')
    read(21, *) cdum
    read(21, *) cdum
    read(21, *) cdum
    ioff = 0
    do k = 1, int(3*dble(na)/5) !read each block of 5.
      do i = 1, 3*na
        irow = i
        read(21, *) idum, xx(1:5)  ! read integer plus 5 entries from file
        do j = 1, 5
          cx%Hessian(irow, j+ioff) = xx(j)
        enddo
      enddo
      ioff = ioff + 5
      read(21, *, END=78) cdum
    enddo
78  continue    


    ! Read remainder row of temp.hess.
    !
    ii = mod(3*na,5)  ! Remaining columns left....
    if (ii > 0) then
      do i = 1, 3*na
        irow = i
        read(21, *) idum, xx(1:ii) ! only reads ii columns...
        do j = 1, ii
          cx%Hessian(irow, j+ioff) = xx(j)
        enddo
      enddo
      close(21)
    endif

    return
  end Subroutine ReadHessian


  Subroutine ReadOptimizedCoordinates(cx, fileName)
    implicit none
    type(cxs) :: cx
    logical :: fileExists
    integer :: i
    real(8) :: xx, yy, zz
    character (len=*) :: fileName
    character (len=20) :: cdum

    inquire(file = fileName, exist=fileExists)
    if (.not. fileExists) then
      stop '* Optimized coordinates file does not exist in pes.f90/ReadOptimizedCoordinates'
    endif

    open(21, file=fileName, status='unknown')
    read(21, *) cdum
    read(21, *) cdum
    do i = 1, cx%na
      read(21, *) cdum, xx, yy, zz
      cx%r(1, i) = xx * ang_to_bohr
      cx%r(2, i) = yy * ang_to_bohr
      cx%r(3, i) = zz * ang_to_bohr
!        print*,'READING OPT COORDINS: ',i,cx%r(1,i),cx%r(2,i),cx%r(3,i)
    enddo
    close(21)

    return
  end Subroutine


  Subroutine CreateFileTemplate(minimize, exec, n, buffer, iline)
    implicit none
    logical :: minimize
    character (len=100) :: exec
    character (len=100), dimension(NLINEMAX) :: buffer
    integer :: iline, n

    if (minimize) then
      exec = pesoptexec
      n = nlineopt
      buffer(:) = pesoptlines(:)
      iline = optcoordsline
    else
      exec = pesexec
      n = nline
      buffer(:) = peslines(:)
      iline = coordsline
    endif

    return
  end Subroutine


  Subroutine ExecuteCalculation(str, success)
    implicit none
    logical :: success
    integer :: cstat, estat
    character (len=100) :: cmsg, str

    estat = 1 ; cstat = 0 ; cmsg = ''

    call execute_command_line(str, EXITSTAT=estat, CMDSTAT=cstat, CMDMSG=cmsg)

    ! print*

    if (estat > 0) then
      success = .false.
      return
    endif

    ! Check for errors!
    !
    if (cstat > 0) then
      print*, "PES calculation failed in pes.f90/ExecuteCalculation with error: ", trim(cmsg)
      stop
    else if (cstat < 0) then
      print*, "PES command unknown in pes.f90/ExecuteCalculation"
      stop
    endif

  end Subroutine


  Subroutine ExtractData(str)
    implicit none
    integer :: cstat, estat
    character (len=100) :: cmsg, str

    call execute_command_line(str, exitstat = estat, cmdstat = cstat, &
        cmdmsg = cmsg)
    if (cstat /= 0) then
      print *, "Error in pes.f90/ExtractData - force/energy grep not working!"
      stop
    endif

  end Subroutine

  !
  !************************************************************************
  !> FireMinimise
  !!
  !! minimises the energy using the FIRE algorithm, incorporating the
  !! graph constraining ppotentials if nebrestend = .true.
  !! the last bit of the minimization, we switch to quickmin, optional
  !!
  !! - cx: Chemical structure object to be optimized
  !! - iter: returns the number of iterations it took to optimize
  !! - all other: parameteres for the graph restraining potential.
  !!
  !************************************************************************
  !
  Subroutine FireMinimise(cx, iterin, gdsrestspring_in, nbstrength_in, nbrange, kradius, nebrestrend)
    implicit none
    integer :: i, j, k, l, m, n, Niter, nd, Nmin, Nmin_count, iter, iterin !Fire variables
    type(cxs) :: cx, cx0
    real(8), dimension(cx%na*3) :: v, F
    real(8) :: f_inc, f_dec, f_alp, alph_start, pfire, alph, delt, Nstep, forg(3), eorg(3), fthresh, dot !Fire variables
    real(8) :: gdsrestspring, nbstrength, nbrange, kradius, Fnorm, nrea, gdsrestspring_in, nbstrength_in, aa, mx, gg, xe, maxf
    logical :: success, nebrestrend

    gdsrestspring = gdsrestspring_in/1.0d0
    nbstrength = nbstrength_in/1.0d0
    ! beware of this - valgrind complains...
    !if (nebrestrend) Call CopyToNewCXS(cx, cx0)

    fthresh = 0.100
    forg = 100.0
    eorg(1) = 100.0; eorg(2) = -520.0 ; eorg(3) = 65.0d0
    nd = cx%na*3
    Nstep = 30.d0
    Niter = 50
    !! by about half way the optimization the contraints should be dissapearing
    !nrea = -log(0.001)/(dble(Niter)/2.0d0)
    maxf = 0.750d0
    mx = 1.0d0-fthresh/4.0d0
    xe = maxf
    aa = atanh(2.0*mx-1.0)/(xe-xe/2.0)

    cx%force = 0.0d0 ; F = 0.0d0 ; v = 0.0d0
    call AbInitio(cx, 'grad', success)
    if (nebrestrend) call GraphConstraints(cx, gdsrestspring, nbstrength, 0.3*nbrange, kradius)
    !if (nebrestrend) Call MatchCXConstraint( cx, cx0, gdsrestspring )
    do j = 1, cx%na
      if (.not. cx%fixedatom(j)) then
        do k = 1, 3
          if (.not. cx%fixeddof((j-1)*3+k)) F((j-1)*3+k) = -cx%dvdr(k, j)
        enddo
      endif
    enddo
    Fnorm = norm2(F)
    Nmin = 5
    f_inc = 1.10d0 ; f_dec = 0.50d0 ; f_alp = 0.990d0 ; alph_start = 0.10d0
    delt = Nstep ; alph = alph_start ; Nmin_count = 0
    cx%p = 0.0d0
    do j = 1, cx%na
      cx%p(1:3, j) = cx%p(1:3, j) + delt * F((j-1)*3+1:(j-1)*3+3)/2.0d0
      v((j-1)*3+1:(j-1)*3+3) = cx%p(1:3, j)/cx%mass(j)
    enddo
    do iter = 1, Niter
      if (sum(forg)/3.0d0 .gt. fthresh*10.0d0) then
        pfire = dot_product(v(:), F(:))
        if (Fnorm .gt. 0.00001) then
          v(:) = (1.0d0-alph)*v(:) + alph*F(:)*norm2(v(:))/Fnorm
        endif
        if (pfire > 0.0d0 .and. Nmin_count > Nmin) then
          delt = min(delt*f_inc, 5*Nstep)
          alph = alph * f_alp
        elseif (pfire <= 0.0d0) then
          Nmin_count = 0
          delt = delt*f_dec
          v(:) = 0.0d0
          cx%p = 0.0d0
          alph = alph_start
        endif
        if (pfire > 0.0d0) Nmin_count = Nmin_count + 1
        do j = 1, cx%na
          cx%r(1:3, j) = cx%r(1:3, j) + delt * v((j-1)*3+1:(j-1)*3+3)
          cx%p(1:3, j) = v((j-1)*3+1:(j-1)*3+3)*cx%mass(j) + delt * F((j-1)*3+1:(j-1)*3+3)
          v((j-1)*3+1:(j-1)*3+3) = cx%p(1:3, j)/cx%mass(j)
        enddo
      else
        ! quickmin once the gradient is quite flat...
        !
        dot = 1.0d0
        if (Fnorm .gt. 0.00000001) then
          dot = dot_product(reshape(cx%p, (/nd/)), F(:)/Fnorm)
        endif
        if (dot .lt. 0.0d0) cx%p = 0.0d0
        if (dot .gt. 0.0d0) cx%p = dot * reshape(F, (/3, cx%na/))/Fnorm
        do j = 1, cx%na
          cx%r(1:3, j) = cx%r(1:3, j) +  Nstep/5.0d0 * ( cx%p(1:3, j) / cx%mass(j))
          cx%p(1:3, j) = v((j-1)*3+1:(j-1)*3+3)*cx%mass(j) + delt * F((j-1)*3+1:(j-1)*3+3)
        enddo
      endif

      call AbInitio(cx, 'grad', success)
      gg = 0.50d0*(1.0d0+tanh((sum(forg)/3.0d0-xe/2.0d0)*aa))
      !if (nebrestrend) Call MatchCXConstraint( cx, cx0, gdsrestspring*gg )
      !if (nebrestrend) Call MatchCXConstraint( cx, cx0, gdsrestspring*exp(-nrea*dble(iter)) )
      ! making the nb range quite small
      if (nebrestrend) call GraphConstraints(cx, gdsrestspring*gg, nbstrength*gg, 0.3*nbrange, kradius*gg)  ! Added *gg
      do j = 1, cx%na
        if (.not. cx%fixedatom(j)) then
          do k = 1, 3
            if (.not. cx%fixeddof((j-1)*3+k)) F((j-1)*3+k) = -cx%dvdr(k, j)
          enddo
        endif
      enddo
      Fnorm = norm2(F)
      eorg(1) = eorg(2) ; eorg(2) = eorg(3) ;  eorg(3) = cx%vcalc
      forg(1) = forg(2) ; forg(2) = forg(3) ;  forg(3) = Fnorm
      !print*, 'ENERGIES = ', eorg(3), eorg(2), eorg(1)
      !print*, 'FNORMS  = ', forg(3), forg(2), forg(1)
      !print*, 'FORCE = ', iter,sum(forg)/3.0d0, abs(eorg(2)*2.0d0-eorg(1)-eorg(3)), fthresh, gg
      if (iter .gt. 3 .and. (sum(forg)/3.0d0 .lt. fthresh .or. abs(eorg(2)*2.0d0-eorg(1)-eorg(3)) .lt. 1.0E-7)) exit
      !if ( sum(forg)/3.0d0 .lt. fthresh  ) exit
    enddo
    iterin = iterin + iter
    !print*,'FINISHED MINIM'

    return
  end Subroutine FireMinimise


  subroutine mincxsenergynebrest(cx, success)
    type(cxs)        :: cx
    logical          :: success

    success = .true.
    call GraphConstraints(cx, gdsrestspring, nbstrength, nbrange, kradius)  ! Added *gg
    !Call GraphConstraints( cx, gdsrestspring * 0.50d0, nbstrength *0.50d0 , 0.50d0 * nbrange, kradius )  ! Added *gg
    call AbInitio(cx, 'ener', success)

    return
  end subroutine

  
  subroutine mincxsenergy(cx,success)
    type(cxs)        :: cx
    logical          :: success

    success = .true.
    call AbInitio(cx, 'ener', success)

    return
  end subroutine

  !
  !*************************************************************************
  !> MatchCXContraints
  !!
  !! provides a force proportional to the difference btween the coordinates of the two structures cx and mcx
  !!
  !! - cx: A chemical structure object.
  !! - kspring: Spring constant for the bonding restraint term.
  !!
  !*************************************************************************
  !
  Subroutine MatchCXConstraint(cx, mcx, kspring)
    implicit none
    type(cxs) :: cx, mcx
    integer :: i, j, na
    real(8) :: kspring

    na = cx%na
    do i = 1, na
      do j = 1, 3
        cx%dvdr(j, i) = cx%dvdr(j, i) + 2.0d0*kspring*(mcx%r(j, i)-cx%r(j, i))
      enddo
    enddo

    return
  end Subroutine MatchCXConstraint


  !
  !*************************************************************************
  !
  !> GetActMolecularEnergies
  !!
  !! Get the energies of each molecule in the reactants and produts of rp for molecules that are "active"
  !!
  !! - cx: chemical structure.
  !! - MolecAct: array with the molecue indecies
  !! - nactmol: number of molecules in MolecAct array
  !! - molen: array with energies
  !! - filenam: The filename to be outputed
  !!
  !*************************************************************************
  !
  Subroutine GetActMolecEnergies(cx, MolecAct, nactmol, molen)
    implicit none
    integer :: i, j, k, l, m, n, nd, nb, naa, imove, ic, nactmol
    type(cxs) :: cx, mcx(nactmol)
    integer :: MolecAct(nactmol), nrx
    double precision, optional :: molen(nactmol)
    logical :: ldum

    do i = 1, nactmol  ! for active molecules in fron and tail
      M = MolecAct(i)  ! index of molecule in larger set
      call CreateMolecularCX(cx, mcx(i), m)
      ! calculate energy
      call AbInitio(mcx(i), 'grad', ldum)
      if (present(molen)) molen(i) = mcx(i)%vcalc
      cx%molen(m) = mcx(i)%vcalc
    enddo

    return
  end Subroutine GetActMolecEnergies


  !
  !*************************************************************************
  !
  !> GetMolecularEnergies
  !!
  !! Get the energies of each molecule in the reactants and produts of rp
  !!
  !! - cx: chemical structure.
  !! - filenam: The filename to be outputed
  !!
  !*************************************************************************
  !
  Subroutine GetMolecularEnergies(cx)
    implicit none
    integer :: i, j, k, l, m, n, nd, nb, naa, imove, ic, nactmol
    type(cxs) :: cx
    type(cxs) :: mcx(cx%nmol)
    logical :: ldum

    do i = 1, cx%nmol  ! for active molecules in fron and tail
      call CreateMolecularCX(cx, mcx(i), i)
      call AbInitio(mcx(i), 'grad', ldum)
      cx%molen(i) = mcx(i)%vcalc
    enddo

    return
  end Subroutine GetMolecularEnergies


end module pes
