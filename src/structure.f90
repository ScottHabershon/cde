!
!***************************************************************************************
!
!> @brief Chemical structure object definition.
!
!> Defines the chemical structure (cxs) object type.
!
!***************************************************************************************
!

Module chemstr

  use constants
  use functions
  use globaldata
  implicit none


  !> Chemical structure object definition.
  !
  type cxs
    sequence
    integer :: na
    real(8), dimension(:, :), allocatable :: r                   !< Coordinates in atomic units (bohr)
    real(8), dimension(:, :), allocatable :: p                   !< Momenta in atomic units (bohr / au)
    real(8), dimension(:, :), allocatable :: dvdr                !< PES derivatives wrt coordinates in atomic units (Eh/bohr)
    real(8), dimension(:, :), allocatable :: force               !< Forces wrt coordinates in atomic units (Eh/bohr)
    real(8), dimension(:, :), allocatable::  Hessian             !< Hessian matrix
    character (len=2), dimension(:), allocatable :: atomlabel    !< Atom labels
    real(8), dimension(:), allocatable :: mass                   !< Atomic masses in au.
    integer :: ndofconstr                                        !< Number of DOF constraints.
    integer :: natomconstr                                       !< Number of atom constraints.
    logical, dimension(:), allocatable :: fixeddof               !< Fixed DOF ids
    logical, dimension(:), allocatable :: fixedatom              !< Fixed atom ids
    real(8) :: vcalc                                             !< Calculate potential energy of structure.
    integer, dimension(:, :), allocatable :: graph               !< Bonding graph for structure.
    real(8), dimension(:), allocatable :: sprintcoords           !< SPRINT coordinates of the structure.
    integer, dimension(:), allocatable :: namol                  !< Atom numbers for each molecule in graph.
    integer, dimension(:, :), allocatable :: molid               !< Atom indices for each molecule in graph.
    integer :: nmol                                              !< Number of molecules defined by structure
    real(8) :: vcon                                              !< Constraint potential energy
    real(8) :: fitness                                           !< Fitness for molecular optimization
    real(8) :: fitness_scaled                                    !< Scaled fitness for molecular optimization.
    character(len=10) :: method = ''                             !< abinitio method used to calculate properties
    integer :: nbonds                                            !< Number of bonds
    integer :: nangles                                           !< Number of angles
    integer :: ntors                                             !< Number of torsions
    real(8), dimension(NBONDMAX):: bondl                         !< Calculated bond lengths
    real(8), dimension(NANGLEMAX) :: angle                       !< Calculated angles
    real(8), dimension(NTORSMAX) :: torsion                      !< Calculated torsion angles
    integer, dimension(NBONDMAX, 2) :: bondid                    !< IDs of atoms in bonds
    integer, dimension(NANGLEMAX, 3) :: angleid                  !< IDs of atoms in angles
    integer, dimension(NTORSMAX, 4) :: torsid                    !< IDs of atoms in torsions
    real(8), dimension(NBONDMAX,2,3) :: dbonddr                  !< Derivatives of bond-lengths wrt xyz of atoms.
    real(8), dimension(NANGLEMAX,3,3) :: dangdr                  !< Derivatives of angles wrt xyz of atoms.
    real(8), dimension(NTORSMAX,4,3) :: dtorsdr                  !< Derivatives of torsions wrt xyz of atoms.

    real(8), dimension(:), allocatable :: molen                  !< energy per molecule
    integer, dimension(:), allocatable :: molspin                !< spin on each molecule
    integer, dimension(:), allocatable :: molcharge              !< Total charge on each molecule

    integer :: itargetmol                                        ! Molecule ID for target molecule (igfunc=4)
  end type cxs

contains

  !
  !************************************************************************
  !> CreateCXS
  !!
  !! Creates a CXS object directly using labels and x,y,z coordinates which
  !! are passed in.
  !!
  !! cx - the created CXS object.
  !! na - number of atoms.
  !! label - atomic labels.
  !! x,y,z, - atomic coordinates in atomic units.
  !
  !************************************************************************
  !
  Subroutine CreateCXS(cx, na, label, x, y, z)
    type(cxs) :: cx
    integer :: na, namax, i, ierr
    character (len=2), dimension(:) :: label
    real(8), dimension(:) :: x, y, z
    ! Allocate space.
    !
    cx%na = na
    allocate(cx%r(3,na), stat = ierr)
    if (ierr /= 0) stop '* Memory allocation error in CreateCXS: r(3,na)'
    cx%r = 0.0d0

    allocate(cx%p(3,na), stat = ierr)
    if (ierr /= 0) stop '* Memory allocation error in CreateCXS: p(3,na)'
    cx%p = 0.0d0

    allocate(cx%dvdr(3,na), stat = ierr)
    if (ierr /= 0) stop '* Memory allocation error in CreateCXS: dvdr(3,na)'
    cx%dvdr = 0.0d0

    allocate(cx%force(3,na), stat = ierr)
    if (ierr /= 0) stop '* Memory allocation error in CreateCXS: forces(3,na)'
    cx%force = 0.0d0

    allocate(cx%Hessian(3*na,3*na), stat=ierr)
    if (ierr /= 0) stop '* Memory allocation error in CreateCXS: Hessian'

    allocate(cx%atomlabel(na), stat = ierr)
    if (ierr /= 0) stop '* Memory allocation error in CreateCXS: atomlabel(na)'
    cx%atomlabel = ''

    allocate(cx%mass(na), stat = ierr)
    if (ierr /= 0) stop '* Memory allocation error in CreateCXS: mass(na)'
    cx%mass = 0.0d0

    allocate(cx%fixeddof(3*na), stat = ierr)
    if (ierr /= 0) stop '* Memory allocation error in CreateCXS: fixeddof(3*na)'
    cx%fixeddof = .false.

    allocate(cx%fixedatom(na), stat = ierr)
    if (ierr /= 0) stop '* Memory allocation error in CreateCXS: fixedatom(na)'
    cx%fixedatom = .false.

    allocate(cx%graph(na,na), stat = ierr)
    if (ierr /= 0) stop '* Memory allocation error in CreateCXS: graph(na,na)'
    cx%graph = 0

    allocate(cx%sprintcoords(na), stat = ierr)
    if (ierr /= 0) stop '* Memory allocation error in CreateCXS: sprintcoords(na)'
    cx%sprintcoords = 0.0d0

    allocate(cx%molid(na,na), stat = ierr)
    if (ierr /= 0) stop '* Memory allocation error in CreateCXS: molid(na,na)'
    cx%molid = 0

    allocate(cx%namol(na), stat = ierr)
    if (ierr /= 0) stop '* Memory allocation error in CreateCXS: namol(na)'
    cx%namol = 0

    allocate(cx%molcharge(nmolmax), stat = ierr)
    if (ierr /= 0) stop '* Memory allocation error in CreateCXS: molcharge(na)'
    cx%molcharge = 1

    allocate(cx%molspin(nmolmax), stat = ierr)
    if (ierr /= 0) stop '* Memory allocation error in CreateCXS: molspin(na)'
    cx%molspin = 0

    allocate(cx%molen(nmolmax), stat = ierr)
    if (ierr /= 0) stop '* Memory allocation error in CreateCXS: molen(na)'
    cx%molen = 0.0d0

    ! Allocate coordinates to r(3,na).
    !
    do i = 1, na
      cx%atomlabel(i) = label(i)
      cx%r(1, i) = x(i)
      cx%r(2, i) = y(i)
      cx%r(3, i) = z(i)
    enddo
    cx%Hessian(:, :) = 0.d0

    ! Allocate masses.
    !
    call SetMass(cx)

    return
  end subroutine CreateCXS


  !
  !************************************************************************
  !> DeleteCXS
  !!
  !! Deletes memory associated with a CXS object.
  !!
  !! cx - the created CXS object.
  !
  !************************************************************************
  !
  Subroutine DeleteCXS(cx)
    type(cxs) :: cx

    ! Deallocate space.
    !
    deallocate(cx%r)
    deallocate(cx%p)
    deallocate(cx%dvdr)
    deallocate(cx%force)
    deallocate(cx%Hessian)
    deallocate(cx%atomlabel)
    deallocate(cx%mass)
    deallocate(cx%fixeddof)
    deallocate(cx%fixedatom)
    deallocate(cx%graph)
    deallocate(cx%sprintcoords)
    deallocate(cx%molid)
    deallocate(cx%namol)
    deallocate(cx%molcharge)
    deallocate(cx%molspin)
    deallocate(cx%molen)

    return
  end subroutine DeleteCXS


  !
  !************************************************************************
  !> CopyToNewCXS
  !!
  !! Creates a new copy of cx1 and puts it in cx2. Note that this routine
  !! allocates new space for cx2 - cx2 must not have already been
  !! allocated.
  !!
  !! - cx1: Input chemical structure object.
  !! - cx2: Nex chemical structure object.
  !!
  !************************************************************************
  !
  Subroutine CopytoNewCXS(cx1, cx2)
    implicit none
    type(cxs) :: cx1, cx2
    integer :: i, na
    real(8), dimension(:), allocatable :: x, y, z
    character (len=2), dimension(:), allocatable :: label

    na = cx1%na
    allocate(x(na), y(na), z(na), label(na))
    do i = 1, na
      label(i) = cx1%atomlabel(i)
      x(i) = cx1%r(1, i)
      y(i) = cx1%r(2, i)
      z(i) = cx1%r(3, i)
    enddo
    call CreateCXS(cx2, na, label, x, y, z)
    deallocate(x, y, z, label)

    ! Copy the rest of the variables from cx1 to cx2.
    !
    cx2%p = cx1%p
    cx2%dvdr = cx1%dvdr
    cx2%force = cx1%force
    cx2%mass = cx1%mass
    cx2%hessian = cx1%hessian
    cx2%ndofconstr = cx1%ndofconstr
    cx2%natomconstr = cx1%natomconstr
    cx2%fixeddof = cx1%fixeddof
    cx2%fixedatom = cx1%fixedatom
    cx2%vcalc = cx1%vcalc
    cx2%graph = cx1%graph
    cx2%sprintcoords = cx1%sprintcoords
    cx2%molid = cx1%molid
    cx2%namol = cx1%namol
    cx2%nmol = cx1%nmol
    cx2%vcon = cx1%vcon
    cx2%fitness = cx1%fitness
    cx2%fitness_scaled = cx1%fitness_scaled
    cx2%molcharge = cx1%molcharge
    cx2%molspin = cx1%molspin
    cx2%molen = cx1%molen
    cx2%itargetmol = cx1%itargetmol
    cx2%hessian = cx1%hessian

    return
  end Subroutine CopytoNewCXS


  !
  !************************************************************************
  !> CopyCXS
  !!
  !! Creates a copy of cx1 and in cx2. Note that cx2 must already
  !! exist as a cxs object.
  !!
  !! - cx1: Input chemical structure object.
  !! - cx2: Nex chemical structure object.
  !!
  !************************************************************************
  !
  Subroutine CopyCXS(cx1, cx2)
    implicit none
    type(cxs) :: cx1, cx2
    integer :: i, na

    na = cx1%na
    cx2%na = na
    do i = 1, na
      cx2%atomlabel(i) = cx1%atomlabel(i)
      cx2%r(1, i) = cx1%r(1, i)
      cx2%r(2, i) = cx1%r(2, i)
      cx2%r(3, i) = cx1%r(3, i)
    enddo

    cx2%p = cx1%p
    cx2%dvdr = cx1%dvdr
    cx2%force = cx1%force
    cx2%hessian = cx1%hessian
    cx2%mass = cx1%mass
    cx2%ndofconstr = cx1%ndofconstr
    cx2%natomconstr = cx1%natomconstr
    cx2%fixeddof = cx1%fixeddof
    cx2%fixedatom = cx1%fixedatom
    cx2%vcalc = cx1%vcalc
    cx2%graph = cx1%graph
    cx2%sprintcoords = cx1%sprintcoords
    cx2%molid = cx1%molid
    cx2%namol = cx1%namol
    cx2%nmol = cx1%nmol
    cx2%vcon = cx1%vcon
    cx2%fitness = cx1%fitness
    cx2%fitness_scaled = cx1%fitness_scaled
    cx2%molcharge = cx1%molcharge
    cx2%molspin = cx1%molspin
    cx2%molen = cx1%molen

    return
  end Subroutine CopyCXS


  !
  !************************************************************************
  !> SetMass
  !!
  !! Sets the atomic masses in a cxs object based on atomic labels.
  !! The mass values are stored in 'constants.f90'.
  !!
  !! - cx: The chemical structure object.
  !************************************************************************
  !
  Subroutine SetMass(cx)
    type (cxs) :: cx
    integer :: i, id
    character (len=2) :: label

    do i = 1, cx%na
      label = cx%atomlabel(i)
      if (len(trim(adjustl(label))) .ne. 0) then
        id = LabelToNumber(label)
        cx%mass(i) = MASS(id)
      else
        cx%mass(i) = 0.0d0
      endif
    enddo

    return
  end Subroutine SetMass


  !
  !************************************************************************
  !> CreateCXSFromXYZ
  !!
  !! Creates empty CXS from XYZ
  !!
  !! - cx: the created CXS object.
  !! - ifile: the input file to read the coordinates from. Note that this must be
  !!         an xyz file format.
  !!
  !************************************************************************
  !
  Subroutine CreateCXSFromXYZ(cx, ifile)
    character (len=*), intent(in) :: ifile
    character (len=100) :: comment
    character (len=2), dimension(NAMAX) :: label
    real(8), dimension(NAMAX) :: x, y, z
    integer :: na, ierr, ios, i
    type(cxs) :: cx
    logical :: there

    ! Check file existence.
    !
    inquire( file = ifile, exist = there )
    if (.not. there) stop '* ERROR in ReadCXS in structure.f90: specified input file does not exist'

    ! Open and read file.
    !
    open(10, file = ifile, status = 'unknown')
    read(10, *, iostat=ios) na
    read(10, '(A)') comment
    do i = 1, na
      read(10, *, iostat=ios) label(i), x(i), y(i), z(i)
      x(i) = x(i) * ang_to_bohr
      y(i) = y(i) * ang_to_bohr
      z(i) = z(i) * ang_to_bohr
    enddo
    call CreateCXS(cx, na, label, x, y, z)
    close(10)

    return
 end Subroutine CreateCXSFromXYZ


  !
  !************************************************************************
  !> ReadXYZtoCXS
  !!
  !! Reads in a set of coordinates from xyz and puts it in CXS
  !! beware: the numer of atoms and indecies must be the same in the XYZ as the
  !! CXS  allocated arrays
  !!
  !! - cx: the created CXS object.
  !! - ifile: the input file to read the coordinates from. Note that this must be
  !!         an xyz file format.
  !!seb
  !************************************************************************
  !
  Subroutine ReadXYZtoCXS(cx, ifile)
    character (len=*), intent(in) :: ifile
    character (len=100) :: comment
    character (len=2) :: label
    real(8) :: x, y, z
    integer :: na, ierr, ios, i
    type(cxs) :: cx
    logical :: there

    ! Check file existence.
    !
    inquire(file = ifile, exist = there )
    if (.not. there) stop '* ERROR in ReadCXS in structure.f90: specified input file does not exist'

    ! Open and read file.
    !
    Open(10, file = ifile, status = 'unknown')
    read(10, *, iostat=ios) na
    if (na .ne. cx%na) then
      print*, 'Number of atoms read from file'//trim(adjustl(ifile))//&
      ' does not match allocated array size in CXS (ReadXYZtoCXS in structure.f90)'
      stop
    endif
    if (ios /= 0) then
      print *, '* Error in ReadCXS: Odd number of atoms in input file - ', ifile
      stop
    endif
    read(10, '(A)') comment
    cx%r = 0.0d0
    ! Read in the atoms.
    !
    do i = 1, na
      read(10, *, iostat=ios) label, x, y, z
      if (ios /= 0)then
        stop '* ERROR reading atom list from input file in ReadCXS'
      endif
      cx%atomlabel(i) = label
      cx%r(1, i) = x * ang_to_bohr
      cx%r(2, i) = y * ang_to_bohr
      cx%r(3, i) = z * ang_to_bohr
    enddo
    close(10)

    return
  end Subroutine ReadXYZtoCXS


  !
  !************************************************************************
  !> ReadCXS
  !!
  !! Reads in a set of coordinates to a CXS object from an xyz file.
  !!
  !! - cx: the created CXS object.
  !! - ifile: the input file to read the coordinates from. Note that this must be
  !!         an xyz file format.
  !!
  !************************************************************************
  !
  Subroutine ReadCXS(cx, ifile)
    character (len=*), intent(in) :: ifile
    character (len=100) :: comment
    character (len=2) :: label
    real(8) :: x, y, z
    integer :: na, ierr, ios, i
    type(cxs) :: cx
    logical :: there

    ! Check file existence.
    !
    inquire(file = ifile, exist = there)
    if (.not. there) stop '* ERROR in ReadCXS in structure.f90: specified input file does not exist'

    ! Open and read file.
    !
    Open(10, file = ifile, status = 'unknown')
    read(10, *, iostat=ios)na
    if (ios /= 0) then
      print *, '* Error in ReadCXS: Odd number of atoms in input file - ', ifile
      stop
    endif
    read(10, '(A)') comment

    ! Assign workspace for na atoms.
    !
    cx%na = na
    allocate( cx%r(3,na), stat = ierr )
    if (ierr /= 0) stop '* Memory allocation error in ReadCXS: r(3,na)'

    allocate( cx%p(3,na), stat = ierr )
    if (ierr /= 0) stop '* Memory allocation error in ReadCXS: p(3,na)'

    allocate( cx%dvdr(3,na), stat = ierr )
    if (ierr /= 0) stop '* Memory allocation error in ReadCXS: dvdr(3,na)'

    allocate( cx%force(3,na), stat = ierr )
    if (ierr /= 0) stop '* Memory allocation error in ReadCXS: forces(3,na)'

    allocate( cx%hessian(3*na,3*na))

    allocate( cx%atomlabel(na), stat = ierr )
    if (ierr /= 0) stop '* Memory allocation error in ReadCXS: atomlabel(na)'

    allocate( cx%mass(na), stat = ierr )
    if (ierr /= 0) stop '* Memory allocation error in ReadCXS: mass(na)'

    allocate( cx%fixeddof(3*na), stat = ierr )
    if (ierr /= 0) stop '* Memory allocation error in ReadCXS: fixeddof(3*na)'

    allocate( cx%fixedatom(na), stat = ierr )
    if (ierr /= 0) stop '* Memory allocation error in ReadCXS: fixedatom(na)'

    allocate( cx%graph(na,na), stat = ierr )
    if (ierr /= 0) stop '* Memory allocation error in ReadCXS: graph(na,na)'

    allocate( cx%sprintcoords(na), stat = ierr )
    if (ierr /= 0) stop '* Memory allocation error in ReadCXS: sprintcoords(na)'

    allocate( cx%molid(na,na), stat = ierr )
    if (ierr /= 0) stop '* Memory allocation error in ReadCXS: molid(na,na)'

    allocate( cx%namol(na), stat = ierr )
    if (ierr /= 0) stop '* Memory allocation error in ReadCXS: namol(na)'

    allocate( cx%molcharge(nmolmax), stat = ierr )
    if (ierr /= 0) stop '* Memory allocation error in ReadCXS: molcharge(nmolmax)'

    allocate( cx%molspin(nmolmax), stat = ierr )
    If (ierr /= 0) stop '* Memory allocation error in ReadCXS: molspin(nmolmax)'

    allocate( cx%molen(na), stat = ierr )
    if (ierr /= 0) stop '* Memory allocation error in ReadCXS: molen(na)'

    ! Read in the atoms.
    !
    do i = 1, na
      read(10, *, iostat=ios) label, x, y, z
      if (ios /= 0) then
        stop '* ERROR reading atom list from input file in ReadCXS'
      endif
      cx%atomlabel(i) = label
      cx%r(1, i) = x * ang_to_bohr
      cx%r(2, i) = y * ang_to_bohr
      cx%r(3, i) = z * ang_to_bohr
    enddo
    close(10)

    return
  end Subroutine ReadCXS


  !
  !************************************************************************
  !> SetCXSConstraints
  !!
  !! sets constraints on DOFs and atoms in a structure.
  !!
  !! - cx: the created CXS object.
  !! - NDOFconstr: Number of DOF constraints.
  !! - FixedDOF: Array containing integer number of constrained DOFs
  !! - Natomconstr: Number of atom constraints
  !! - Fixedatom: Array containing integer ids of fixed atoms.
  !!
  !************************************************************************
  !
  Subroutine SetCXSconstraints(cx, NDOFconstr, FixedDOF, Natomconstr, FixedAtom)
    implicit none
    type(cxs) :: cx
    integer :: i, j
    integer, intent(in) :: NDOFconstr, Natomconstr
    integer, dimension(:), intent(in) :: FixedDOF, FixedAtom

    cx%fixedatom(:) = .false.
    cx%fixeddof(:) = .false.
    cx%ndofconstr = ndofconstr
    cx%natomconstr = natomconstr

    do i = 1, cx%NDOFconstr
      j = FixedDOF(i)
      cx%fixeddof(j) = .true.
    enddo

    do i = 1, cx%Natomconstr
      j = FixedAtom(i)
      cx%fixedatom(j) = .true.
      cx%fixeddof(3*j-2) = .true.
      cx%fixeddof(3*j-1) = .true.
      cx%fixeddof(3*j) = .true.
    enddo

    return
  end Subroutine SetCXSconstraints


  !
  !************************************************************************
  !> GetProjectedMomenta
  !!
  !! Calculates the momenta projected along the direction of the
  !! force array - this is usually the forces from (CI)NEB.
  !!
  !! - cx: A chemical structure object.
  !!
  !************************************************************************
  !
  Subroutine GetProjectedMomenta(cx)
    implicit none
    type(cxs) :: cx
    integer :: i, j, idof, k
    real(8) :: dot, norm, nfac, nfac2

    ! Get the normalization for the force.
    !
    norm = 0.0
    idof = 0
    do i = 1, cx%na
      if (.not. cx%Fixedatom(i)) then
        do k = 1, 3
          idof = idof + 1
          if (.not. cx%fixeddof(idof)) then
            norm = norm + cx%force(k, i) * cx%force(k, i)
          endif
        enddo
      else
        idof = idof + 3
      endif
    enddo
  !  if (sqrt(norm) .gt. epsil) then
     nfac = 1.0 / sqrt(norm)
  !  else
  !   nfac = 0.0d0
  !  endif


    ! Calculate dot product.
    !
    dot = 0.0
    idof = 0
    do i = 1, cx%na
      ! print*, 'MASS = ', cx%mass(i), cx%atomlabel(i)
      if (.not. cx%fixedatom(i)) then
        do k = 1, 3
          idof = idof + 1
          if (.not. cx%fixeddof(idof)) then
            dot = dot + (cx%p(k, i)/cx%mass(i)) * cx%force(k, i) * nfac
          !    dot = dot + (cx%p(k,i)) * cx%force(k,i) * nfac
          endif
        enddo
      else
        idof = idof + 3
      endif
    enddo

   ! print*,'DOT = ',dot

    ! Zero the velocity if it is pointing in the opposite direction to force.
    !
    if (dot < 0.0) then
      cx%p(:, :) = 0.d0
    else
      idof = 0
      do i = 1, cx%na
        if (.not. cx%fixedatom(i)) then
          do k = 1, 3
            idof = idof + 1
            if (.not. cx%fixeddof(idof)) then
              cx%p(k, i) = dot * cx%force(k, i) * nfac * cx%mass(i)
            endif
          enddo
        else
          idof = idof + 3
        endif
      enddo
    endif

    return
  end Subroutine GetProjectedMomenta


  !
  !************************************************************************
  !> GetGraph
  !!
  !! Calculates the bonding graph for the system. The values which define
  !! whether or not atoms are bonded are based on covalent radii, and
  !! are stored in constants.f90 in an atomic-number-labelled array
  !! called BondingCutoff(:,:).
  !!
  !! - cx: A chemical structure object.
  !!
  !************************************************************************
  !
  Subroutine GetGraph(cx)
    type(cxs) :: cx
    integer :: id1, id2, i, j
    real(8) :: dx, dy, dz
    real(8) :: rsq, rr, rcut

    do i = 1, cx%na
      cx%graph(i,i) = 0
      do j = i + 1, cx%na
        dx = cx%r(1, i) - cx%r(1, j)
        dz = cx%r(2, i) - cx%r(2, j)
        dy = cx%r(3, i) - cx%r(3, j)
        rsq = dx*dx + dy*dy + dz*dz
        rr = sqrt(rsq)

        ! Turn atomic labels into integer values.
        !
        id1 = LabelToNumber(cx%atomlabel(i))
        id2 = LabelToNumber(cx%atomlabel(j))

        ! Get cutoff.
        !
        Rcut = (CovRad(id1) + covrad(id2)) * BONDINGSF

        ! Check that this pair-type have actually been included in
        ! constants.f90
        !
        if (CovRad(id1) < 1d-3 .or. CovRad(id2) < 1d-3) then
          print *
          print *, '* Error: undefined BondingCutoff in structure.f90 for', cx%atomlabel(i), cx%atomlabel(j)
          print *
          stop
        endif

        ! Assess whether bonded or not.
        !
        if (rr <= rcut) then
          cx%graph(i, j) = 1
          cx%graph(j, i) = 1
          !print*,'Binding: ',  i,j,cx%atomlabel(i),cx%atomlabel(j),rr,rcut,CovRad(id1),covrad(id2)
        else
          cx%graph(i, j) = 0
          cx%graph(j, i) = 0
          !PRint*,'NON-Binding: ',  i,j,cx%atomlabel(i),cx%atomlabel(j),rr,rcut,CovRad(id1),covrad(id2)
        endif
      enddo
    enddo

    return
  end Subroutine GetGraph


  !
  !********************************************************************************
  !> GetMols
  !!
  !! Uses the connectivity graph to determine the molecules defined by the
  !! chemical structure cx. To do so, we use Floyd-Warshall shortest-path walks
  !! on the chemical structure graph. Molen/Molcharge arrays will be reset, so beware!
  !!
  !! - cx: Input chemical structure object.
  !
  !********************************************************************************
  !
  Subroutine GetMols( cx )
    implicit none
    type(cxs) :: cx
    integer :: na, i, j, k, l
    real(8), dimension(:, :), allocatable :: dist, dsp
    integer :: ifinder
    integer, dimension(NAMAX) :: ifound

    ! Allocate distance workspaces.
    na = cx%na

    allocate(dist(na, na), dsp(na, na))

    ! Set flags determining whether each atom
    ifound(:) = 0

    ! Zero molecule details.
    !
    cx%namol(:) = 0
    cx%molid(:, :) = 0
    cx%nmol = 0
    cx%molen(:) = 0

    ! Calculate a distance matrix - the distance is defined here
    ! as 1 if two atoms are bonded. Otherwise it is assigned the
    ! value BIG (from constants.f90).
    !
    do i = 1, na
      do j = 1, na
        if (cx%graph(i, j) == 1) then
          dist(i, j) = 1.d0
          dist(j, i) = 1.d0
        else
          dist(i, j) = BIG
          dist(j, i) = BIG
        endif
      enddo
    enddo

    ! Get shortest paths.
    !
    call GetShortestPaths(na, dist, dsp)

    ! Identify molecules.
    !
    ifinder = 1
    do while (ifinder == 1)

      inner: do i = 1, na
        if (ifound(i) == 0) then
          cx%nmol = cx%nmol + 1
          ifound(i) = 1
          cx%namol(cx%nmol) = 1
          cx%molid(cx%nmol, 1) = i
          exit inner
        endif
      enddo inner

      outer: do j = 1, na
        if (ifound(j) == 0) then
          do k = 1, cx%nmol
            do l = 1, cx%namol(k)
              if (dsp(j, cx%molid(k, l)) < BIG-1) then
                cx%namol(k) = cx%namol(k) + 1
                cx%molid(k, cx%namol(k)) = j
                ifound(j) = 1
                cycle outer
              endif
            enddo
          enddo
        endif
      enddo outer

      ifinder = 0
      do i = 1, na
        if (ifound(i) == 0) then
          ifinder = 1
        endif
      enddo

    enddo

    deallocate( dist, dsp )
    ! molen and molcharge arrays will be set to zero, so beware!
    if (allocated(cx%molen) ) then
      deallocate(cx%molen)
      allocate(cx%molen(cx%nmol))
      cx%molen = 0.0d0
    endif
    if (allocated(cx%molcharge) ) then
      deallocate(cx%molcharge)
      allocate(cx%molcharge(cx%nmol))
      cx%molcharge = 0
    endif
    if (allocated(cx%molspin) ) then
      deallocate(cx%molspin)
      allocate(cx%molspin(cx%nmol))
      cx%molspin = 1
    endif

    return
  end Subroutine GetMols


  !
  !===================================================================================
  !> GetShortestPaths
  !!
  !! Calculates the matrix dsp(i,j) which contains the shortest path between points i
  !! and j as determined by Floyd-Warshall algorithm.
  !!
  !! - n:  number of points
  !! - dg(:,:):  distance between points.
  !! - dsp(:,:): Shortest-path between points.
  !!
  !===================================================================================
  !
  Subroutine GetShortestPaths(N,dg,dsp)
    implicit none
    integer :: N, i, j, k
    real(8), dimension(:, :) :: dg, dsp
    real(8), dimension(:, :), allocatable :: dgstore

    allocate(dgstore(N, N))
    dgstore(:, :) = dg(:, :)

    do k = 1, N
      do i = 1, N
        do j = 1, N
          dg(i, j) = min(dg(i, j), dg(i, k) + dg(j, k))
        enddo
      enddo
    enddo

    dsp(:, :) = dg(:, :)
    dg(:, :) = dgstore(:, :)
    deallocate(dgstore)

    return
  end Subroutine GetShortestPaths


  !
  !*************************************************************************
  !> GraphConstraints
  !!
  !! Calculates the potential energy and forces arising due to graph
  !! constraints.
  !!
  !! - cx: A chemical structure object.
  !! - kspring: Spring constant for the bonding restraint term.
  !! - nbstrength: Repulsive exponential interaction strength in au for
  !!               atom pairs with g(i,j) = 0
  !! - nbrange: Repulsive exponential interaction range in au for
  !!            atom pairs with g(i,j) = 0
  !! - kradius: Harmonic restraint term strength for different molecules.
  !!
  !*************************************************************************
  !
  Subroutine GraphConstraints(cx, kspring, nbstrength, nbrange, kradius)
    implicit none
    type(cxs) :: cx
    integer :: i, j, na, f, ii
    integer :: jj, i1, i2, j1, id1, id2
    real(8) :: kspring, rr, rsq, t1, dx, dy, dz, kradius
    real(8) :: nbstrength, nbrange, dr, eterm, onr, rx
    real(8), dimension(:, :), allocatable :: rmin, rmax
    real(8) :: M1, M2
    real(8), dimension(3) :: xcom1, xcom2

    ! Set local variables.
    !
    na = cx%na

    ! Allocate and assign the min / max constraints.
    !
    allocate(rmin(na, na), rmax(na, na))
    do i = 1, na
      do j = i, na
        id1 = LabelToNumber(cx%atomlabel(i))
        id2 = LabelToNumber(cx%atomlabel(j))
        rx = BONDINGSF * (COVRAD(id1) + COVRAD(id2))
        rmin(i, j) = rx - BONDINGRANGE1
        rmax(i, j) = rx + BONDINGRANGE2
        rmin(j, i) = rmin(i, j)
        rmax(j, i) = rmax(i, j)
      enddo
    enddo

    ! Calculate constraints based on current graph.
    !
    cx%vcon = 0.d0
    do i = 1, na-1
      do j = i+1, na
        dx = cx%r(1, i) - cx%r(1, j)
        dy = cx%r(2, i) - cx%r(2, j)
        dz = cx%r(3, i) - cx%r(3, j)
        !print*, dx,dy,dz
        rsq = dx*dx + dy*dy + dz*dz
        rr = sqrt(rsq)

        if (cx%graph(i, j) == 1) then    !! BONDED

          if (rr < rmin(i, j)) then
            !print*, i,j,cx%atomlabel(i)//' AND '//cx%atomlabel(j)//' apparently are too close!'//&
            !' dis ', rr, rmin(i,j)
            dr = rr - rmin(i, j)
            cx%vcon = cx%vcon + kspring * dr**2
            t1 = 2.0 * kspring * dr

            call AccumulateDerivatives(cx, t1, i, j)

          else if (rr > rmax(i, j)) then
            !print*, i,i,cx%atomlabel(i)//' AND '//cx%atomlabel(j)//' apparently are too far!'//&
            !' dis ', rr, rmax(i,j)
            dr = rr - rmax(i, j)
            cx%vcon = cx%vcon + kspring * dr**2
            t1 = 2.0 * kspring * dr

            call AccumulateDerivatives(cx, t1, i, j)

          endif
!!$
!!$
!!$                dr = rr - (0.5*(rmin(i,j) + rmax(i,j)))
!!$                cx%vcon = cx%vcon + kspring * dr**2
!!$                t1 = 2.0 * kspring * dr
!!$                Call AccumulateDerivatives(cx, t1, i, j)

        else if (cx%graph(i, j) == 0) then     !! NOT BONDED

          ! Exponential version...
          !
!             eterm = exp(-rr / nbrange)
!             cx%vcon = cx%vcon + nbstrength * eterm
!             t1 = nbstrength * eterm * (-1.d0/nbrange)

          ! Gaussian version...
          !
          eterm = exp(-(rr*rr) / (2.0 * nbrange**2))
          cx%vcon = cx%vcon + nbstrength * eterm
          t1 = nbstrength * (-2.0*rr / (2.0 * nbrange**2)) * eterm
          !print*, cx%atomlabel(i)//' AND '//cx%atomlabel(j)//'NON BONDING STRENGTH = ',  eterm

          call AccumulateDerivatives(cx, t1, i, j)

        endif
      enddo
    enddo

    !!
    !! Repulsion between SEPARATE molecules...
    !!
    !if ( cx%nmol > 1 ) Call COMGraphConstraints( cx, kradius )
    !
    ! Repulsion between SEPARATE molecules...
    !
    if (cx%nmol > 1) then

      ! COM version.
      !
      do i = 1, cx%nmol-1
        do j = i+1, cx%nmol

          M1 = 0.d0
          xcom1(:) = 0.d0
          do i1 = 1, cx%namol(i)
            ii = cx%molid(i, i1)
!                M1 = M1 + cx%mass(ii)
            M1 = M1 + 1.d0
            xcom1(1) = xcom1(1) + cx%r(1, ii) !* cx%mass(ii)
            xcom1(2) = xcom1(2) + cx%r(2, ii) !* cx%mass(ii)
            xcom1(3) = xcom1(3) + cx%r(3, ii) !* cx%mass(ii)
          enddo
          xcom1(:) = xcom1(:) / M1

          M2 = 0.d0
          xcom2(:) = 0.d0
          do i1 = 1, cx%namol(j)
            ii = cx%molid(j, i1)
!                M2 = M2 + cx%mass(ii)
            M2 = M2 + 1.d0
            xcom2(1) = xcom2(1) + cx%r(1, ii) !* cx%mass(ii)
            xcom2(2) = xcom2(2) + cx%r(2, ii) !* cx%mass(ii)
            xcom2(3) = xcom2(3) + cx%r(3, ii) !* cx%mass(ii)
          enddo
          xcom2(:) = xcom2(:) / M2

          dx = xcom1(1) - xcom2(1)
          dy = xcom1(2) - xcom2(2)
          dz = xcom1(3) - xcom2(3)
          rsq = dx*dx + dy*dy + dz*dz
          rr = sqrt(rsq)

!!$             if (rr > RADIUS_MAX) then
!!$
!!$                onr = 1.d0 / rr
!!$                dr = rr - RADIUS_MAX
!!$                cx%vcon = cx%vcon + kradius * dr**2
!!$                t1 = 2.0 * kradius * dr
!!$
!!$                do i1 = 1, cx%namol(i)
!!$                   ii = cx%molid(i,i1)
!!$!                   cx%dvdr(1,ii) = cx%dvdr(1,ii) + t1 * dx * onr * cx%mass(ii)/M1
!!$!                   cx%dvdr(2,ii) = cx%dvdr(2,ii) + t1 * dy * onr * cx%mass(ii)/M1
!!$!                   cx%dvdr(3,ii) = cx%dvdr(3,ii) + t1 * dz * onr * cx%mass(ii)/M1
!!$                   cx%dvdr(1,ii) = cx%dvdr(1,ii) + t1 * dx * onr /M1
!!$                   cx%dvdr(2,ii) = cx%dvdr(2,ii) + t1 * dy * onr /M1
!!$                   cx%dvdr(3,ii) = cx%dvdr(3,ii) + t1 * dz * onr /M1
!!$                enddo
!!$
!!$                do i1 = 1, cx%namol(j)
!!$                   jj = cx%molid(j,i1)
!!$!                   cx%dvdr(1,jj) = cx%dvdr(1,jj) - t1 * dx * onr* cx%mass(jj)/M2
!!$!                   cx%dvdr(2,jj) = cx%dvdr(2,jj) - t1 * dy * onr* cx%mass(jj)/M2
!!$!                   cx%dvdr(3,jj) = cx%dvdr(3,jj) - t1 * dz * onr* cx%mass(jj)/M2
!!$
!!$                   cx%dvdr(1,jj) = cx%dvdr(1,jj) - t1 * dx * onr/M2
!!$                   cx%dvdr(2,jj) = cx%dvdr(2,jj) - t1 * dy * onr/M2
!!$                   cx%dvdr(3,jj) = cx%dvdr(3,jj) - t1 * dz * onr/M2
!!$                enddo

!!$             else if (rr < RADIUS_MIN) then

          if (rr < RADIUS_MIN) then

            onr = 1.d0 / rr
            dr = rr - RADIUS_MIN
            cx%vcon = cx%vcon + kradius * dr**2
            t1 = 2.0 * kradius * dr

            do i1 = 1, cx%namol(i)
              ii = cx%molid(i, i1)
!                   cx%dvdr(1,ii) = cx%dvdr(1,ii) + t1 * dx * onr * cx%mass(ii)/M1
!                   cx%dvdr(2,ii) = cx%dvdr(2,ii) + t1 * dy * onr * cx%mass(ii)/M1
!                   cx%dvdr(3,ii) = cx%dvdr(3,ii) + t1 * dz * onr * cx%mass(ii)/M1
              cx%dvdr(1, ii) = cx%dvdr(1, ii) + t1 * dx * onr /M1
              cx%dvdr(2, ii) = cx%dvdr(2, ii) + t1 * dy * onr /M1
              cx%dvdr(3, ii) = cx%dvdr(3, ii) + t1 * dz * onr /M1
            enddo

            do i1 = 1, cx%namol(j)
              jj = cx%molid(j, i1)
!                   cx%dvdr(1,jj) = cx%dvdr(1,jj) - t1 * dx * onr* cx%mass(jj)/M2
!                   cx%dvdr(2,jj) = cx%dvdr(2,jj) - t1 * dy * onr* cx%mass(jj)/M2
!                   cx%dvdr(3,jj) = cx%dvdr(3,jj) - t1 * dz * onr* cx%mass(jj)/M2

              cx%dvdr(1, jj) = cx%dvdr(1, jj) - t1 * dx * onr /M2
              cx%dvdr(2, jj) = cx%dvdr(2, jj) - t1 * dy * onr /M2
              cx%dvdr(3, jj) = cx%dvdr(3, jj) - t1 * dz * onr /M2
            enddo

!!$                do i1 = 1, cx%namol(i)
!!$                   ii = cx%molid(i,i1)
!!$                   cx%dvdr(1,ii) = cx%dvdr(1,ii) + t1 * dx * onr * cx%mass(ii)/M1
!!$                   cx%dvdr(2,ii) = cx%dvdr(2,ii) + t1 * dy * onr * cx%mass(ii)/M1
!!$                   cx%dvdr(3,ii) = cx%dvdr(3,ii) + t1 * dz * onr * cx%mass(ii)/M1
!!$                enddo
!!$
!!$                do i1 = 1, cx%namol(j)
!!$                   jj = cx%molid(j,i1)
!!$                   cx%dvdr(1,jj) = cx%dvdr(1,jj) - t1 * dx * onr* cx%mass(jj)/M2
!!$                   cx%dvdr(2,jj) = cx%dvdr(2,jj) - t1 * dy * onr* cx%mass(jj)/M2
!!$                   cx%dvdr(3,jj) = cx%dvdr(3,jj) - t1 * dz * onr* cx%mass(jj)/M2
!!$                enddo

          endif

!!$             do i1 = 1, cx%namol(i)
!!$                ii = cx%molid(i,i1)
!!$                do j1 = 1, cx%namol(j)
!!$                   jj = cx%molid(j,j1)
!!$
!!$                   dx = cx%r(1,ii) - cx%r(1,jj)
!!$                   dy = cx%r(2,ii) - cx%r(2,jj)
!!$                   dz = cx%r(3,ii) - cx%r(3,jj)
!!$                   rsq = dx*dx + dy*dy + dz*dz
!!$                   rr = sqrt(rsq)
!!$
!!$                   if (rr < RADIUS_MIN) then
!!$                      onr = 1.d0 / rr
!!$                      dr = rr - RADIUS_MIN
!!$                      cx%vcon = cx%vcon + kradius * dr**2
!!$                      t1 = 2.0 * kradius * dr
!!$                      cx%dvdr(1,ii) = cx%dvdr(1,ii) + t1 * dx * onr
!!$                      cx%dvdr(2,ii) = cx%dvdr(2,ii) + t1 * dy * onr
!!$                      cx%dvdr(3,ii) = cx%dvdr(3,ii) + t1 * dz * onr
!!$                      cx%dvdr(1,jj) = cx%dvdr(1,jj) - t1 * dx * onr
!!$                      cx%dvdr(2,jj) = cx%dvdr(2,jj) - t1 * dy * onr
!!$                      cx%dvdr(3,jj) = cx%dvdr(3,jj) - t1 * dz * onr
!!$                   endif
!!$                enddo
!!$             enddo

        enddo
      enddo

    endif

    deallocate(rmin, rmax)

    return
  end Subroutine GraphConstraints


  !
  !*************************************************************************
  !> GraphConstraints_DoubleEnded
  !!
  !! Calculates the potential energy and forces arising due to graph
  !! constraints.
  !!
  !! - cx: A chemical structure object.
  !! - kspring: Spring constant for the bonding restraint term.
  !! - nbstrength: Repulsive exponential interaction strength in au for
  !!               atom pairs with g(i,j) = 0
  !! - nbrange: Repulsive exponential interaction range in au for
  !!            atom pairs with g(i,j) = 0
  !! - kradius: Harmonic restraint term strength for different molecules.
  !!
  !*************************************************************************
  !
  Subroutine GraphConstraints_DoubleEnded(cx, cxstart, kspring, nbstrength, nbrange, kradius)
    implicit none
    type(cxs) :: cx, cxstart
    integer :: i, j, na, f, ii
    integer :: jj, i1, i2, j1, id1, id2
    real(8) :: kspring, rr, rsq, t1 ,dx, dy, dz, kradius
    real(8) :: nbstrength, nbrange, dr, eterm, onr, rx
    real(8) :: dx0, dy0, dz0, rr0, rsq0, factor
    real(8), dimension(:, :), allocatable :: rmin, rmax
    real(8) :: M1, M2
    real(8), dimension(3) :: xcom1, xcom2

    ! Set local variables.
    !
    na = cx%na

    ! Allocate and assign the min / max constraints.
    !
    allocate(rmin(na,na), rmax(na,na))
    do i = 1, na
      do j = i, na
        id1 = LabelToNumber(cx%atomlabel(i))
        id2 = LabelToNumber(cx%atomlabel(j))
        !rx = ( COVRAD(id1) + COVRAD(id2) )
        rx = BONDINGSF * (COVRAD(id1) + COVRAD(id2))
        !rmin(i,j) = rx/2.0d0
        !rmax(i,j) = rx  - BONDINGRANGE1
        rmin(i, j) = rx - BONDINGRANGE1
        rmax(i, j) = rx + BONDINGRANGE2
        rmin(j, i) = rmin(i, j)
        rmax(j, i) = rmax(i, j)
      enddo
    enddo

    !
    ! Calculate constraints based on current graph of cx and the origin-point cxstart.
    !
    cx%vcon = 0.d0
    do i = 1, na-1
      do j = i+1, na

        ! Calculate distance in cx.
        !
        dx = cx%r(1, i) - cx%r(1, j)
        dy = cx%r(2, i) - cx%r(2, j)
        dz = cx%r(3, i) - cx%r(3, j)
        rsq = dx*dx + dy*dy + dz*dz
        rr = sqrt(rsq)

        ! Calculate distance in cxstart.
        !
        dx0 = cxstart%r(1, i) - cxstart%r(1, j)
        dy0 = cxstart%r(2, i) - cxstart%r(2, j)
        dz0 = cxstart%r(3, i) - cxstart%r(3, j)
        rsq0 = dx0*dx0 + dy0*dy0 + dz0*dz0
        rr0 = sqrt(rsq0)

        if (cx%graph(i, j) == 1 .and. cxstart%graph(i, j) == 0) then    !! BOND FORMATION !!

          if (rr < rmin(i, j)) then

            dr = rr - rmin(i, j)
            cx%vcon = cx%vcon + kspring * dr**2
            t1 = 2.0 * kspring * dr

            call AccumulateDerivatives(cx, t1, i, j)

          else if (rr > rmax(i, j) ) then
            dr = rr - rmax(i, j)
            cx%vcon = cx%vcon + kspring * dr**2
            t1 = 2.0 * kspring * dr

            call AccumulateDerivatives(cx, t1, i, j)

          endif

        else if (cx%graph(i, j) == 0 .and. cxstart%graph(i, j) == 1) then     !! BOND BREAKING !!

          ! Gaussian repulsion
          !
          eterm = exp(-(rr*rr) / (2.0 * nbrange**2))
          cx%vcon = cx%vcon + nbstrength * eterm
          t1 = nbstrength * (-2.0*rr / (2.0 * nbrange**2)) * eterm
          call AccumulateDerivatives(cx, t1, i, j)

        else if (cx%graph(i, j) == 1 .and. cxstart%graph(i, j) == 1) then

          factor = 5.0
          cx%vcon = cx%vcon + factor*kspring*(rr - rr0)*(rr-rr0)
          cx%dvdr(1, i) = cx%dvdr(1, i) + factor*2.0 * kspring * (rr-rr0) * (dx/rr)
          cx%dvdr(2, i) = cx%dvdr(2, i) + factor*2.0 * kspring * (rr-rr0) * (dy/rr)
          cx%dvdr(3, i) = cx%dvdr(3, i) + factor*2.0 * kspring * (rr-rr0) * (dz/rr)

          cx%dvdr(1, j) = cx%dvdr(1, j) + factor*2.0 * kspring * (rr-rr0) * (-dx/rr)
          cx%dvdr(2, j) = cx%dvdr(2, j) + factor*2.0 * kspring * (rr-rr0) * (-dy/rr)
          cx%dvdr(3, j) = cx%dvdr(3, j) + factor*2.0 * kspring * (rr-rr0) * (-dz/rr)

        else if (cx%graph(i, j) == 0 .and. cxstart%graph(i, j) == 0) then

          ! Gaussian repulsion with reduced range....
          !
          factor = 0.8
          eterm = exp(-(rr*rr) / (2.0 * (factor*nbrange)**2) )
          cx%vcon = cx%vcon + nbstrength * eterm
          t1 = nbstrength * (-2.0*rr / (2.0 * (factor*nbrange)**2)) * eterm
          call AccumulateDerivatives(cx, t1, i, j)

        endif
      enddo
    enddo

    !!
    !! Repulsion between SEPARATE molecules...
    !!
    !if ( cx%nmol > 1 ) Call COMGraphConstraints( cx, kradius )
    !
    ! Repulsion between SEPARATE molecules...
    !
    if (cx%nmol > 1) then
      ! COM version.
      !
      do i = 1, cx%nmol-1
        do j = i+1, cx%nmol

          M1 = 0.d0
          xcom1(:) = 0.d0
          do i1 = 1, cx%namol(i)
            ii = cx%molid(i, i1)
!                M1 = M1 + cx%mass(ii)
            M1 = M1 + 1.d0
            xcom1(1) = xcom1(1) + cx%r(1, ii) !* cx%mass(ii)
            xcom1(2) = xcom1(2) + cx%r(2, ii) !* cx%mass(ii)
            xcom1(3) = xcom1(3) + cx%r(3, ii) !* cx%mass(ii)
          enddo
          xcom1(:) = xcom1(:) / M1

          M2 = 0.d0
          xcom2(:) = 0.d0
          do i1 = 1, cx%namol(j)
            ii = cx%molid(j, i1)
!                M2 = M2 + cx%mass(ii)
            M2 = M2 + 1.d0
            xcom2(1) = xcom2(1) + cx%r(1, ii) !* cx%mass(ii)
            xcom2(2) = xcom2(2) + cx%r(2, ii) !* cx%mass(ii)
            xcom2(3) = xcom2(3) + cx%r(3, ii) !* cx%mass(ii)
          enddo
          xcom2(:) = xcom2(:) / M2

          dx = xcom1(1) - xcom2(1)
          dy = xcom1(2) - xcom2(2)
          dz = xcom1(3) - xcom2(3)
          rsq = dx*dx + dy*dy + dz*dz
          rr = sqrt(rsq)

          if (rr < RADIUS_MIN) then

            onr = 1.d0 / rr
            dr = rr - RADIUS_MIN
            cx%vcon = cx%vcon + kradius * dr**2
            t1 = 2.0 * kradius * dr

            do i1 = 1, cx%namol(i)
              ii = cx%molid(i, i1)
              cx%dvdr(1, ii) = cx%dvdr(1, ii) + t1 * dx * onr /M1
              cx%dvdr(2, ii) = cx%dvdr(2, ii) + t1 * dy * onr /M1
              cx%dvdr(3, ii) = cx%dvdr(3, ii) + t1 * dz * onr /M1
            enddo

            do i1 = 1, cx%namol(j)
              jj = cx%molid(j, i1)
              cx%dvdr(1, jj) = cx%dvdr(1, jj) - t1 * dx * onr /M2
              cx%dvdr(2, jj) = cx%dvdr(2, jj) - t1 * dy * onr /M2
              cx%dvdr(3, jj) = cx%dvdr(3, jj) - t1 * dz * onr /M2
            enddo
          endif
        enddo
      enddo
    endif

    deallocate( rmin, rmax )

    return
  end Subroutine GraphConstraints_DoubleEnded


  !
  !*************************************************************************
  !> ProjActMolRotTransDVDR2
  !!   projects the rotational and translational degrees of freedom of each nactmol active molecule in
  !!   molecact from the dcerivative vector dvdr. WARNING - THE PROJECTION OF
  !!   OF DVDR ONLY WORKS PARTIALLY FOR ROTATIONAL DOF, IF DVDR IS VERY LARGE,
  !!   IT NO LONGER WORKS APPROPRIATELY..
  !!
  !!
  !! - cx:       Chemical structure
  !! - bx:       Chemical structure object
  !! - rad_min:  distance between two molecules hard spheres
  !! - nactmol:  number of active molecules in cx and bx
  !! - MolecAct: the ID of the molecules that are acive in cx and bx (size of array = nactmol)
  !! - gm:       matrix of molecules with atoms bonding/ breaking during a reaction
  !!
  !*************************************************************************
  !
  Subroutine ProjActMolRotTransDVDR2(cx, nactmol, MolecAct, frozlistin)
    implicit none
    type(cxs)          :: cx
    integer            :: i, k, info, i1, n, nn, mm, nactmol, n1
    integer            :: MolecAct(nactmol)
    double precision   :: rr(3*cx%na,3*cx%na), DD(3*cx%na,3*cx%na), lam2(3*cx%na), sig2(3*cx%na,3*cx%na), WORK2(500)
    double precision   :: PPT(cx%na*3,cx%na*3), PP(cx%na*3,cx%na*3), com(3), WORK(500), small
    double precision   :: sig(3*cx%na,3*cx%na),lam(3*cx%na), DD2(3*cx%na,3*cx%na)
    logical            :: frozlist(nactmol)
    logical, optional  :: frozlistin(nactmol)

    if (present(frozlistin)) then
      frozlist = frozlistin
    else
      frozlist = .false.
    endif

    PPT = 0.0d0 ; small = 0.00000010d0
    rr = 0.0d0 ; nn = 0

    do n1 = 1, nactmol
      if (frozlist(n1)) cycle
      n = MolecAct(n1)
      com = 0.0d0
      mm = 0
      do i1 = 1, cx%namol(n)
        i = cx%molid(n, i1)
        com = com + cx%r(1:3, i)
        mm = mm + 1
      enddo
      com = com/dble(mm)
      do i1 = 1, cx%namol(n) ; i = cx%molid(n, i1)
        cx%r(1:3, i) = cx%r(1:3, i) - com
      enddo
      if (cx%namol(n) .gt. 1) then
        nn = nn + 1
        do i1 = 1, cx%namol(n) ; i = cx%molid(n, i1)
         rr(nn, (i-1)*3+2) = +1.0d0*cx%r(3, i)
         rr(nn, (i-1)*3+3) = -1.0d0*cx%r(2, i)
        enddo
        nn = nn + 1
        do i1 = 1, cx%namol(n) ; i = cx%molid(n, i1)
         rr(nn, (i-1)*3+3) = +1.0d0*cx%r(1, i)
         rr(nn, (i-1)*3+1) = -1.0d0*cx%r(3, i)
        enddo
        nn = nn + 1
        do i1 = 1, cx%namol(n) ; i = cx%molid(n, i1)
         rr(nn, (i-1)*3+1) = +1.0d0*cx%r(2, i)
         rr(nn, (i-1)*3+2) = -1.0d0*cx%r(1, i)
        enddo
      endif
      nn = nn + 1
      do i1 = 1, cx%namol(n) ; i = cx%molid(n, i1)
       rr(nn, (i-1)*3+1) = +1.0d0
      enddo
      nn = nn + 1
      do i1 = 1, cx%namol(n) ; i = cx%molid(n, i1)
       rr(nn, (i-1)*3+2) = +1.0d0
      enddo
      nn = nn + 1
      do i1 = 1, cx%namol(n) ; i = cx%molid(n, i1)
       rr(nn, (i-1)*3+3) = +1.0d0
      enddo
      do i1 = 1, cx%namol(n) ; i = cx%molid(n, i1)
        cx%r(1:3, i) = cx%r(1:3, i)+com
      enddo
    enddo
    if (nn .eq. 0) return
    !DD2(1:nn,1:nn) = matmul(transpose(rr(1:nn,:)),(rr(1:nn,:)))
    DD2(1:nn, 1:nn) = matmul((rr(1:nn, :)), transpose(rr(1:nn, :)))
    call dsyev('V', 'U', nn, DD2(1:nn, 1:nn), nn, lam(1:nn), WORK, 500, info)
    !print*, 'LAM = ', lam(1:nn)
    k = 1
    !do while (abs(lam(k)) .le. 1.0d-5)
    ! k = k + 1
    !enddo
    sig = 0.0d0
    do i = k, nn
       if (abs(lam(i)) .lt. 0.0001) lam(i) = lam(i) + 0.0001
       !if (abs(lam(i)) .lt. 0.00000001) lam(i) = lam(i) + 0.00000001
       sig(i, i) = 1.0d0/lam(i)
    enddo
    DD(k:nn, k:nn) = matmul((DD2(k:nn, k:nn)), matmul(sig(k:nn ,k:nn), transpose(DD2(k:nn, k:nn))))
    PP = matmul(transpose(rr(1:nn, :)), matmul(DD(1:nn, 1:nn), rr(1:nn, :)))
    !print*, 'ZEROP? = ', sum(PP-matmul(PP,PP))
    cx%dvdr = reshape(matmul(PP, reshape(cx%dvdr, (/3*cx%na/))), (/3,cx%na/))
    return
  end Subroutine ProjActMolRotTransDVDR2


  !
  !*************************************************************************
  !> ProjOutActMolRotTransDVDR
  !!   projects out the rotational and translational degrees of freedom of the system of active molecules
  !!   from the dcerivative vector dvdr
  !!
  !!
  !! - cx:       Chemical structure
  !! - bx:       Chemical structure object
  !! - rad_min:  distance between two molecules hard spheres
  !! - nactmol:  number of active molecules in cx and bx
  !! - MolecAct: the ID of the molecules that are acive in cx and bx (size of array = nactmol)
  !! - mcx (optional): if another cx is provided, then it uses the molid and namol arrays from that one
  !!
  !*************************************************************************
  !
  Subroutine ProjOutActMolRotTransDVDR(cx, nactmol, MolecAct, mcx)
    implicit none
    type(cxs)          :: cx
    type(cxs), optional  :: mcx
    integer            :: i, k, info, i1, n, nn, mm, nactmol, n1
    integer            :: MolecAct(nactmol)
    double precision   :: rr(3*cx%na,3*cx%na), DD(3*cx%na,3*cx%na), lam2(3*cx%na), sig2(3*cx%na,3*cx%na), WORK2(500)
    double precision   :: PPT(cx%na*3,cx%na*3), PP(cx%na*3,cx%na*3), com(3), WORK(500), small
    double precision   :: sig(3*cx%na,3*cx%na),lam(3*cx%na), DD2(3*cx%na,3*cx%na), unt(3*cx%na,3*cx%na)
    integer, allocatable :: molid(:,:), namol(:), nmol


    if (present(mcx)) then
       allocate(molid(mcx%na, mcx%na), namol(mcx%na))
       molid = mcx%molid
       namol = mcx%namol
       nmol = mcx%nmol
    else
       allocate(molid(cx%na, cx%na), namol(cx%na))
       molid = cx%molid
       namol = cx%namol
       nmol = cx%nmol
    endif
    PPT = 0.0d0 ; small = 0.00000010d0
    rr = 0.0d0 ; nn = 0
    mm = 0
    com = 0.0d0
    !do n1 = 1, nactmol ; n = MolecAct(n1)
    ! print*,n, MolecularFormula(cx,n)
    !enddo
    do n1 = 1, nactmol ; n = MolecAct(n1)
      do i1 = 1, namol(n) ; i = molid(n, i1)
        com = com + cx%r(1:3, i)
        mm = mm + 1
      enddo
    enddo
    com = com/dble(mm)
    do n1 = 1, nactmol ; n = MolecAct(n1)
      do i1 = 1, namol(n) ; i = molid(n, i1)
        cx%r(1:3, i) = cx%r(1:3, i)-com
      enddo
    enddo
    do n1 = 1, nactmol ; n = MolecAct(n1)
      !print*, 'FORMULA = ', MolecularFormula(cx,n)
      do i1 = 1, namol(n) ; i = molid(n, i1)
       !print*, 'ATOM = ', i,cx%atomlabel(i)
       rr(1,(i-1)*3+2) = -1.0d0*cx%r(3,i)
       rr(1,(i-1)*3+3) = +1.0d0*cx%r(2,i)
       rr(2,(i-1)*3+3) = +1.0d0*cx%r(1,i)
       rr(2,(i-1)*3+1) = -1.0d0*cx%r(3,i)
       rr(3,(i-1)*3+1) = +1.0d0*cx%r(2,i)
       rr(3,(i-1)*3+2) = -1.0d0*cx%r(1,i)
       rr(4,(i-1)*3+1) = 1.0d0
       rr(5,(i-1)*3+2) = 1.0d0
       rr(6,(i-1)*3+3) = 1.0d0
      enddo
    enddo
    nn = 6
    do n1 = 1, nactmol ; n = MolecAct(n1)
      do i1 = 1, namol(n) ; i = molid(n, i1)
        cx%r(1:3, i) = cx%r(1:3, i)+com
      enddo
    enddo
    DD2(1:nn, 1:nn) = matmul((rr(1:nn, :)), transpose(rr(1:nn, :)))
    call dsyev('V', 'U', nn, DD2(1:nn, 1:nn), nn, lam(1:nn), WORK, 500, info)
    !print*, 'LAM = ', lam(1:nn)
    !print*, 'INFO = ', info
    k = 1
    do while (lam(k) .le. 0.0d0)
      k = k + 1
    enddo
    !print*, 'K = ', k
    sig = 0.0d0 ;  unt = 0.0d0
    do i = k, nn
      if (lam(i) .eq. 0.0d0) lam(i) = lam(i) + 0.000001
      sig(i, i) = 1.0d0/lam(i)
    enddo
    do i = 1, cx%na*3
      unt(i, i) = 1.0d0
    enddo
    !print*, 'NA = ', cx%na
    DD(k:nn, k:nn) = matmul((DD2(k:nn, k:nn)), matmul(sig(k:nn, k:nn), transpose(DD2(k:nn, k:nn))))
    !print*, 'ZERO UN =? ', norm2(matmul(DD(1:nn,1:nn),matmul((rr(1:nn,:)),transpose(rr(1:nn,:))))-unt)
    PP = (matmul(transpose(rr(k:nn, :)), matmul(DD(k:nn, k:nn), rr(k:nn, :))))
    !print*, 'ZEROP? = ', sum(PP-matmul(PP,PP))
    cx%dvdr = reshape(matmul(unt-PP, reshape(cx%dvdr, (/3*cx%na/))), (/3, cx%na/))

    return
  end Subroutine ProjOutActMolRotTransDVDR


  !
  !*************************************************************************
  !> ProjOutActMolRotTransDVDRPair
  !!
  !! removes the overal translation and rotation DOF of the superset of a pair of cx(1:2) from the DVDR,
  !! for molecules that are "active" (break/form bonds based on graphs in images at endpoints)
  !! it uses the molecule id's of the image ic
  !! note that it does not check weather the molecules are frozen or not...
  !!
  !! - rp: Reaction Path Object
  !! - if:
  !!
  !*************************************************************************
  !
  Subroutine ProjOutActMolRotTransDVDRPair(cx, MolecAct, nactmol, ic)
    implicit none
    type(cxs) :: cx(2)
    integer :: i, j, k, l, n, m, i1, im, info, ai, bi, n1, id, nactmol, infos, nd, j1, ic, nn, mm, na
    integer :: MolecAct(cx(ic)%nmol), molid(cx(ic)%na,cx(ic)%na), namol(cx(ic)%na), nmol
    real(8) :: PP(cx(1)%na*3*2,cx(1)%na*3*2), com(3), WORK(1000)
    real(8) :: DD(cx(1)%na*3*2,cx(1)%na*3*2)
    real(8) :: DD2(cx(1)%na*3*2,cx(1)%na*3*2)
    real(8) :: rr(6,cx(1)%na*3*2), unt(cx(1)%na*3*2,cx(1)%na*3*2)
    real(8) :: lam(cx(1)%na*3*2),sig(cx(1)%na*3*2,cx(1)%na*3*2)
    real(8) :: dvdr(cx(1)%na*3*2)

    na = cx(ic)%na
    molid = cx(ic)%molid
    namol = cx(ic)%namol
    nmol  = cx(ic)%nmol
    PP = 0.0d0
    rr = 0.0d0 ; nn = 0
    mm = 0
    com = 0.0d0
    do im= 1, 2
      do n1 = 1, nactmol ; n = MolecAct(n1)
        do i1 = 1, namol(n) ; i = molid(n, i1)
          com = com + cx(im)%r(1:3, i)
          mm = mm + 1
        enddo
      enddo
    enddo
    com = com/dble(mm)
    do im= 1, 2
      do n1 = 1, nactmol ; n = MolecAct(n1)
        do i1 = 1, namol(n) ; i = molid(n, i1)
          cx(im)%r(1:3, i) = cx(im)%r(1:3, i)-com
        enddo
      enddo
    enddo
    do im= 1, 2
      do n1 = 1, nactmol ; n = MolecAct(n1)
        mm = + 1
        do i1 = 1, namol(n) ; i = molid(n, i1)
          rr(mm, (im-1)*na*3 + (i-1)*3+2) = + 1.0d0*cx(im)%r(3, i)
          rr(mm, (im-1)*na*3 + (i-1)*3+3) = - 1.0d0*cx(im)%r(2, i)
        enddo
        mm = + 2
        do i1 = 1, namol(n) ; i = molid(n, i1)
          rr(mm, (im-1)*na*3 + (i-1)*3+3) = +1.0d0*cx(im)%r(1, i)
          rr(mm, (im-1)*na*3 + (i-1)*3+1) = -1.0d0*cx(im)%r(3, i)
        enddo
        mm = + 3
        do i1 = 1, namol(n) ; i = molid(n, i1)
          rr(mm, (im-1)*na*3 + (i-1)*3+1) = +1.0d0*cx(im)%r(2, i)
          rr(mm, (im-1)*na*3 + (i-1)*3+2) = -1.0d0*cx(im)%r(1, i)
        enddo
        mm = + 4
        do i1 = 1, namol(n) ; i = molid(n, i1)
          rr(mm, (im-1)*na*3 + (i-1)*3+1) = 1.0d0
        enddo
        mm = + 5
        do i1 = 1, namol(n) ; i = molid(n, i1)
          rr(mm, (im-1)*na*3 + (i-1)*3+2) = 1.0d0
        enddo
        mm = + 6
        do i1 = 1, namol(n) ; i = molid(n, i1)
          rr(mm, (im-1)*na*3 + (i-1)*3+3) = 1.0d0
        enddo
      enddo
    enddo
    do im= 1, 2
      do n1 = 1, nactmol ; n = MolecAct(n1)
        do i1 = 1, namol(n) ; i = molid(n, i1)
          cx(im)%r(1:3, i) = cx(im)%r(1:3, i)+com
        enddo
      enddo
    enddo
    nn = mm
    DD2(1:nn, 1:nn) = matmul((rr(1:nn, :)), transpose(rr(1:nn, :)))
    call dsyev('V', 'U', nn, DD2(1:nn, 1:nn), nn, lam(1:nn) ,WORK, 1000, info)
    !print*, 'LAM = ', lam(1:nn)
    k = 1
    !do while (abs(lam(k) .le. 0.0d0)
    ! k = k + 1
    !enddo
    sig = 0.0d0 ;  unt = 0.0d0
    do i = k, nn
      if (lam(i) .eq. 0.0d0) lam(i) = lam(i) + 0.000001
      sig(i, i) = 1.0d0/lam(i)
    enddo
    do i = 1, cx(1)%na*3*2
      unt(i, i) = 1.0d0
    enddo
    DD(k:nn, k:nn) = matmul((DD2(k:nn, k:nn)), matmul(sig(k:nn, k:nn), transpose(DD2(k:nn, k:nn))))
    PP = matmul(transpose(rr(k:nn, :)),matmul(DD(k:nn, k:nn), rr(k:nn, :)))
    print *, 'ZEROP? = ', sum(PP-matmul(PP, PP))
    do im = 1, 2
      dvdr((im-1)*na*3+1:(im-1)*na*3+na*3) = reshape(cx(im)%dvdr, (/na*3/))
    enddo
    dvdr = matmul(unt-PP, dvdr)
    do im = 1, 2
      cx(im)%dvdr = reshape(dvdr((im-1)*na*3+1:(im-1)*na*3+na*3), (/3,na/))
    enddo
    return
  end subroutine ProjOutActMolRotTransDVDRPair


  !
  !*************************************************************************
  !> ProjOutActMolRotTransDVDR2
  !!
  !! removes the overal translation and rotation DOF of the entire rp coordinates from the DVDR,
  !! for molecules that are "active" (break/form bonds based on graphs in images at endpoints)
  !! it uses the molecule id's of the image ic
  !! note that it does not check weather the molecules are frozen or not...
  !!
  !! - rp: Reaction Path Object
  !! - if:
  !!
  !*************************************************************************
  !
  Subroutine ProjOutActMolRotTransDVDR2(cx,MolecAct, nactmol, ic)
    implicit none
    type(cxs) :: cx(2)
    integer :: i, j, k, l, n, m, i1, im, info, ai, bi, n1, id, nactmol, infos, nd, j1, ic, nn, mm, na
    integer :: MolecAct(nactmol), molid(cx(ic)%na,cx(ic)%na), namol(cx(ic)%na), nmol
    real(8) :: PP(cx(1)%na*3*2,cx(1)%na*3*2), com(3), WORK(1000)
    real(8) :: DD(cx(1)%na*3*2,cx(1)%na*3*2)
    real(8) :: DD2(cx(1)%na*3*2,cx(1)%na*3*2)
    real(8) :: rr(6,cx(1)%na*3*2), unt(cx(1)%na*3*2,cx(1)%na*3*2)
    real(8) :: lam(cx(1)%na*3*2),sig(cx(1)%na*3*2,cx(1)%na*3*2)
    real(8) :: dvdr(cx(1)%na*3*2)

    na = cx(ic)%na
    molid = cx(ic)%molid
    namol = cx(ic)%namol
    nmol  = cx(ic)%nmol
    PP = 0.0d0
    rr = 0.0d0 ; nn = 0
    mm = 0
    com = 0.0d0
    do im= 1, 2
      do n1 = 1, nactmol ; n = MolecAct(n1)
        do i1 = 1, namol(n) ; i = molid(n, i1)
          com = com + cx(im)%r(1:3, i)
          mm = mm + 1
        enddo
      enddo
    enddo
    com = com/dble(mm)
    do im= 1, 2
      do n1 = 1, nactmol ; n = MolecAct(n1)
        do i1 = 1, namol(n) ; i = molid(n, i1)
          cx(im)%r(1:3, i) = cx(im)%r(1:3, i)-com
        enddo
      enddo
    enddo
    do im= 1, 2
      do n1 = 1, nactmol ; n = MolecAct(n1)
        do i1 = 1, namol(n) ; i = molid(n,i1)
          rr(1, (im-1)*na*3 + (i-1)*3+2) = -1.0d0*cx(im)%r(3, i)
          rr(1, (im-1)*na*3 + (i-1)*3+3) = +1.0d0*cx(im)%r(2, i)
          rr(2, (im-1)*na*3 + (i-1)*3+3) = +1.0d0*cx(im)%r(1, i)
          rr(2, (im-1)*na*3 + (i-1)*3+1) = -1.0d0*cx(im)%r(3, i)
          rr(3, (im-1)*na*3 + (i-1)*3+1) = +1.0d0*cx(im)%r(2, i)
          rr(3, (im-1)*na*3 + (i-1)*3+2) = -1.0d0*cx(im)%r(1, i)
          rr(4, (im-1)*na*3 + (i-1)*3+1) = 1.0d0
          rr(5, (im-1)*na*3 + (i-1)*3+2) = 1.0d0
          rr(6, (im-1)*na*3 + (i-1)*3+3) = 1.0d0
        enddo
      enddo
    enddo
    do im= 1, 2
      do n1 = 1, nactmol ; n = MolecAct(n1)
        do i1 = 1, namol(n) ; i = molid(n, i1)
          cx(im)%r(1:3, i) = cx(im)%r(1:3, i)+com
        enddo
      enddo
    enddo
    mm = 6
    nn = mm
    DD2(1:nn, 1:nn) = matmul((rr(1:nn, :)), transpose(rr(1:nn, :)))
    call dsyev('V', 'U', nn, DD2(1:nn, 1:nn), nn, lam(1:nn), WORK, 1000, info)
    !print*, 'LAM = ', lam(1:nn)
    k = 1
    !do while (lam(k) .le. 0.0d0)
    ! k = k + 1
    !enddo
    sig = 0.0d0 ;  unt = 0.0d0
    do i = k, nn
      if (lam(i) .eq. 0.0d0) lam(i) = lam(i) + 0.000001
      sig(i, i) = 1.0d0/lam(i)
    enddo
    do i = 1, na*3*2
      unt(i, i) = 1.0d0
    enddo
    !DD(k:nn,k:nn) = matmul(transpose(DD2(k:nn,k:nn)),matmul(sig(k:nn,k:nn),(DD2(k:nn,k:nn))))
    !PP = matmul(transpose(rr(1:nn,:)),matmul(transpose(DD(1:nn,1:nn)),rr(1:nn,:)))
    DD(k:nn, k:nn) = matmul((DD2(k:nn, k:nn)), matmul(sig(k:nn, k:nn), transpose(DD2(k:nn, k:nn))))
    !print*, 'ZERO1 = ', norm2(unt(k:nn,k:nn)-matmul(DD(k:nn,k:nn),matmul((rr(1:nn,:)),transpose(rr(1:nn,:)))))
    PP = matmul(transpose(rr(1:nn, :)), matmul(DD(1:nn, 1:nn), rr(1:nn, :)))
    !print*, 'ZEROP? = ', norm2(PP-matmul(PP,PP)), norm2(PP)
    do im = 1, 2
      dvdr((im-1)*na*3+1:(im-1)*na*3+na*3) = reshape(cx(im)%dvdr, (/na*3/))
    enddo
    !l = 0
    !do i = 1, na
    ! do j = 1, 3
    !   l = l + 1
    !   print*, 'DIFF = ', cx(2)%dvdr(j,i) - dvdr(na*3+ l)
    ! enddo
    !enddo
    !print*, 'ZEROP2? = ', norm2((unt-PP)-matmul(unt-PP,unt-PP)), norm2(unt-PP)
    !
    !print*, 'LEN 1 = ', norm2(dvdr)
    !print*, 'DVDR 1 = ', dvdr(1:na*3)
    dvdr = matmul(unt-transpose(PP), dvdr)
    !print*, 'LEN 2 = ', norm2(dvdr(1:na*3))
    !print*, 'DVDR 1 = ', dvdr(1:na*3)
    do im = 1, 2
      cx(im)%dvdr = reshape(dvdr((im-1)*na*3+1:(im-1)*na*3+na*3), (/3 ,na/))
    enddo
    return
  end subroutine


  !
  !*************************************************************************
  !> ProjMolRotTransDVDR2
  !!   projects the rotational and translational degrees of freedom of each molecule from the dcerivative vector dvdr
  !!
  !! - cx:       Chemical structure
  !! - bx:       Chemical structure object
  !! - rad_min:  distance between two molecules hard spheres
  !! - nactmol:  number of active molecules in cx and bx
  !! - MolecAct: the ID of the molecules that are acive in cx and bx (size of array = nactmol)
  !! - gm:       matrix of molecules with atoms bonding/ breaking during a reaction
  !!
  !*************************************************************************
  !
  Subroutine ProjMolRotTransDVDR2(cx)
    implicit none
    type(cxs)          :: cx
    integer            :: i, k, info, i1, n, nn, mm
    double precision   :: rr(3*cx%na,3*cx%na), DD(3*cx%na,3*cx%na), lam2(3*cx%na), sig2(3*cx%na,3*cx%na), WORK2(500), rr2(1,3*cx%na)
    double precision   :: PPT(cx%na*3,cx%na*3), PP(cx%na*3,cx%na*3), com(3), WORK(500), small
    double precision   :: sig(3*cx%na,3*cx%na),lam(3*cx%na), DD2(3*cx%na,3*cx%na)

    PPT = 0.0d0 ; small = 0.00000010d0
    rr = 0.0d0 ; nn = 0 ; mm = 0
    do n = 1, cx%nmol
      com = 0.0d0
      do i1 = 1, cx%namol(n)
        i = cx%molid(n,i1)
        com = com + cx%r(1:3,i)
        mm = mm + 1
      enddo
      com = com/dble(mm)
      do i1 = 1, cx%namol(n)
        i = cx%molid(n,i1)
        cx%r(1:3,i) = cx%r(1:3,i)-com
      enddo
      if (cx%namol(n) .gt. 1) then
        nn = nn + 1
        do i1 = 1, cx%namol(n)
          i = cx%molid(n,i1)
          rr(nn,(i-1)*3+2) = +1.0d0*cx%r(3,i)
          rr(nn,(i-1)*3+3) = -1.0d0*cx%r(2,i)
        enddo
        nn = nn + 1
        do i1 = 1, cx%namol(n)
          i = cx%molid(n,i1)
          rr(nn,(i-1)*3+3) = +1.0d0*cx%r(1,i)
          rr(nn,(i-1)*3+1) = -1.0d0*cx%r(3,i)
        enddo
        nn = nn + 1
        do i1 = 1, cx%namol(n)
          i = cx%molid(n,i1)
          rr(nn,(i-1)*3+1) = +1.0d0*cx%r(2,i)
          rr(nn,(i-1)*3+2) = -1.0d0*cx%r(1,i)
        enddo
      endif
      nn = nn + 1
      do i1 = 1, cx%namol(n)
        i = cx%molid(n,i1)
        rr(nn,(i-1)*3+1) = +1.0d0
      enddo
      nn = nn + 1
      do i1 = 1, cx%namol(n)
        i = cx%molid(n,i1)
        rr(nn,(i-1)*3+2) = +1.0d0
      enddo
      nn = nn + 1
      do i1 = 1, cx%namol(n)
        i = cx%molid(n,i1)
        rr(nn,(i-1)*3+3) = +1.0d0
      enddo

      do i1 = 1, cx%namol(n)
        i = cx%molid(n,i1)
        cx%r(1:3,i) = cx%r(1:3,i)+com
      enddo
    enddo

    DD2(1:nn, 1:nn) = matmul((rr(1:nn, :)),transpose(rr(1:nn, :)))
    call dsyev('V', 'U', nn, DD2(1:nn, 1:nn), nn, lam(1:nn), WORK, 500, info)
    !print*, 'LAM = ', lam(1:nn)
    k = 1
    !do while (lam(k) .le. 0.0d0)
    ! k = k + 1
    !enddo
    sig = 0.0d0
    do i = k, nn
      if (lam(i) .eq. 0.0d0) lam(i) = lam(i) + 0.000001
      sig(i, i) = 1.0d0/lam(i)
    enddo
    DD(k:nn, k:nn) = matmul((DD2(k:nn, k:nn)), matmul(sig(k:nn, k:nn), transpose(DD2(k:nn, k:nn))))
    PP = matmul(transpose(rr(1:nn, :)), matmul(DD(1:nn, 1:nn), rr(1:nn, :)))
    !print*, 'ZEROP? = ', sum(PP-matmul(PP,PP))
    cx%dvdr = reshape(matmul(PP, reshape(cx%dvdr, (/3*cx%na/))), (/3,cx%na/))
    return
  end Subroutine ProjMolRotTransDVDR2

  !
  !*************************************************************************
  !> ProjMolRotTransDVDR
  !!   projects the rotational and translational degrees of freedom of each molecule from the dcerivative vector dvdr
  !!   not sure if this one works...
  !!
  !! - cx:       Chemical structure
  !! - bx:       Chemical structure object
  !! - rad_min:  distance between two molecules hard spheres
  !! - nactmol:  number of active molecules in cx and bx
  !! - MolecAct: the ID of the molecules that are acive in cx and bx (size of array = nactmol)
  !! - gm:       matrix of molecules with atoms bonding/ breaking during a reaction
  !!
  !*************************************************************************
  !
  Subroutine ProjMolRotTransDVDR(cx)
    implicit none
    type(cxs)          :: cx
    integer            :: i, k, info, i1, n, rr3(3*cx%na,3*cx%na)
    double precision   :: rr(6,3*cx%na), DD(6,6), lam2(3*cx%na), sig2(3*cx%na,3*cx%na), WORK2(500), rr2(1,3*cx%na)
    double precision   :: PPT(cx%na*3,cx%na*3), PP(cx%na*3,cx%na*3), com(3), WORK(50), sig(6,6),lam(6), DD2(6,6)

    PPT = 0.0d0
    rr = 0.0d0
    do n = 1, cx%nmol
      if (cx%namol(n) .eq. 1) cycle
      com = 0.0d0
      do i1 = 1, cx%namol(n) ; i = cx%molid(n,i1)
        com = com + cx%r(1:3,i)
      enddo
      com = com/dble(cx%namol(n))
      do i1 = 1, cx%namol(n) ; i = cx%molid(n,i1)
        cx%r(1:3,i) = cx%r(1:3,i)-com
      enddo
      !rr = 0.0d0
      dO i1 = 1, cx%namol(n) ; i = cx%molid(n,i1)

       rr(1,(i-1)*3+2) = +1.0d0*cx%r(3,i)
       rr(1,(i-1)*3+3) = -1.0d0*cx%r(2,i)

       rr(2,(i-1)*3+3) = +1.0d0*cx%r(1,i)
       rr(2,(i-1)*3+1) = -1.0d0*cx%r(3,i)

       rr(3,(i-1)*3+1) = +1.0d0*cx%r(2,i)
       rr(3,(i-1)*3+2) = -1.0d0*cx%r(1,i)

       rr(4,(i-1)*3+1) = +1.0d0
       rr(5,(i-1)*3+2) = +1.0d0
       rr(6,(i-1)*3+3) = +1.0d0
      enddo
      do i1 = 1, cx%namol(n) ; i = cx%molid(n,i1)
        cx%r(1:3,i) = cx%r(1:3,i)+com
      enddo
      !DD2 = matmul((rr),transpose(rr))
      !Call dsyev('V','U',6,DD2(1:6,1:6),6,lam(1:6),WORK,50,info)
      !print*, 'LAM = ', lam
      !sig = 0.0d0
      !do i = 1, 6
      !   sig(i,i) = 1.0d0/lam(i)
      !enddo
      !!if (cx%namol(n) .eq. 2) then
      ! !DD(2:6,2:6) = matmul((DD2(2:6,2:6)),matmul(sig(2:6,2:6),transpose(DD2(2:6,2:6))))
      ! !PP = matmul(transpose(rr(2:6,:)),matmul(DD(2:6,2:6),rr(2:6,:)))
      !!else
      !DD = matmul((DD2),matmul(sig,transpose(DD2)))
      !PP = matmul(transpose(rr),matmul(DD,rr))
      PP = matmul(transpose(rr),rr)
      !!endif
      !print*, 'NAMOL =- ', cx%namol(n)
      !print*, 'ZEROP? = ', sum(PP-matmul(PP,PP))
      PPT = PPT + PP
    enddo
    ! also permit the rotation of the overall system...
    com = 0.0d0
    do i = 1, cx%na
      com = com + cx%r(1:3,i)
    enddo
    com = com/dble(cx%na)
    do i = 1, cx%na
      cx%r(1:3,i) = cx%r(1:3,i)-com
    enddo
    rr = 0.0d0
    do i = 1, cx%na
     rr(1,(i-1)*3+2) = rr(1,(i-1)*3+2) +1.0d0*cx%r(3,i)
     rr(1,(i-1)*3+3) = rr(1,(i-1)*3+3) -1.0d0*cx%r(2,i)

     rr(2,(i-1)*3+3) = rr(2,(i-1)*3+3) +1.0d0*cx%r(1,i)
     rr(2,(i-1)*3+1) = rr(2,(i-1)*3+1) -1.0d0*cx%r(3,i)

     rr(3,(i-1)*3+1) = rr(3,(i-1)*3+1) +1.0d0*cx%r(2,i)
     rr(3,(i-1)*3+2) = rr(3,(i-1)*3+2) -1.0d0*cx%r(1,i)
     rr(4,(i-1)*3+1) = rr(4,(i-1)*3+1) +1.0d0
     rr(5,(i-1)*3+2) = rr(5,(i-1)*3+2) +1.0d0
     rr(6,(i-1)*3+3) = rr(6,(i-1)*3+3) +1.0d0
    enddo
    do i = 1, cx%na
      cx%r(1:3,i) = cx%r(1:3,i)+com
    enddo
    !DD2(1:6,1:6) = matmul((rr(1:6,:)),transpose(rr(1:6,:)))
    !Call dsyev('V','U',6,DD2(1:6,1:6),6,lam(1:6),WORK,50,info)
    !print*, 'LAM = ', lam
    !sig = 0.0d0
    !do i = 1, 6
    !   sig(i,i) = 1.0d0/lam(i)
    !enddo
    !DD(1:6,1:6) = matmul((DD2(1:6,1:6)),matmul(sig(1:6,1:6),transpose(DD2(1:6,1:6))))
    !PP = matmul(transpose(rr(1:6,:)),matmul(DD(1:6,1:6),rr(1:6,:)))
    !DD2(1:3,1:3) = matmul((rr(1:3,:)),transpose(rr(1:3,:)))
    !Call dsyev('V','U',3,DD2(1:3,1:3),3,lam(1:3),WORK,50,info)
    !print*, 'LAM = ', lam
    !sig = 0.0d0
    !do i = 1, 3
    !   sig(i,i) = 1.0d0/lam(i)
    !enddo
    !DD(1:3,1:3) = matmul((DD2(1:3,1:3)),matmul(sig(1:3,1:3),transpose(DD2(1:3,1:3))))
    !PP = matmul(transpose(rr(1:3,:)),matmul(DD(1:3,1:3),rr(1:3,:)))
    PP = matmul(transpose(rr),rr)
    PPT = PPT + PP
    PP = matmul(transpose(PPT),(PPT))
    rr2 = 1.0d0
    rr2 = transpose(matmul(PP,transpose(rr2)))
    PP = matmul(rr2,transpose(rr2))
    !PP = matmul(PPT,transpose(PPT))
    n = cx%na*3
    Call dsyev('V','U',n,PP,n,lam2,WORK2,500,info)
    sig2 = 0.0d0 ; n = 0
    do i = 1, cx%na*3
       if (abs(lam2(i)) .gt. 0.001) then
        n = n + 1
        sig2(:,n) = PP(:,i)
       endif
    enddo
    PPT = matmul(sig2(:,1:n),transpose(sig2(:,1:n)))
    !print*, 'ZEROP? = ', sum(PP-matmul(PP,PP))
    !print*, 'TOTAL ZEROP? = ', sum(PPT-matmul(PPT,PPT))
    cx%dvdr = reshape(matmul(PP,reshape(cx%dvdr,(/3*cx%na/))),(/3,cx%na/))
    return
  end Subroutine ProjMolRotTransDVDR


  !
  !************************************************************************
  !> MolIdFromAtom
  !!
  !! Tells you which molecule the atom belogs to. Returns the ID of that molecule
  !!
  !! - cx: chemical structure
  !! - atid: atom index
  !!
  !************************************************************************
  !
  function MolIdOfAtom(atid, cx)
    implicit none
    integer   :: j, k, i, atid
    type(cxs) :: cx
    integer   :: MolIdOfAtom

    MolIdOfAtom = 0
    do i = 1, cx%nmol
      do j = 1, cx%namol(i)
        if (cx%molid(i, j) == atid) then
          MolIdOfAtom = i
          return
        endif
      enddo
    enddo
    if (MolIdOfAtom .eq. 0) then
      print*, 'MolIdOfAtom found that atom: ', atid
      print*, 'does not belong to a "molecule"\rpath.f90'
      print*, 'MOLID = '
      do i = 1, cx%nmol
        do j = 1, cx%namol(i)
          print*, i, j, cx%molid(i, j)
        enddo
      enddo
      stop
    endif
    return
  end function MolIdOfAtom


  !
  !*************************************************************************
  !
  !> MoleculeMomentOfInertiaTensor
  !!
  !! Returns the Moment of inertia tensor matrix of a molecule mi in chemical structure cx
  !!
  !! mi - the index of the molecule
  !! mit - returns moment of inertia tensor
  !!
  !*************************************************************************
  !
  Function MolMomentOfInertiaTensor(cx, mi) result(mit)
    implicit none
    type(cxs) :: cx
    integer :: l, mi, k, m, i, j
    double precision :: mit(3, 3), com(3), mmas

    com = 0.0d0 ; mmas = 0.0d0
    do i = 1, cx%namol(mi)
      com(1) = com(1) + cx%mass(i) * cx%r(1, i)
      com(2) = com(2) + cx%mass(i) * cx%r(2, i)
      com(3) = com(3) + cx%mass(i) * cx%r(3, i)
      mmas = mmas + cx%mass(i)
    enddo

    com = com / mmas
    mit = 0.0d0

    do i = 1, cx%namol(mi)
      mit(1, 1) = mit(1, 1) + cx%mass(i) * ((cx%r(2, i) - com(2))**2 + (cx%r(3, i) - com(3))**2)
      mit(2, 2) = mit(2, 2) + cx%mass(i) * ((cx%r(1, i) - com(1))**2 + (cx%r(3, i) - com(3))**2)
      mit(3, 3) = mit(3, 3) + cx%mass(i) * ((cx%r(1, i) - com(1))**2 + (cx%r(2, i) - com(2))**2)
      mit(1, 2) = mit(1, 2) - cx%mass(i) * ((cx%r(1, i) - com(1)) * (cx%r(2, i) - com(2)))
      mit(1, 3) = mit(1, 3) - cx%mass(i) * ((cx%r(1, i) - com(1)) * (cx%r(3, i) - com(3)))
      mit(2, 3) = mit(2, 3) - cx%mass(i) * ((cx%r(2, i) - com(2)) * (cx%r(3, i) - com(3)))
    enddo
    mit(2, 1) = mit(1, 2)
    mit(3, 1) = mit(1, 3)
    mit(3, 2) = mit(2, 3)
    return
  end function MolMomentOfInertiaTensor


  !
  !*************************************************************************
  !> MolecularCOM
  !!
  !! returns the centre of mass of  molecule m from reaction strucuture cx
  !!
  !! - cx: Chemical structure object.
  !! - m: molecule in cxs
  !!
  !*************************************************************************
  !
  Function MolecularCOM(cx, m)
    implicit none
    type(cxs) :: cx
    integer :: i, j, im, k, m
    real(8) :: MolecularCOM(3), MM

    MolecularCOM = 0.0d0 ; MM = 0.0d0
    do i = 1, cx%namol(m)
      MolecularCOM = MolecularCOM + cx%r(1:3, cx%molid(m, i)) * cx%mass(i)
      MM = MM + cx%mass(i)
    enddo
    if (cx%namol(m) .eq. 0 .or. cx%namol(m) .ne. cx%namol(m) ) &
      print*, 'HMM...  NAMOL ZERO ? = ', cx%namol(m), m
    MolecularCOM = MolecularCOM/MM
    return
  end Function MolecularCOM


  !
  !*************************************************************************
  !> AccumulateDerivatives
  !!
  !! Updates dvdr with additional forces.
  !!
  !! - cx: Chemical structure object.
  !! - t1: Force term to be applied.
  !! - i, j: Atom indexes.
  !!
  !*************************************************************************
  !
  Subroutine AccumulateDerivatives(cx, t1 , i, j)
    implicit none
    type(cxs) :: cx
    real(8) :: dx, dy, dz, rr, rsq, t1, onr
    integer :: i, j

    dx = cx%r(1, i) - cx%r(1, j)
    dy = cx%r(2, i) - cx%r(2, j)
    dz = cx%r(3, i) - cx%r(3, j)
    rsq = dx*dx + dy*dy + dz*dz
    rr = sqrt(rsq)
    if (rr .gt. epsil) then
      onr = 1.d0 / rr
    else
      onr = 0.0d0
    endif

    cx%dvdr(1, i) = cx%dvdr(1, i) + t1 * dx * onr
    cx%dvdr(2, i) = cx%dvdr(2, i) + t1 * dy * onr
    cx%dvdr(3, i) = cx%dvdr(3, i) + t1 * dz * onr

    cx%dvdr(1, j) = cx%dvdr(1, j) - t1 * dx * onr
    cx%dvdr(2, j) = cx%dvdr(2, j) - t1 * dy * onr
    cx%dvdr(3, j) = cx%dvdr(3, j) - t1 * dz * onr

    return
  end Subroutine


  !
  !*************************************************************************
  !> Get3rings
  !!
  !! Counts the number of 3-rings in the cx graph.
  !!
  !! - cx: Chemical structure object.
  !! - nrings: Number of 3-rings identified.
  !!
  !*************************************************************************
  !
  Subroutine Get3Rings(cx, nrings)
    implicit none
    type(cxs) :: cx
    integer :: nrings, i, j, k, na

!!$    do i = 1, cx%na-1
!!$       do j = i, cx%na
!!$          print*,'GRAPH NOW: ',cx%graph(i,j)
!!$       enddo
!!$    enddo

    na = cx%na
    nrings = 0
    do i = 1, na-2
      do j = i+1,na-1
        if (cx%graph(i, j) == 1) then
          do k = j+1, na
            if (k /= i .and. k /=j .and. cx%graph(i,k) == 1 .and.cx%graph(j,k) == 1) then
              nrings = nrings + 1
            endif
          enddo
        endif
      enddo
    enddo
    return
  end Subroutine Get3Rings


  !
  !*************************************************************************
  !> Get4rings
  !!
  !! Counts the number of 4-rings in the cx graph.
  !!
  !! - cx: Chemical structure object.
  !! - nrings: Number of 4-rings identified.
  !!
  !*************************************************************************
  !
  Subroutine Get4Rings(cx, nrings)
    implicit none
    type(cxs) :: cx
    integer :: nrings, i, j, k, l, na

    na = cx%na
    nrings = 0
    do i = 1, na-1
      do j = i+1, na
        if (cx%graph(i, j) == 1) then
          do k = 1, na

            if (k /= i .and. k /=j .and. cx%graph(i, k) == 1.and.cx%graph(j, k)==0) then
              do l = 1, na
                if (l/=i .and. l/=j .and. l /= k) then
                  if (cx%graph(k, l) == 1 .and. cx%graph(l, i) == 1) then
                    nrings = nrings + 1
                  endif
                endif
              enddo
            endif

            if (k /= i .and. k /=j .and. cx%graph(j, k) == 1.and.cx%graph(i, k)==0) then
              do l = 1, na
                if (l/=i .and. l/=j .and. l /= k) then
                  if (cx%graph(k, l) == 1 .and. cx%graph(l, i) == 1) then
                    nrings = nrings + 1
                  endif
                endif
              enddo
            endif

          enddo
        endif
      enddo
    enddo
    return
  end Subroutine Get4Rings


  !
  !************************************************************************
  !> PrintCXSToFile
  !!
  !! Outputs cxs coordinates to an xyz file.
  !!
  !! - cx: The CXS object.
  !! - filename: the output xyz file.
  !! - value: A real value which can be output as the comment line.
  !!          For example, this could be the energy or fitness.
  !!
  !************************************************************************
  !
  Subroutine PrintCXSToFile(cx, filename, value)
    type(cxs) :: cx
    character (len=*) :: filename
    real*8 :: x, y, z, value
    integer :: i, j

    open(13, file=filename, status='unknown')
    write(13, '(i5)') cx%na
    write(13, *) value
    do j = 1, cx%na
      x = cx%r(1, j) * bohr_to_ang
      y = cx%r(2, j) * bohr_to_ang
      z = cx%r(3, j) * bohr_to_ang
      write(13,'(a2, 2x, 3(f15.11, 2x))') cx%atomlabel(j), x, y, z
    enddo
    close(13)

    return
  end Subroutine PrintCXSToFile


  !
  !*************************************************************************
  !
  !> CreateMolecularCX
  !!
  !! Creates a CXS mcx for a single molecule, indexed m in CX
  !!
  !! - cx: the parent chemical structure from which mcx is generated
  !! - mcx: the CXS of a the molecule m in cx
  !! - m: the molecule's index in cx, which we wish to make mcx from..
  !!
  !*************************************************************************
  !
  Subroutine CreateMolecularCX(cx, mcx, m)
    implicit none
    type(cxs) :: cx,mcx
    integer :: i, j, k, l, m
    if ( allocated(mcx%r )        ) deallocate(mcx%r )
    if ( allocated(mcx%p )        ) deallocate(mcx%p )
    if ( allocated(mcx%dvdr )     ) deallocate(mcx%dvdr )
    if ( allocated(mcx%force )    ) deallocate(mcx%force )
    if ( allocated(mcx%mass )     ) deallocate(mcx%mass )
    if ( allocated(mcx%atomlabel )) deallocate(mcx%atomlabel )
    if ( allocated(mcx%fixeddof ) ) deallocate(mcx%fixeddof )
    if ( allocated(mcx%fixedatom )) deallocate(mcx%fixedatom )
    if ( allocated(mcx%graph )    ) deallocate(mcx%graph )
    if ( allocated(mcx%molid )    ) deallocate(mcx%molid )
    if ( allocated(mcx%namol )    ) deallocate(mcx%namol )
    if ( allocated(mcx%molen )    ) deallocate(mcx%molen )
    if ( allocated(mcx%molcharge )) deallocate(mcx%molcharge )
    if ( allocated(mcx%molspin   )) deallocate(mcx%molspin )

    allocate(mcx%r(1:3, cx%namol(m)))
    allocate(mcx%p(1:3, cx%namol(m)))
    allocate(mcx%dvdr(1:3, cx%namol(m)))
    allocate(mcx%force(1:3, cx%namol(m)))
    allocate(mcx%mass(cx%namol(m)))
    allocate(mcx%atomlabel(cx%namol(m)))
    allocate(mcx%fixeddof(3*cx%namol(m)))
    allocate(mcx%fixedatom(cx%namol(m)))
    allocate(mcx%graph(cx%namol(m), cx%namol(m)))
    ! in case the molecule fragments, make namol array > 1
    allocate(mcx%molid(20, cx%namol(m)))
    allocate(mcx%namol(20))
    allocate(mcx%molen(20))
    allocate(mcx%molcharge(nmolmax))
    allocate(mcx%molspin(nmolmax))
    mcx%fixeddof = .false.
    mcx%fixedatom = .false.
    mcx%ndofconstr = 0
    mcx%natomconstr = 0
    mcx%vcalc = 0.0d0
    mcx%vcon = 0.0d0
    mcx%nmol = 1
    mcx%na = cx%namol(m)
    mcx%namol(1) = cx%namol(m)
    if (allocated(cx%molen)) then
      mcx%molen(1) = cx%molen(m)
      mcx%vcalc = cx%molen(m)
    endif
    if (allocated(cx%molcharge)) mcx%molcharge(1) = cx%molcharge(m)
    if (allocated(cx%molspin)) mcx%molspin(1) = cx%molspin(m)
    do j = 1, cx%namol(m)  ! for atoms in molecule
      mcx%molid(1, j) = j
      mcx%r(1:3, j) = cx%r(1:3, cx%molid(m, j))
      mcx%dvdr(1:3, j) = cx%dvdr(1:3, cx%molid(m, j))
      do l = 1, 3
       if (cx%fixeddof((cx%molid(m, j)-1)*3+l)) then
          mcx%ndofconstr = mcx%ndofconstr + 1
          mcx%fixeddof((j-1)*3+l) = .true.
       endif
      enddo
      if (cx%fixedatom(cx%molid(m, j))) then
        mcx%fixedatom(j) = .true.
        mcx%natomconstr = mcx%natomconstr + 1
      endif
      mcx%atomlabel(j) = cx%atomlabel(cx%molid(m, j))
      mcx%mass(j) = MASS(LabelToNumber(mcx%atomlabel(j)))
      do l = 1, cx%namol(m)
        mcx%graph(j, l) = cx%graph(cx%molid(m, j),cx%molid(m, l))
      enddo
    enddo
    mcx%method = cx%method
    return
  end Subroutine CreateMolecularCX


  subroutine CreateSubStructureCX(cx, mcx, id, na)
    implicit none
    integer :: i, j, k, l, m, na, j1, l1
    integer :: id(na)
    type(cxs) :: cx,mcx
    if ( allocated(mcx%r )        ) deallocate(mcx%r )
    if ( allocated(mcx%p )        ) deallocate(mcx%p )
    if ( allocated(mcx%dvdr )     ) deallocate(mcx%dvdr )
    if ( allocated(mcx%force )    ) deallocate(mcx%force )
    if ( allocated(mcx%mass )     ) deallocate(mcx%mass )
    if ( allocated(mcx%atomlabel )) deallocate(mcx%atomlabel )
    if ( allocated(mcx%fixeddof ) ) deallocate(mcx%fixeddof )
    if ( allocated(mcx%fixedatom )) deallocate(mcx%fixedatom )
    if ( allocated(mcx%graph )    ) deallocate(mcx%graph )
    if ( allocated(mcx%molid )    ) deallocate(mcx%molid )
    if ( allocated(mcx%namol )    ) deallocate(mcx%namol )
    if ( allocated(mcx%molen )    ) deallocate(mcx%molen )
    if ( allocated(mcx%molcharge )) deallocate(mcx%molcharge )
    if ( allocated(mcx%molspin   )) deallocate(mcx%molspin )

    allocate(mcx%r(1:3, na))
    allocate(mcx%p(1:3, na))
    allocate(mcx%dvdr(1:3, na))
    allocate(mcx%force(1:3, na))
    allocate(mcx%mass(na))
    allocate(mcx%atomlabel(na))
    allocate(mcx%fixeddof(3*na))
    allocate(mcx%fixedatom(na))
    allocate(mcx%graph(na, na))
    allocate(mcx%molcharge(nmolmax))
    allocate(mcx%molspin(nmolmax))
    mcx%fixeddof = .false.
    mcx%fixedatom = .false.
    mcx%ndofconstr = 0
    mcx%natomconstr = 0
    mcx%vcalc = 0.0d0
    mcx%vcon = 0.0d0
    mcx%na = na
    mcx%method = cx%method
    do j1 = 1, na  ! for atoms in molecule
      j = id(j1)
      mcx%r(1:3, j1) = cx%r(1:3, j)
!      print*, 'j, = ', j, 'id= ', id(j1)
      do l = 1, 3
        if (cx%fixeddof((j-1)*3+l)) then
          mcx%ndofconstr = mcx%ndofconstr + 1
          mcx%fixeddof((j1-1)*3+l) = .true.
        endif
      enddo
      if (cx%fixedatom(j)) then
        mcx%fixedatom(j1) = .true.
        mcx%natomconstr = mcx%natomconstr + 1
      endif
      mcx%atomlabel(j1) = cx%atomlabel(j)
      mcx%mass(j1) = MASS(LabelToNumber(mcx%atomlabel(j1)))
      do l1 = 1, na ; l = id(l1)
        mcx%graph(j1, l1) = cx%graph(j, l)
      enddo
    enddo
    return
  end subroutine


  !
  !*************************************************************************
  !
  !> MolecularGraph
  !!
  !! Returns graph of molecule in cx
  !!
  !! - cx: the parent chemical structure from which mcx is generated
  !! - m: the molecule's index in cx, which we wish to make the graph from..
  !!
  !*************************************************************************
  !
  function MolecularGraph(cx, m) result(G)
    implicit none
    integer              :: i, j, k, l
    integer, intent(in)  :: m
    type(cxs)      :: cx
    integer        :: G(cx%namol(m), cx%namol(m))

    do j = 1, cx%namol(m)  ! for atoms in molecule
      do l = 1, cx%namol(m)
        G(j, l) = cx%graph(cx%molid(m, j),cx%molid(m, l))
      enddo
    enddo
    return
  end function


  !
  !*************************************************************************
  !
  !> MolecularAtomicLabel
  !!
  !! Returns the atomic labels for the molecle
  !!
  !! - cx: the parent chemical structure from which the atomic labels is generated
  !! - m: the molecule's index in cx, which we wish to make the atomic labels from..
  !!
  !*************************************************************************
  !
  function MolecularALabel(cx, m) result(lab)
    implicit none
    integer              :: i, j, k, l
    integer, intent(in)  :: m
    type(cxs)            :: cx
    character(len=2)     :: lab(cx%namol(m))

    do j = 1, cx%namol(m)  ! for atoms in molecule
      lab(j) = cx%atomlabel(cx%molid(m, j))
    enddo
    return
  end function


  !
  !*************************************************************************
  !> RemoveHydrogens
  !!
  !! Removes all hydrogen atoms from a cxs object, and returns a
  !! new cxs object with only heteroatoms.
  !
  !! - cxin: Input chemical structure object.
  !! - cxout: Output chemical structure object, with no hydrogens.
  !!
  !*************************************************************************
  !
  Subroutine RemoveHydrogens(cxin, cxout)
    implicit none
    type(cxs) :: cxin, cxout
    integer :: i, j, ic, n
    real(8), dimension(:), allocatable :: xt, yt, zt
    character(len=2), dimension(:), allocatable :: label


    ! Count heteroatoms in cxin.
    !
    n = 0
    do i = 1, cxin%na
      if (trim(cxin%atomlabel(i)) /= 'H') then
        n = n + 1
      endif
    enddo

    allocate(xt(n), yt(n), zt(n), label(n))

    ic = 0
    do i = 1, cxin%na
      if (trim(cxin%atomlabel(i)) /= 'H') then
        ic = ic + 1
        xt(ic) = cxin%r(1, i)
        yt(ic) = cxin%r(2, i)
        zt(ic) = cxin%r(3, i)
        label(ic) = cxin%atomlabel(i)
      endif
    enddo

    call CreateCXS(cxout, n, label, xt, yt, zt)

    deallocate(xt, yt, zt, label)

    return
  end Subroutine RemoveHydrogens


  !
  !*************************************************************************
  !> SetCXSLattice
  !!
  !! Initializes the coordinates of the chemical structure object cxs
  !! on a simple cubic lattice, with lattice spacing determined by the
  !! parameter LATTICESTEP.
  !!
  !! - cx: Input chemical structure object.
  !!
  !*************************************************************************
  !
  Subroutine SetCXSLattice(cx)
    implicit none
    type(cxs) :: cx
    integer :: nlat, ic, i, j, k

    nlat = nint(dble(cx%na)**(1.d0/3.d0))
    if (nlat**3 < cx%na) nlat = nlat + 1

    ic = 0
    outer: do i = 1, nlat
      do j = 1, nlat
        do k = 1, nlat
          ic = ic+1
          if (ic <= cx%na) then
            cx%r(1, ic) = dble(i-1) * LATTICESTEP
            cx%r(2, ic) = dble(j-1) * LATTICESTEP
            cx%r(3, ic) = dble(k-1) * LATTICESTEP
          else
            exit outer
          endif
        enddo
      enddo
    enddo outer

    return
  end Subroutine SetCXSLattice


  !
  !*************************************************************************
  !> OptCXSAgainstGraph
  !!
  !! Optimizes the coordinates of a chemical structure object under
  !! the action of the graph-restraining potential only.
  !!
  !! - cx: Input chemical structure object.
  !! - gdsrestspring: GDS spring restraint strength (au)
  !! - nbstrength: Repulsion interaction strength for non-bonded atoms.
  !! - nbrange: Repulsion interaction range for non-bonded atoms.
  !! - kradius: Spring constant for repulsion between deiscrete molecules
  !!            (au)
  !!
  !*************************************************************************
  !
  Subroutine OptCXSAgainstGraph(cx, gdsrestspring, nbstrength, nbrange, kradius, &
      nrelax, step)
    implicit none
    type(cxs) :: cx
    real(8) :: gdsrestspring, nbstrength, nbrange, kradius
    real(8) :: step, x, y, z
    integer :: it, nrelax, idof, i, k, j

    cx%dvdr(:, :) = 0.d0
    call GraphConstraints(cx, gdsrestspring, nbstrength, nbrange, kradius)
    do it = 1, nrelax
      idof = 0
      do i = 1, cx%na
        if (.not. cx%fixedatom(i)) then
          do k = 1, 3
            idof = idof + 1
            if (.not. cx%fixeddof(idof)) then
                cx%r(k, i) = cx%r(k, i) - step * cx%dvdr(k, i)
            endif
          enddo
        else
          idof = idof + 3
        endif
      enddo
      cx%dvdr(:, :) = 0.d0
      call GraphConstraints(cx, gdsrestspring, nbstrength, nbrange, kradius)

   !    write(43,'(i5)')cx%na
   !    write(43,*)
   !    do j = 1, cx%na
   !       x = cx%r(1,j) * bohr_to_ang
   !       y = cx%r(2,j) * bohr_to_ang
   !       z = cx%r(3,j) * bohr_to_ang
   !       write(43,'(a2,2x,3(f14.8,2x))')cx%atomlabel(j),x,y,z
   !    enddo

    enddo
!    stop

    return
  end Subroutine OptCXSAgainstGraph


  !
  !*************************************************************************
  !> OptimizeGRP
  !!
  !! Minimises the molecules according to the Graph restraining potential
  !!
  !! - cx: The input chemical structure
  !! - success: weather we satisfied the contraints after the monimization or not
  !! - rest: parameters for the GRP and minimization
  !!
  !*************************************************************************
  !
  Subroutine OptimizeGRP(cx, success, gdsrestspring, nbstrength, nbrange, kradius, ngdsrelax, gdsdtrelax)
    implicit none
    type(cxs) :: cx, cxtmp1, cxtmp2
    integer :: i, j, k, l, n, m, it, cc1, cc2, idof, isum, ngdsrelax, iact
    real(8) :: gdsrestspring, nbstrength, nbrange, kradius, gdsdtrelax,sum,rmax
    logical :: success, debug

    debug = .false.
    ! Make a copy of the current structure - cx%graph(:,:) should already contain
    ! the target graph.
    !
    call CopyToNewCXS(cx, cxtmp1)
    call CopyToNewCXS(cx, cxtmp2)
    cc2 = 0 ; cc1 = 0
    cx%dvdr(:, :) = 0.0D0
    call GraphConstraints(cx, gdsrestspring, nbstrength, nbrange, kradius)
  !  if (debug) open(22,file='banana.xyz',status='unknown')
    do it = 1, ngdsrelax
      idof = 0
      sum = 0.d0
      rmax = -1d6
      iact = 0
      do i = 1, cx%na
        if (.not. cx%fixedatom(i)) Then
          do k = 1, 3
            idof = idof + 1
            if (.not. cx%fixeddof(idof)) then
              iact = iact + 1
              cx%r(k, i) = cx%r(k, i) - gdsdtrelax * cx%dvdr(k, i)
              sum = sum + cx%dvdr(k, i)**2
              if (abs( cx%dvdr(k, i)) > rmax) rmax = abs(cx%dvdr(k, i))
            endif
          enddo
        else
          idof = idof + 3
        endif
      enddo
      cx%dvdr(:, :) = 0.0D0
      cc1 = cc1 + 1
      if (mod(cc1, 200) .eq. 0) then
        call CopyCXS(cx, cxtmp2)
        call GetGraph(cxtmp2)
        isum = 0
        do k = 1, cx%na
          do i = k, cx%na
            isum = isum + abs(cxtmp2%graph(k, i) - cxtmp1%graph(k, i))
            !if (abs(cxtmp2%graph(k,i) - cxtmp1%graph(k,i)) .ne. 0) &
            ! print*, i,k,'ATOMS '//trim(cxtmp2%atomlabel(i))//' '//&
            !   trim(cxtmp2%atomlabel(k))//' are = ',cxtmp2%graph(k,i), 'SHOULD BE ', cxtmp1%graph(k,i)
          enddo
        enddo
        if (isum == 0 ) cc2  = cc2 + 1
      !write(22,*) cx%na
      !write(22,*) ' '
      !do i = 1, cx%na
      !  write(22,'("'//cxtmp2%atomlabel(i)//'",3(X,F15.7))') cx%r(1:3,i)*bohr_to_ang
      !enddo
      !flush(22)
      endif
  !   if (debug) then
  !    write(22,*) cx%na
  !    write(22,*) 'FORCE? = ', norm2(reshape(cx%dvdr(1:3,1:cx%na),(/cx%na*3/)))
  !    do i = 1, cx%na
  !      write(22,'("'//cxtmp2%atomlabel(i)//'",3(X,F15.7))') cx%r(1:3,i)*bohr_to_ang
  !    enddo
  !    flush(22)
      !print*, 'FORCE? = ', norm2(reshape(cx%dvdr(1:3,1:cx%na),(/cx%na*3/)))
  !   endif
      call GraphConstraints(cx, gdsrestspring, nbstrength, nbrange, kradius)
      !print*, 'CCX = ', cc2
      if (cc2 > 3 ) then
      success = .true.
      return
      endif
    enddo
    success = .false.
    !Call CopyCXS(cx, cxtmp2)
    !Call GetGraph( cxtmp2 )
    !isum = 0
    !do k = 1, cx%na
    !   do i = 1, cx%na
    !      isum = isum + abs(cxtmp2%graph(k,i) - cxtmp1%graph(k,i))
    !      if (abs(cxtmp2%graph(k,i) - cxtmp1%graph(k,i)) .ne. 0) then
    !       print*, 'isum = ', isum
    !       print*, 'ATOM ', cxtmp2%atomlabel(k), 'AND ', cxtmp2%atomlabel(i), 'SHOULD BE', &
    !       cxtmp2%graph(k,i), cxtmp1%graph(k,i)
    !      endif
    !   enddo
    !enddo
    !if (isum == 0 ) cc2  = cc2 + 1
    !stOP
    !close(22)
    return
  end Subroutine OptimizeGRP


  !
  !*************************************************************************
  !> OptimizeGRP2
  !!
  !! Minimises the molecules according to the Graph restraining potential
  !!
  !! - cx: The input chemical structure
  !! - success: weather we satisfied the contraints after the monimization or not
  !! - rest: parameters for the GRP and minimization
  !!
  !*************************************************************************
  !
  Subroutine OptimizeGRP2(cx, cxstart, success, gdsrestspring, nbstrength, nbrange, kradius, ngdsrelax, gdsdtrelax)
    implicit none
    type(cxs) :: cx, cxstart, cxtmp1, cxtmp2
    integer :: i, j, k, l, n, m, it, cc1, cc2, idof, isum, ngdsrelax, iact
    real(8) :: gdsrestspring, nbstrength, nbrange, kradius, gdsdtrelax, sum, rmax
    logical :: success, debug

    debug = .false.
    ! Make a copy of the current structure - cx%graph(:,:) should already contain
    ! the target graph.
    !
    call CopyToNewCXS(cx, cxtmp1)
    call CopyToNewCXS(cx, cxtmp2)
    cc2 = 0 ; cc1 = 0
    cx%dvdr(:, :) = 0.0D0
    !Call GraphConstraints( cx, gdsrestspring, nbstrength, nbrange, kradius )
    call GraphConstraints_DoubleEnded(cx, cxstart, gdsrestspring,nbstrength, nbrange, kradius)
    if (debug) open(22, file='banana.xyz', status='unknown', position='append')
    do it = 1, ngdsrelax
      idof = 0
      sum = 0.d0
      rmax = -1d6
      iact = 0
      do i = 1, cx%na
        if (.not. cx%fixedatom(i)) Then
          do k = 1, 3
            idof = idof + 1
            if (.not. cx%fixeddof(idof)) then
              iact = iact + 1
              cx%r(k, i) = cx%r(k, i) - gdsdtrelax * cx%dvdr(k, i)
              sum = sum + cx%dvdr(k, i)**2
              if (abs(cx%dvdr(k, i)) > rmax)rmax = abs(cx%dvdr(k, i))
            endif
          enddo
        else
          idof = idof + 3
        endif
      enddo
      cx%dvdr(:, :) = 0.0D0
      cc1 = cc1 + 1
      if ( mod(cc1,200) .eq. 0 ) then
        call CopyCXS(cx, cxtmp2)
        call GetGraph(cxtmp2)
        isum = 0
        do k = 1, cx%na
          do i = k, cx%na
            isum = isum + abs(cxtmp2%graph(k, i) - cxtmp1%graph(k, i))
            !if (abs(cxtmp2%graph(k,i) - cxtmp1%graph(k,i)) .ne. 0) &
            ! print*, i,k,'ATOMS '//trim(cxtmp2%atomlabel(i))//' '//&
            !   trim(cxtmp2%atomlabel(k))//' are = ',cxtmp2%graph(k,i), 'SHOULD BE ', cxtmp1%graph(k,i)
          enddo
        enddo
        if (isum == 0 ) cc2  = cc2 + 1
      !write(22,*) cx%na
      !write(22,*) ' '
      !do i = 1, cx%na
      !  write(22,'("'//cxtmp2%atomlabel(i)//'",3(X,F15.7))') cx%r(1:3,i)*bohr_to_ang
      !enddo
      !flush(22)
      endif
  !   if (debug) then
  !    write(22,*) cx%na
  !    write(22,*) 'FORCE? = ', norm2(reshape(cx%dvdr(1:3,1:cx%na),(/cx%na*3/)))
  !    do i = 1, cx%na
  !      write(22,'("'//cxtmp2%atomlabel(i)//'",3(X,F15.7))') cx%r(1:3,i)*bohr_to_ang
  !    enddo
  !    flush(22)
      !print*, 'FORCE? = ', norm2(reshape(cx%dvdr(1:3,1:cx%na),(/cx%na*3/)))
  !   endif
      !Call GraphConstraints( cx, gdsrestspring, nbstrength, nbrange, kradius )
      call GraphConstraints_DoubleEnded(cx, cxstart, gdsrestspring,nbstrength, nbrange, kradius)
      !print*, 'CCX = ', cc2
      if (cc2 > 3 ) then
      success = .true.
      return
      endif
    enddo
    success = .false.
    !Call CopyCXS(cx, cxtmp2)
    !Call GetGraph( cxtmp2 )
    !isum = 0
    !do k = 1, cx%na
    !   do i = 1, cx%na
    !      isum = isum + abs(cxtmp2%graph(k,i) - cxtmp1%graph(k,i))
    !      if (abs(cxtmp2%graph(k,i) - cxtmp1%graph(k,i)) .ne. 0) then
    !       print*, 'isum = ', isum
    !       print*, 'ATOM ', cxtmp2%atomlabel(k), 'AND ', cxtmp2%atomlabel(i), 'SHOULD BE', &
    !       cxtmp2%graph(k,i), cxtmp1%graph(k,i)
    !      endif
    !   enddo
    !enddo
    !if (isum == 0 ) cc2  = cc2 + 1
    !stOP
    !close(22)
    return
  end Subroutine OptimizeGRP2


  !
  !*************************************************************************
  !> OptimizeGRP_DoubleEnded
  !!
  !! Optimizes the geometry of cx according to a modified Graph restraining
  !! potential which tries to keep geometry close to a starting point cxstart.
  !!
  !! - cx: The input chemical structure
  !! - success: weather we satisfied the contraints after the monimization or not
  !! - rest: parameters for the GRP and minimization
  !!
  !*************************************************************************
  !
  Subroutine OptimizeGRP_DoubleEnded(cx, cxstart, success, gdsrestspring, nbstrength, nbrange, kradius, ngdsrelax, gdsdtrelax)
    implicit none
    type(cxs) :: cx, cxstart
    integer :: i, j, k, l, n, m, it, cc1, cc2, idof, isum, ngdsrelax, iact
    real(8) :: gdsrestspring, nbstrength, nbrange, kradius, gdsdtrelax, sum, rmax
    logical :: success, debug

    debug = .false.
    if (debug) open(21, file='ori.xyz', status='unknown', position='append')
    if (debug) then
      write(21, *) cx%na
      write(21, *) 'FORCE? = ', norm2(reshape(cx%dvdr(1:3,1:cx%na), (/cx%na*3/)))
      do i = 1, cx%na
        write(21, '("'//cx%Atomlabel(i)//'", 3(X, F15.7))') cx%r(1:3, i)*bohr_to_ang
      enddo
      close(21)
    endif
    !if (debug) open(22,file='banana.xyz',status='unknown',position='append')

    ! SD optimization.
    !
    success = .false.
    outer: do it = 0, ngdsrelax
      idof = 0
      sum = 0.d0
      rmax = -1d6
      iact = 0

      if (it > 0) then

        do i = 1, cx%na
          if (.not. cx%fixedatom(i)) Then
            do k = 1, 3
              idof = idof + 1
              if (.not. cx%fixeddof(idof)) then
                iact = iact + 1
                cx%r(k, i) = cx%r(k, i) - gdsdtrelax * cx%dvdr(k, i)
                sum = sum + cx%dvdr(k, i)**2
                if (abs( cx%dvdr(k, i) ) > rmax)rmax = abs(cx%dvdr(k, i))
              endif
            enddo
          else
            idof = idof + 3
          endif
        enddo

        ! Check convergence based on RMS forces and maximum force.
        !
        sum = dsqrt(sum / dble(idof) )

        if (sum < GRPMINTHRESH .and. rmax < GRPMAXTHRESH) then
          success = .true.
          exit outer
        endif
      endif

      ! Update derivatives.
      !
      cx%dvdr(:, :) = 0.d0
      ! the nbrange should be shorter than with other GRP approaches...hence the
      ! 0.5 factor..?
      call GraphConstraints_DoubleEnded(cx, cxstart, gdsrestspring, nbstrength, nbrange, kradius)
      if (debug) then
        write(22, *) cx%na
        write(22, *) 'FORCE? = ', norm2(reshape(cx%dvdr(1:3, 1:cx%na), (/cx%na*3/)))
        do i = 1, cx%na
          write(22, '("'//cx%Atomlabel(i)//'", 3(X, F15.7))') cx%r(1:3, i)*bohr_to_ang
        enddo
        flush(22)
        !print*, 'FORCE? = ', norm2(reshape(cx%dvdr(1:3,1:cx%na),(/cx%na*3/)))
      endif

    enddo outer
    if (debug) close(22)

    return
  end Subroutine OptimizeGRP_DoubleEnded


  !
  !************************************************************************
  !> MolecularFormula
  !!
  !! Returns a variable length string with the molecular formula of molecule
  !! molid, from chemical structure cx.
  !!
  !! - cx: chemical structure
  !! - molecule index
  !!
  !************************************************************************
  !
  function MolecularFormula(cx, molid)
    implicit none
    integer      :: j, k, i, molid, itmp
    integer,parameter :: nelem = 100
    type(cxs)    :: cx
    integer      :: MolComp(nelem)
    character(2) :: AtomLab(nelem), ctmp
    character(20):: ctmp2, MolForm
    character(:), allocatable :: MolecularFormula
    logical      :: swapped = .false.
    
    MolComp = 0
    AtomLab = ''
    do i = 1, cx%namol(molid)
      j = LabelToNumber(cx%atomlabel(cx%molid(molid, i)))
      MolComp(j) = MolComp(j) + 1
      if (AtomLab(j) .eq. '' ) AtomLab(j) = cx%atomlabel(cx%molid(molid, i))
    enddo
    ! Bubble sorts the Molecular Formula with the Highest Number of that element last
    !
    do j = nelem-1, 1, -1
      swapped = .false.
      do i = 1, j
        if (MolComp(i) > MolComp(i+1)) then
          itmp = MolComp(i)
          ctmp = AtomLab(i)
          MolComp(i) = MolComp(i+1)
          AtomLab(i) = AtomLab(i+1)
          MolComp(i+1) = itmp
          AtomLab(i+1) = ctmp
          swapped = .true.
        endif
      enddo
      if (.not. swapped) exit
    enddo
    MolForm = ''
    i = nelem
    do while ( MolComp(i) .ne. 0 )
      write(ctmp2,'(I8)') int(MolComp(i))
      MolForm = trim(adjustl(MolForm))//trim(adjustl(AtomLab(i)))//trim(adjustl(ctmp2))
      i = i - 1
    enddo
    allocate(character(len=len(MolForm)) :: MolecularFormula)
    MolecularFormula = MolForm
    return
  end function MolecularFormula


  !
  !*************************************************************************
  !> GetAtomValency
  !!
  !! returns the valency of every atom in cx, the first index (2,:)
  !! only counts bonds which can actually change during the gds calculation
  !!
  !! - cx: The input chemical structure
  !!
  !*************************************************************************
  !
  Function GetAtomValency(cx, fixedbonds) result(val)
    implicit none
    type(cxs) :: cx
    logical  :: fixedbonds(:, :)
    integer :: i, j, k, l, n, m, idof, sum, naa
    integer :: val(2,cx%na)

    val = 0
    do i = 1, cx%na-1
      do j = i+1, cx%na
        if (cx%graph(i, j) .eq. 1) then
          if (.not. fixedbonds(i, j) ) then
            if (cx%fixedatom(i) .and. .not. cx%fixedatom(j)) val(2, i) = val(2, i) + 1
            if (cx%fixedatom(j) .and. .not. cx%fixedatom(i)) val(2, j) = val(2, j) + 1
          endif
          val(1, i) = val(1, i) + 1
          val(1, j) = val(1, j) + 1
        endif
      enddo
    enddo
    return
  end Function GetAtomValency


  !
  !*************************************************************************
  !> GetElementPairValency
  !!
  !! returns a matrix of size val(na,na) with the number of bonds of a particular pair
  !! of elements (eg H-C) that corresponds to that pair of atoms va(i,j), where the
  !! first index says how many of those there are bonded to atom i of type j
  !! used for rxnvalatom tests
  !!
  !! - cx: The input chemical structure
  !!
  !*************************************************************************
  !
  Function GetElementPairValency(cx) result(val)
    implicit none
    type(cxs) :: cx
    integer :: i, j, k, l, n, m, idof, naa
    integer  :: val(cx%na, cx%na)
    val = 0
    naa = cx%na
    do k = 1, nrxval
      do i = 1, naa
        if (cx%atomlabel(i) == trim(rxvalatom(k, 1))) then
          do j = 1, naa ; if (val(i, j) .ne. 0) cycle
            if (cx%atomlabel(j) == trim(rxvalatom(k, 2))) then
              do l = 1, naa ; if (val(i, l) .ne. 0 .or. i .eq. l) cycle
                if(cx%atomlabel(j) == cx%atomlabel(l)) val(i, j) = val(i, j) + cx%graph(i, j)
              enddo
              do l = j+1, naa ; if (val(i,l) .ne. 0 .or. i .eq. l) cycle
                if(cx%atomlabel(j) == cx%atomlabel(l)) val(i, l) = val(i, j)
              enddo
            endif
          enddo
        endif
      enddo
    enddo
    return
  end Function


  !
  !*************************************************************************
  !> AllowedCXSValenceRange
  !!
  !! Checks wheather the valence of elements falls withtin the allowed range
  !!
  !! - cx: The input chemical structure
  !!
  !*************************************************************************
  !
  Function AllowedCXSValenceRange(cx)
    use globaldata
    implicit none
    type(cxs) :: cx
    integer :: i, j, k, l, n, m, idof, sum, naa
    logical :: AllowedCXSValenceRange
    naa = cx%na
    do k = 1, nvalcon
      do i = 1, naa
        if (trim(cx%atomlabel(i)) == trim(valatom(k))) then
          if (.not. valfz(k)) then
            sum = 0
            do j = 1, naa
              if (i/=j) then
                  sum = sum + cx%graph(i, j)
              endif
            enddo
          else
            ! this only counts bonds which can actually change during the gds calculation
            sum = 0
            do j = 1, naa
              if (i/=j) then
                if (cx%fixedatom(i) .and. .not. cx%fixedatom(j) .and. &
                    .not. fixedbonds(i, j) ) then
                  sum = sum + cx%graph(i, j)
                endif
              endif
            enddo
          endif
          if (sum < valrange(k, 1) .or. sum > valrange(k, 2))then
            AllowedCXSValenceRange = .false.
            !print*, 'FAILED  - ', sum,  valrange(k,1),  valrange(k,2),trim(valatom(k)), trim(cx%atomlabel(i))
            return
          endif
        endif
      enddo
    enddo
    AllowedCXSValenceRange = .true.
    return
  end Function AllowedCXSValenceRange


  !
  !*************************************************************************
  !> AllowedCXSReactiveValence
  !!
  !! Checks wheather the valence of an element to some other falls withing the allowed values
  !!
  !! - cx: The input chemical structure
  !!
  !*************************************************************************
  !
  Function AllowedCXSReactiveValence(cx)
    use globaldata
    implicit none
    type(cxs) :: cx
    integer :: i, j, k, l, n, m, idof, sum, naa
    logical :: AllowedCXSReactiveValence

    naa = cx%na
    do i = 1, naa
      do k = 1, nrxval
        if (cx%atomlabel(i) == trim(rxvalatom(k, 1))) then
          ! Sum up bonded atoms.
          !
          sum = 0
          do j = 1, naa
            if (cx%atomlabel(j) == trim(rxvalatom(k, 2))) then
              if (i/=j) then
                sum = sum + cx%graph(i, j)
              endif
            endif
          enddo
          if (sum < rxvalrange(k, 1) .or. sum > rxvalrange(k, 2))then
            AllowedCXSReactiveValence = .false.
            return
          endif
        endif
      enddo
    enddo
    AllowedCXSReactiveValence = .true.
    return
  end Function AllowedCXSReactiveValence


  !
  !*************************************************************************
  !> AllowedCXSBondsMax
  !!
  !! Checks wheather the bonds present are allowed as prescibed by allowedbonds keyword
  !!
  !! - cx: The input chemical structure
  !!
  !*************************************************************************
  !
  Function AllowedCXSBondsMax(cx)
    use globaldata
    implicit none
    type(cxs) :: cx
    integer :: i, j, k, l, n, m, idof, sum, naa
    logical :: AllowedCXSBondsMax

    naa = cx%na
    do k = 1, nallowbonds
      do i = 1, naa
        if (trim(cx%atomlabel(i)) == trim( allowbondsatom(k,1) )) then
          sum = 0
          do j = 1, naa
            if (i/=j .and. (trim(cx%atomlabel(j)) == trim(allowbondsatom(k, 2))) ) then
              sum = sum + cx%graph(i, j)
            endif
          enddo
          if (sum > allowbondsmax(k)) then
            AllowedCXSBondsMax = .false.
            return
          endif
        endif
      enddo
    enddo
    AllowedCXSBondsMax = .true.
    return
  end Function AllowedCXSBondsMax


  !
  !***************************************************************************
  !> SetReactiveIndices
  !!
  !! Sets the bondchange(:,:) and atomchange(:) indices, which indicate
  !! whether the specified bonds or atoms are allowed to change as a
  !! result of a proposed graph move.
  !!
  !! Returns the bondchange(:,:) and atomchange(:) variables, as well
  !! as the number of reactive atoms (nrx) and a list of the reactive
  !! atoms (in rxindex(1:nrx)).
  !!
  !! bondchange(:,:) - Logical array indicating whether bond (i,j) is
  !!                   allowed to change or not.
  !! atomchange(:) - Logical array indicating whether or not atom i is
  !!                 allowed to change its bonding pattern.
  !! na - Number of atoms
  !! cx - A cxs object, primarily used to provide atom labels.
  !! rxindex - List of reactive atoms.
  !! nrx - Number of reactive atoms.
  !!
  !***************************************************************************
  !
  Subroutine SetReactiveIndices(bondchange, atomchange, naa, cx, rxindex, nrx)
    use globaldata
    implicit none
    integer :: i, j, k
    integer :: naa, rxindex(NAMAX), nrx, sum
    type(cxs) :: cx
    logical :: atomchange(naa), bondchange(naa,naa)
    character (len=2) :: id1, id2, i1, i2

    ! Set initial flags for atomchange and bondchange.
    !
    atomchange(:) = .FALSE.
    bondchange(:,:) = .TRUE.

    ! Use the information from the input file to define which atoms
    ! can react, and which bonds can change.
    !
    ! First, check reactive atom types...
    !
    do i = 1, naa
      do j = 1, nreactivetypes
        if (reactivetype(j) == cx%atomlabel(i)) then
          atomchange(i) = .TRUE.
          !print*, 'TRUE = ', i
        endif
      enddo
    enddo

    ! Second, check reactive atom ranges and ids.
    !
    do i = 1, naa
       atomchange(i) = reactive(i)
    enddo

    ! Third, fix indicated bonds.
    !
    do i = 1, naa
      id1 = cx%atomlabel(i)
      do j = 1, naa
        id2 = cx%atomlabel(j)
        do k = 1, nfixtype
          i1 = fixedbondtype(k, 1)
          i2 = fixedbondtype(k, 2)
          if (trim(id1) == trim(i1) .and. trim(id2) == trim(i2)) then
            bondchange(i, j) = .FALSE.
          else if (trim(id1) == trim(i2) .and. trim(id2) == trim(i1)) then
            bondchange(i, j) = .FALSE.
          endif
        enddo

        ! NEW - fix bonds
        if (fixedbonds(i, j)) then
          bondchange(i, j) = .FALSE.
        endif

      enddo
      bondchange(i, i) = .FALSE.
    enddo

    ! THIS LEADS TO ATOMS GETTING STUCK!
    !! Finaally, check the valence constraints of the atoms..
    !!
    !do i = 1, naa
    !
    !   do k = 1, nrxval
    !      if (cx%atomlabel(i) == trim(rxvalatom(k,1))) then
    !
    !         ! Sum up bonded atoms.
    !         !
    !         sum = 0
    !         do j = 1, naa
    !            if (cx%atomlabel(j) == trim(rxvalatom(k,2))) then
    !               if (i/=j) then
    !                  sum = sum + cx%graph(i,j)
    !               endif
    !            endif
    !         enddo
    !         if (sum < rxvalrange(k,1) .or. sum > rxvalrange(k,2))then
    !            atomchange(i) = .FALSE.
    !         endif
    !      endif
    !   enddo
    !enddo


    ! Set up an array of reactive atoms to help with random selection.
    !
    nrx = 0
    do i = 1, naa
      if (atomchange(i)) then
        nrx = nrx + 1
        rxindex(nrx) = i
      endif
    enddo

    return
  end Subroutine SetReactiveIndices


  !
  !***************************************************************************
  !> CheckForbidden
  !!
  !! Checks whether the current graph, generated during the GraphMoves,
  !! routine, matches one of the forbidden patterns given in the
  !! fobidfile.
  !!
  !! This subroutine is a nightmare-ish hellscape of enddo and endif...
  !! nobody should have to deal with this...
  !!
  !! Returns iflag = 0 if the graph is forbidden.
  !!
  !***************************************************************************
  !
  Subroutine CheckForbidden(cx, nforbid, naforbid, gforbid, forbidlabel, iflag)
    implicit none
    type(cxs) :: cx
    integer :: i, j, k, l, m, n, p, iflag, na
    integer :: nforbid, naforbid(NFORBIDMAX), gforbid(NFORBIDMAX,NAMOVEMAX,NAMOVEMAX)
    character (len=4) :: forbidlabel(NFORBIDMAX,NAMOVEMAX)

    na = cx%na
    ! Loop over each forbidden graph
    !
    do i = 1, nforbid
      ! Decide if any of the bonding graphs match the forbidden pattern.
      !
      !
      !** 2-atom graphs.
      !
      if (naforbid(i) == 2) then
        a1: do j = 1, na
          if (trim(cx%atomlabel(j)) == trim(forbidlabel(i, 1))) then
            a2: do k = 1, na
              if (j/=k)cycle a2
              if (trim(cx%atomlabel(k)) == trim(forbidlabel(i, 2))) then
                if (cx%graph(j, k) == gforbid(i, 1, 2)) then
                  iflag = 0
                  return
                endif
              endif
            enddo a2
          endif
        enddo a1

      !** 3-atom graphs
      !
      else if (naforbid(i) == 3) then
        b1: do j = 1, na
          if (trim(cx%atomlabel(j)) == trim(forbidlabel(i, 1))) then
            b2: do k = 1, na
              if (k == j) cycle b2
              if (trim(cx%atomlabel(k)) == trim(forbidlabel(i, 2))) then
                b3: do l = 1, na
                  if (j == l .or. k == l) cycle b3
                  if (trim(cx%atomlabel(l)) == trim(forbidlabel(i, 3))) then
                    if (cx%graph(j, k) == gforbid(i, 1, 2)) then
                        if (cx%graph(j, l) == gforbid(i, 1, 3)) then
                          if (cx%graph(k, l) == gforbid(i, 2, 3)) then
                            iflag = 0
                            return
                          endif
                        endif
                    endif
                  endif
                enddo b3
              endif
            enddo b2
          endif
        enddo b1

      !** 4-atom graphs
      !
      else if (naforbid(i) == 4) then
        c1: do j = 1, na
          if (trim(cx%atomlabel(j)) == trim(forbidlabel(i, 1))) then

            c2: do k = 1, na
              if (j == k)cycle c2
              if (trim(cx%atomlabel(k)) == trim(forbidlabel(i, 2))) then

                c3: do l = 1, na
                  if (j == l .or. k == l) cycle c3
                  if (trim(cx%atomlabel(l)) == trim(forbidlabel(i, 3))) then

                    c4: do m = 1, na
                      if (j==m .or. k ==m .or. l == m)cycle c4
                      if (trim(cx%atomlabel(m)) == trim(forbidlabel(i, 4))) then

                        if (cx%graph(j, k) == gforbid(i, 1, 2)) then
                          if (cx%graph(j, l) == gforbid(i, 1, 3)) then
                            if (cx%graph(j, m) == gforbid(i, 1, 4)) then
                              if (cx%graph(k, l) /= gforbid(i, 2, 3)) then
                                if (cx%graph(k, m) /= gforbid(i, 2, 4)) then
                                  if (cx%graph(l, m) /= gforbid(i, 3, 4)) then
                                    iflag = 0
                                    return
                                  endif
                                endif
                              endif
                            endif
                          endif
                        endif
                      endif
                    enddo c4
                  endif
                enddo c3
              endif
            enddo c2
          endif
        enddo c1

      !** 5-atom graphs
      !
      else if (naforbid(i) == 5) then
        d1: do j = 1, na
          if (trim(cx%atomlabel(j)) == trim(forbidlabel(i,1))) then

            d2: do k = 1, na
              if (j == k)cycle d2
              if (trim(cx%atomlabel(k)) == trim(forbidlabel(i,2))) then

                d3: do l = 1, na
                  if (j == l .or. k == l) cycle d3
                  if (trim(cx%atomlabel(l)) == trim(forbidlabel(i,3))) then

                    d4: do m = 1, na
                      if (j==m .or. k ==m .or. l == m)cycle d4
                      if (trim(cx%atomlabel(m)) == trim(forbidlabel(i,4))) then

                        d5: do p = 1, na
                            if (j == p .or. k == p.or.l==p.or.m==p)cycle d5
                            if (trim(cx%atomlabel(m)) == trim(forbidlabel(i,5))) then

                              if (cx%graph(j,k) == gforbid(i,1,2)) then
                                if (cx%graph(j,l) == gforbid(i,1,3)) then
                                  if (cx%graph(j,m) == gforbid(i,1,4)) then
                                    if (cx%graph(j,p) == gforbid(i,1,5)) then
                                      if (cx%graph(k,l) /= gforbid(i,2,3)) then
                                        if (cx%graph(k,m) /= gforbid(i,2,4)) then
                                          if (cx%graph(k,p) /= gforbid(i,2,5)) then
                                            if (cx%graph(l,m) /= gforbid(i,3,4)) then
                                              if (cx%graph(l,p) /= gforbid(i,3,5)) then
                                                if (cx%graph(m,p) /= gforbid(i,4,5)) then
                                                  iflag = 0
                                                  return
                                                endif
                                              endif
                                            endif
                                          endif
                                        endif
                                      endif
                                    endif
                                  endif
                                endif
                              endif
                            endif
                          enddo d5
                        endif
                      enddo d4
                    endif
                enddo d3
              endif
            enddo d2
          endif
        enddo d1

      else
        stop '* ERROR: Forbidden graph patterns only implemented for up to 4 atoms'
      endif

    enddo

    return
  end Subroutine CheckForbidden


  !
  !***************************************************************************
  !> PrintCXSGraphInfo()
  !!
  !! Outputs the bonding graph information for a chemical structure object.
  !!
  !! - cx: The chemical structure object
  !! - iunit: Integer ID of the output destination.
  !! - message: Character string (80) identifying CXS details.
  !!
  !***************************************************************************
  !
  Subroutine PrintCXSGraphInfo(cx, iunit, message)
    use globaldata
    implicit none
    type(cxs) :: cx          !< Chemical structure object
    integer :: iunit         !< Output unit ID
    character(len=*) :: message  !< A message identifying this CXS.
    integer :: i, j, isum
    character*3 :: ichar1, ichar2

    write(iunit, '(/"================================================================")')
    write(iunit, '("# Graph Info for CXS:",2x,a/)') adjustl(trim(message))

    write(iunit, '("* Bonded atoms *")')
    do i = 1, cx%na-1
      do j = i+1, cx%na
        if (cx%graph(i, j) == 1) then
          write(ichar1, '(i3)') i
          write(ichar2, '(i3)') j
          write(iunit, '(1x, a, a, 1x, "-", 1x, a, a)') cx%atomlabel(i), adjustl(ichar1), cx%atomlabel(j), adjustl(ichar2)
        endif
      enddo
    enddo

    write(iunit, '(/"* Atomic valencies *")')
    do i = 1, cx%na
      isum = 0
      do j = 1, cx%na
        if (i /= j) then
          isum = isum + cx%graph(i, j)
        endif
      enddo
      write(ichar1, '(i3)') i
      write(iunit, '(1x, a, a, 2x, ":", 1x, i4)') cx%atomlabel(i), adjustl(ichar1), isum
    enddo
    write(iunit, '(/"================================================================")')

    return
  end Subroutine PrintCXSGraphInfo


  !
  !***************************************************************************
  !> SetValenceCoords()
  !!
  !! Calculates the set of internal coordinates (bond-lengths, bond-angles and
  !! torsion angles) for cx. Note that the derivatives of each internal
  !! coordinates with respect to the atomic coordinates is also calculated.
  !!
  !! Important: The identification of bonds, angle and torsion angles is
  !! based on the cx%graph(:,:).
  !!
  !! - cx: The chemical structure object
  !!
  !***************************************************************************
  !
  Subroutine SetValenceCoords(cx)
    implicit none
    type(cxs) :: cx          !< Chemical structure object

    ! Zero everything.
    !
    cx%nangles = 0
    cx%ntors = 0

    cx%angle(:) = 0.d0
    cx%torsion(:) = 0.d0

    cx%dangdr(:, :, :) = 0.d0
    cx%dtorsdr(:, :, :) = 0.d0

    cx%angleid(:, :) = 0
    cx%torsid(:, :) = 0

    ! First, get all of the bond-lengths
    !
    call GetBonds(cx)

    ! Next, Get all angles.
    !
    call GetAngles(cx)

    ! Finally, Get all torsion angles.
    !
  !  Call GetTorsions(cx)

    return
  end Subroutine SetValenceCoords


  !
  !***************************************************************************
  !> GetBonds()
  !!
  !! Calculates the set of bond-lengths for cx.
  !! Note that the derivatives of each bond-length
  !! with respect to the atomic coordinates is also calculated.
  !!
  !! Important: The identification of bonds is
  !! based on the cx%graph(:,:).
  !!
  !! - cx: The chemical structure object
  !!
  !***************************************************************************
  !
  Subroutine GetBonds(cx)
    implicit none
    type(cxs) :: cx
    integer :: i, j, na
    real(8) :: dd(3), rr, onr

    ! Zero everything.
    !
    cx%nbonds = 0
    cx%bondl(:) = 0.d0
    cx%dbonddr(:, :, :) = 0.d0
    cx%bondid(:, :) = 0

    ! Loop over all atoms and identify bonds.
    !
    na = cx%na
    do i = 1, na - 1
      do j = i + 1, na

        if (cx%graph(i, j) == 1) then
          cx%nbonds = cx%nbonds + 1
          cx%bondid(cx%nbonds, 1) = i
          cx%bondid(cx%nbonds, 2) = j

          dd(1:3) = cx%r(1:3, i) - cx%r(1:3, j)
          rr = dsqrt(dd(1)*dd(1) + dd(2)*dd(2) + dd(3)*dd(3))
          onr = 1.d0 / rr
          cx%bondl(cx%nbonds) = rr

          cx%dbonddr(cx%nbonds, 1, 1:3) = dd(1:3) * onr
          cx%dbonddr(cx%nbonds, 2, 1:3) = -dd(1:3) * onr
        endif
      enddo
    enddo

    return
  end Subroutine GetBonds


  !
  !***************************************************************************
  !> GetAngles()
  !!
  !! Calculates the set of bond-angles for cx.
  !! Note that the derivatives of each bond-angle
  !! with respect to the atomic coordinates is also calculated.
  !!
  !! Important: The identification of bonds is
  !! based on the cx%graph(:,:).
  !!
  !! - cx: The chemical structure object
  !!
  !***************************************************************************
  !
  Subroutine GetAngles(cx)
    implicit none
    type(cxs) :: cx
    integer :: i, j, na, k, ii, jj, kk, ifound
    real(8) :: dd(3), rr, onr

    ! Zero everything
    !
    cx%nangles = 0
    cx%angle(:) = 0.d0
    cx%dangdr(:, :, :) = 0.d0
    cx%angleid(:, :) = 0

    ! Set number of atoms.
    !
    na = cx%na

    ! Loop over all bonded pairs.
    !
    do i = 1, na-1
      do j = i+1, na

        if (cx%graph(i, j) == 1) then

          ! Loop over atoms, looking for third partner in bond-angle.
          !
          do k = 1, na
            if (k /= i .and. k /= j) then
              ifound = 0

              ! Case where i is bonded to k and j, i in middle.
              !
              if (cx%graph(i, k) == 1) then
                ifound = 1
                cx%nangles = cx%nangles + 1
                cx%angleid(cx%nangles, 1) = j
                cx%angleid(cx%nangles, 2) = i
                cx%angleid(cx%nangles, 3) = k

              ! Case where j is bonded to k, j in middle.
              else if (cx%graph(j, k) == 1) then
                ifound = 1
                cx%nangles = cx%nangles + 1
                cx%angleid(cx%nangles, 1) = i
                cx%angleid(cx%nangles, 2) = j
                cx%angleid(cx%nangles, 3) = k

              endif

              ! Work out the angle and the derivatives.
              !
              if (ifound == 1) then
                ii = cx%angleid(cx%nangles, 1)
                jj = cx%angleid(cx%nangles, 2)
                kk = cx%angleid(cx%nangles, 3)
              endif

            endif
          enddo
        endif
      enddo
    enddo
    return
  end Subroutine GetAngles


  !
  !***************************************************************************
  !> GetMolecularGraphVector
  !!
  !! fobidfile.
  !!
  !! Returns iflag = 0 if the graph is forbidden.
  !!
  !***************************************************************************
  !
  function GetMolecularGraphVector(cx, mi) result(vec)
    implicit none
    integer      :: j, k, i, molid, itmp, mi
    type(cxs)    :: cx
    integer      :: vec(cx%namol(mi)*(cx%namol(mi)-1)/2), mat(cx%namol(mi), cx%namol(mi))
    
    do i = 1, cx%namol(mi)
      do j = 1, cx%namol(mi)
        mat(i, j) = cx%graph(cx%molid(mi, i),cx%molid(mi, j))
      enddo
    enddo
    k = 0
    do i = 1, cx%namol(mi)-1
      do j = i+1, cx%namol(mi)
        k = k + 1
        vec(k) = mat(i, j)
      enddo
    enddo
    return
  end function


  !
  !************************************************************************
  !> GetMolecularSmiles
  !!
  !! Returns MolSmiles, an array of characters with the Smiles strip for each molecule in MolecAct
  !!
  !! - cx: Chemical structure
  !! - nactmol: number of active molecules in cx
  !! - MolecAct: array with indecies of the nactmol molecules that are active in this problem
  !! - MolSmiles: array of nactmol Smiles strips
  !
  !************************************************************************
  !
  Function GetMolecularSmiles(cx, mi) result(MolSmiles)
    implicit none
    type(cxs) :: cx
    integer :: i, j, k, mi
    character(250) :: MolSmiles
    double precision :: x, y, z

    open(14, file='tmp.xyz', status='unknown')
    write(14, '(i5)') cx%namol(mi)
    write(14, '(A10, i2, A1)') '(Molecule ', mi, ')'
    do j = 1, cx%namol(mi)
      x = cx%r(1, cx%molid(mi, j)) * bohr_to_ang
      y = cx%r(2, cx%molid(mi, j)) * bohr_to_ang
      z = cx%r(3, cx%molid(mi, j)) * bohr_to_ang
      write(14, '(a2, 2x, 3(f14.8, 2x))') cx%atomlabel(cx%molid(mi,j)), x, y, z
    enddo
    flush(14)
    !if (i .eq. 1) call system("touch tmp_smile ; rm tmp_smile ; touch tmp_smile")
    call system("babel -ixyz tmp.xyz -osmi --canonical -O tmp &>/dev/null ; cat tmp >> tmp_smile ")
    close(14, status="delete")
    open(14, file='tmp_smile', status='unknown')
    read(14, *) MolSmiles
    close(14, status="delete")
    return
  end Function GetMolecularSmiles


  !
  !************************************************************************
  !> MolecularMass
  !!
  !! Returns the mass of a molecule molid in cx
  !!
  !! - cx: Chemical structure
  !! - molid: id of the molecule in cx
  !!
  !************************************************************************
  !
  Function MolecularMass(cx, molid)
    implicit none
    type(cxs) :: cx
    integer :: i, j, k, molid
    double precision :: MolecularMass

    MolecularMass = 0.0
    do i = 1, cx%namol(molid)
      MolecularMass = MolecularMass + cx%mass(cx%molid(molid, i))
    enddo
    return
  end Function MolecularMass


end module chemstr
