!
!***************************************************************************************
!
!> @brief Reaction-path object definition.
!
!> Defines the reaction-path object type. This object is
!! a string of nimage snapshots along a reaction-path.
!
!***************************************************************************************
!

Module rpath

  use globaldata
  use pes
  use chemstr
  use constants
  use functions
  implicit none

  !> Type definition for the reaction-path object.
  !
  type rxp
    sequence
    integer :: nimage                                   !< Number of images
    integer :: na                                       !< Number of atoms per image
    type(cxs), dimension(:), allocatable :: cx          !< Allocated chemical structure objects in the path.
    real(8), dimension(:,:,:), allocatable :: coeff     !< Fourier coefficients of path (3,na,nb)
    real(8), dimension(:,:,:), allocatable :: pcoeff    !< Momentum of Fourier coefficients of path (3,na,nb)
    real(8), dimension(:,:,:), allocatable :: dcoeff    !< Derivatives of Fourier coefficients of path (3,na,nb)
    real(8), dimension(:), allocatable :: ks            !< intra bead force constants
    real(8) :: vspring                                  !< Spring potential term.
  end type rxp

contains

  !
  !************************************************************************
  !> NewPath
  !!
  !! Initializes a path object.
  !!
  !! - rp: Initialized path object.
  !! - startfrompath: Flag to indicate that we're starting from a full path (logical)
  !! - startfile: filename for starting point structure (xyz file)
  !! - endfile: filename for end point structure (xyz file)
  !! - pathfile: If startfrompath == .TRUE., read the full path from pathfile
  !! - nimage: number of images along path in total (read from input file OR
  !!          from pathfile)
  !! - pathinit: path initialization type when startfrompath = .FALSE..
  !!            Allowed values are 'LINEAR' and ....?
  !! - emptypath: Flag indicating whether to start with an empty path.
  !! - natom: number of atoms.
  !!
  !************************************************************************
  !
  Subroutine NewPath(rp, startfrompath, startfile, endfile, pathfile, nimage, pathinit, &
       emptypath, natom)
    implicit none
    real(8), dimension(:), allocatable :: x, y, z
    type(rxp) :: rp
    integer :: nimage, natom, i
    character (len=*), intent(in) :: startfile, endfile, pathfile, pathinit
    logical, intent(in) :: startfrompath, emptypath
    character (len=2), dimension(:), allocatable :: label

    ! if emptypath = .TRUE., then we just initialize an 'empty' path object with enough
    ! space for nimage images and natom atoms.
    !

    if (emptypath) then
      rp%nimage = nimage
      allocate(rp%cx( rp%nimage ))
      rp%na = natom
      allocate(label(natom), x(natom), y(natom), z(natom))
      label(:) = 'C'
      x(:) = 0.0
      y(:) = 0.0
      z(:) = 0.0
      do i = 1, nimage
        call CreateCXS(rp%cx(i), rp%na, label, x, y, z)
      enddo
      deallocate(label, x, y, z)

    else
      ! Initialize path by either reading entire path from file....
      !
      if (startfrompath) then
        call ReadPathFromFile(rp, pathfile)

        ! ..or by reading end-points then using pathinit method to fill in internals..
        !
      else
        call StartFromEndPoints(rp, startfile, endfile, pathinit, nimage)
      endif
    endif

    ! Allocate space for Fourier coefficients.
    !
    if(.not. allocated(rp%coeff)) allocate(rp%coeff(3, rp%na, rp%nimage))
    if(.not. allocated(rp%pcoeff)) allocate(rp%pcoeff(3, rp%na, rp%nimage))
    if(.not. allocated(rp%dcoeff)) allocate(rp%dcoeff(3, rp%na, rp%nimage))
    if(.not. allocated(rp%ks)) allocate(rp%ks(nimage-1))

    return
  end Subroutine NewPath


  !
  !************************************************************************
  !> CopyPath
  !!
  !! Creates an empty reaction path, rpout, using parameter sizes from rpin
  !!
  !! - rpin: The reaction path object used to get size information.
  !! - rpout: The new reaction path object.
  !!
  !************************************************************************
  !
  Subroutine CopyPath(rpin,rpout)
    implicit none
    type(rxp) :: rpin, rpout
    character (len=5) :: cdum
    integer :: i

    call NewPath(rpout, .FALSE. , cdum, cdum, cdum, rpin%nimage, &
        cdum, .TRUE., rpin%na)

    do i = 1, rpin%nimage
      rpout%cx(i) = rpin%cx(i)
    enddo

    return
  end Subroutine CopyPath


  !
  !************************************************************************
  !> DeletePath
  !!
  !! Deletes memory associated with a path.
  !!
  !! - rp: The reaction path object to be deleted.
  !!
  !************************************************************************
  !
  Subroutine DeletePath(rp)
    implicit none
    type(rxp) :: rp
    integer :: i

    do i = 1, rp%nimage
      if (allocated(rp%cx(i)%r )         ) deallocate( rp%cx(i)%r )
      if (allocated(rp%cx(i)%p )         ) deallocate( rp%cx(i)%p )
      if (allocated(rp%cx(i)%dvdr )      ) deallocate( rp%cx(i)%dvdr )
      if (allocated(rp%cx(i)%force )     ) deallocate( rp%cx(i)%force )
      if (allocated(rp%cx(i)%atomlabel ) ) deallocate( rp%cx(i)%atomlabel )
      if (allocated(rp%cx(i)%mass )      ) deallocate( rp%cx(i)%mass )
      if (allocated(rp%cx(i)%fixeddof )  ) deallocate( rp%cx(i)%fixeddof )
      if (allocated(rp%cx(i)%fixedatom ) ) deallocate( rp%cx(i)%fixedatom )
      if (allocated(rp%cx(i)%graph)      ) deallocate( rp%cx(i)%graph)
      if (allocated(rp%cx(i)%molid )     ) deallocate( rp%cx(i)%molid )
      if (allocated(rp%cx(i)%namol )     ) deallocate( rp%cx(i)%namol )
      if (allocated(rp%cx(i)%molen )     ) deallocate( rp%cx(i)%molen )
    enddo
    deallocate( rp%cx )
    deallocate( rp%coeff )
    deallocate( rp%pcoeff )
    deallocate( rp%dcoeff )

    return
  end Subroutine DeletePath


  !
  !************************************************************************
  !> SetInternalCoords
  !!
  !! When building a path from end-points only, this subroutine sets the
  !! coordinates of the nimage-2 internal beads according to the
  !! method defined in pathinit.
  !!
  !! - rp: The reaction path object.
  !! - pathinit: string defining interpolation type. Allowed values
  !!             are 'LINEAR' for linear interpolation..
  !!
  !************************************************************************
  !
  Subroutine SetInternalCoords(rp, pathinit)
    Integer :: i, j, na
    type(rxp) :: rp
    character (len=*) :: pathinit
    character (len=2), dimension(NAMAX) :: label
    real(8) :: lambda
    real(8), dimension(NAMAX) :: x, y, z

    ! Linear interpolation of the path.
    !
    if (pathinit == 'linear') then

      do i = 2, rp%nimage-1
        lambda = dble(i-1)/dble(rp%nimage-1)
        do j = 1, rp%na
          x(j) = (1.0-lambda) * rp%cx(1)%r(1, j) + lambda * rp%cx(rp%nimage)%r(1, j)
          y(j) = (1.0-lambda) * rp%cx(1)%r(2, j) + lambda * rp%cx(rp%nimage)%r(2, j)
          z(j) = (1.0-lambda) * rp%cx(1)%r(3, j) + lambda * rp%cx(rp%nimage)%r(3, j)
          label(j) = rp%cx(1)%atomlabel(j)
        enddo
        call CreateCXS(rp%cx(i), rp%na, label, x, y, z)
      enddo

    else
      stop '* ERROR in SetInternalCoords in rpath.f90: unknown pathinit type'
    endif

    return
  end Subroutine SetInternalCoords


  !
  !************************************************************************
  !> PrintPathToFile
  !!
  !! Outputs an entire reaction path to an xyz file.
  !!
  !! - rp: The reaction path object.
  !! - filename: the output xyz file.
  !! - do_append: (optional) whether to append to existing or overwrite/write new.
  !!
  !************************************************************************
  !
  Subroutine PrintPathToFile(rp, filename, do_append)
    implicit none
    type(rxp), intent(in) :: rp
    character(len=*), intent(in) :: filename
    logical, intent(in), optional :: do_append
    real*8 :: x, y, z
    integer :: i, j
    logical :: append

    if (present(do_append)) then
      append = do_append
    else
      append = .false.
    endif

    if (append) then
      open(13, file=filename, status='old', access='append')
    else
      open(13, file=filename, status='unknown')
    endif
    do i = 1, rp%nimage
      write(13, '(i5)') rp%na
      write(13, *)
      do j = 1, rp%na
        x = rp%cx(i)%r(1, j) * bohr_to_ang
        y = rp%cx(i)%r(2, j) * bohr_to_ang
        z = rp%cx(i)%r(3, j) * bohr_to_ang
        write(13, '(a2,2x,3(f14.8,2x))') rp%cx(i)%atomlabel(j), x, y, z
      enddo
    enddo
    close(13)

    return
  end Subroutine PrintPathToFile

 
  !
  !************************************************************************
  !> FindIDPPPath
  !!
  !! finds an interpolated path using NEB on a Image-dependent PairPotential,
  !! as described in Smidstrup et al J Chem phys 2014, 140
  !!
  !! - rp: Reaction path
  !!
  !************************************************************************
  !
  subroutine FindIDPPPath(rp, IDPPIter, IDPPConv, IDPPStep, IDPPspring)
    type(rxp) :: rp
    integer :: n, i, iter, j, k, idof, IDPPIter, Iterout, nactmol
    integer, dimension(rp%cx(1)%nmol) :: MolecAct
    real(8) :: fnorm1, fnorm2, fmax, ks(rp%nimage-1), IDPPConv, IDPPStep, IDPPspring, rnd, ir
    real(8), allocatable :: rtarg(:,:,:), rend(:,:), rst(:,:), dx(:)
    real(8) :: lambda

    ! Assign NEB spring in this case.
    !
    rp%ks = IDPPspring

    ! Generate IDPP target distances by linear interpolation.
    !
    allocate(rtarg(rp%nimage, rp%na, rp%na), rst(rp%na,rp%na), rend(rp%na, rp%na), dx(3))
    
    n = rp%nimage
    do i = 1, rp%na
      do j = 1, rp%na
        dx(1:3) = rp%cx(1)%r(1:3, i) - rp%cx(1)%r(1:3, j)
        rst(i, j) = sqrt(dx(1)*dx(1) + dx(2) * dx(2) + dx(3)*dx(3))
        dx(1:3) = rp%cx(n)%r(1:3, i) - rp%cx(n)%r(1:3, j)
        rend(i, j) = sqrt(dx(1)*dx(1) + dx(2) * dx(2) + dx(3)*dx(3))
      enddo
    enddo
    do i = 2, rp%nimage-1
      do j = 1, rp%na
        do k = 1, rp%na
          lambda = dble(i-1) / dble(rp%nimage)
          rtarg(i, j, k) = rst(j, k) + lambda * (rend(j ,k) - rst(j, k))
        enddo
      enddo
    enddo
	
    ! Calculate IDPP forces on all internal images.
    !
    call GetIDPPForces(rp, rtarg)
       
    ! Calculate projected forces.
    !
    select case (projforcetype)
      case(1)
       call GetProjForces1(rp, .true., .FALSE., .FALSE.)
      case(2)
       call GetProjForces2(rp, .true., .FALSE., .FALSE.)
      case(3)
       call GetProjForces3(rp, .true., .FALSE., .FALSE.)
    end select

    ! Loop over SD iterations.
    !
    do iter = 1, IDPPIter
      ! Update positions.
      !
      do i = 2, rp%nimage-1
        idof = 0
        do j = 1, rp%na
          if (.not. rp%cx(i)%fixedatom(j)) then
            do k = 1, 3
              idof = idof + 1
              if (.not. rp%cx(i)%FixedDOF(idof)) then
                rp%cx(i)%r(k, j) = rp%cx(i)%r(k, j) + IDPPstep * rp%cx(i)%force(k, j)
              endif
            enddo
          else
            idof = idof + 3
          endif
        enddo
      enddo
       
      ! Calculate new forces.
      !
      call GetIDPPForces(rp, rtarg)
       
      ! Project forces.
      !
      select case (projforcetype)
        case(1)
          call GetProjForces1(rp, .false., .FALSE., optendsduring)
        case(2)
          call GetProjForces2(rp, .false., .FALSE., optendsduring)
        case(3)
          call GetProjForces3(rp, .false., .FALSE., optendsduring)
        case default
      end select

      call GetForceNorm(rp, fnorm1, fmax, 2, rp%nimage-1)

		  ! Check convergence.
      if (fnorm1 <= IDPPConv) then
        write(logfile, '("* |Forces| < NEBconv :: IDPP NEB CONVERGED, Iter = ",I5)') iter
        exit
      endif
    enddo

    return
  end subroutine FindIDPPPath


  !
  !************************************************************************
  !
  ! GetIDPPForces
  !
  ! Calculates the IDPP forces.
  !
  !************************************************************************
  !
  Subroutine GetIDPPForces(rp, rtarg)
    implicit none
    integer :: i, j, k
    type(rxp) :: rp
    real(8), dimension(rp%nimage, rp%na, rp%na) :: rtarg
    real(8), dimension(3) :: d, dxr
    real(8) :: dr, onr, onr4, dEdr, IDPPenergy, t1
    
    IDPPenergy = 0.d0
    do k = 2, rp%nimage - 1
      rp%cx(k)%force = 0.d0
      do i = 1, rp%na-1
        do j = i+1, rp%na
          d(1:3) = rp%cx(k)%r(1:3, i) - rp%cx(k)%r(1:3, j)
          dr = sqrt(d(1)*d(1) + d(2)*d(2) + d(3)*d(3))
          onr = 1.0 / dr
          onr4 = onr**4
          dxr(1:3) = d(1:3) * onr
          t1 = dr - rtarg(k, i, j)
          IDPPenergy = IDPPenergy + onr4 * t1 * t1
          dEdr = 2.0 * t1 * onr4 - 4.0 * t1 * t1 * (onr4*onr)
          rp%cx(k)%force(1:3, i) = rp%cx(k)%force(1:3, i) - dEdr * dxr(1:3)
          rp%cx(k)%force(1:3, j) = rp%cx(k)%force(1:3, j) + dEdr * dxr(1:3)
        enddo
      enddo
      rp%cx(k)%dvdr(:, :) = -rp%cx(k)%force(:, :)
    enddo

    return
  end Subroutine GetIDPPForces
  

  !
  !************************************************************************
  !
  !> GetForceNorm
  !!
  !! Calculate the norm of the projected forces for the internal beads.
  !!
  !! - rp: The reaction-path object to be refined.
  !! - forcenorm: Returned force norm
  !! - n1, n2: Integer range of images to include in calculation.
  !!
  !************************************************************************
  !
  Subroutine GetForceNorm(rp, forcenorm, fmax, n1, n2)
    implicit none
    type(rxp) :: rp
    real(8) :: forcenorm,fmax
    integer :: i, j, k, isum, idof, n1, n2

    forcenorm = 0.0
    fmax = -1d0
    isum = 0
    do i = n1, n2
      idof = 0
      do j = 1, rp%na
        if (.not. rp%cx(i)%Fixedatom(j)) then
          do k = 1, 3
            idof = idof + 1
            if (.not. rp%cx(i)%FixedDof(idof)) then
              isum = isum + 1
              forcenorm = forcenorm + rp%cx(i)%force(k, j)**2
              if (abs(rp%cx(i)%force(k, j) ) > fmax) then
                fmax = abs(rp%cx(i)%force(k, j))
              endif
            endif
          enddo
        else
          idof = idof + 3
        endif
      enddo
    enddo

    if (dble(isum) .gt. 0.0d0) then
      forcenorm = sqrt(forcenorm / dble(isum))
    else
      forcenorm = 0.0d0
    endif

    return
  end Subroutine GetForceNorm


  !
  !************************************************************************
  !> ReadXYZFrame
  !!
  !! reads a single frame from an xyz file.
  !!
  !! - fileid: Integer file label.
  !! - na: Number of atoms
  !! - label: atom labes read from file.
  !! - x,y,z: Coordinates of each atom read from file
  !!
  !************************************************************************
  !
  Subroutine ReadXYZFrame(fileid, na, label, x, y, z)
    integer :: na, fileid, i
    character (len=2), dimension(:) :: label
    character (len=100) :: comment
    real(8), dimension(:) :: x, y, z

    read(fileid, *) na
    read(fileid, '(A100)') comment

    do i = 1, na
      read(fileid, *) label(i), x(i), y(i), z(i)
      x(i) = x(i) * ang_to_bohr
      y(i) = y(i) * ang_to_bohr
      z(i) = z(i) * ang_to_bohr
    enddo
    return
  end Subroutine ReadXYZFrame


  !
  !************************************************************************
  !> StartFromEndPoints
  !!
  !! Initializes a reaction path using the path end-points.
  !!
  !! - rp: The reaction path object.
  !! - startfile: the xyz file containing the startpoint.
  !! - endfile: the xyz file containing the endpoint.
  !! - pathinit: character*5 string containing identifier for method used
  !!            to interpolate internal image coordinates.
  !! - nimage: Number of images in path.
  !!
  !************************************************************************
  !
  Subroutine StartFromEndPoints(rp, startfile, endfile, pathinit, nimage)
    type(rxp) :: rp
    integer :: ios, na, nastart, i ,j, nimage
    character (len=*) :: pathinit, startfile, endfile
    logical :: there
    real(8), dimension(NAMAX) :: x, y, z
    real(8) :: lambda
    character (len=2), dimension(NAMAX) :: label
    character (len=100) :: comment

    ! Allocate workspace for chemical structures.
    !
    rp%nimage = nimage
    allocate(rp%cx(rp%nimage))

    ! Initialize startpoint from file.
    !
    inquire(file = startfile, exist = there)
    if (.not. there) stop '* ERROR in StartFromEndPoints in rpath.f90: specified startfile does not exist'
    open(12, file = startfile, status = 'unknown')
    call ReadXYZFrame(12, na, label, x, y, z)
    call CreateCXS(rp%cx(1), na, label, x, y, z)
    close(12)

    ! Store number of atoms to check end-point compatibility....
    !
    nastart = na
    rp%na = na

    ! Initialize endpoint from file.
    !
    inquire( file = endfile, exist = there )
    if (.not. there) stop '* ERROR in StartFromEndPoints in rpath.f90: specified endfile does not exist'
    open(12, file = endfile, status = 'unknown')
    call ReadXYZFrame(12, na, label, x, y, z)
    call CreateCXS(rp%cx(nimage), na, label, x, y, z)
    close(12)

    ! Check that numbers of atoms make sense.
    !
    if (na /= nastart) then
      stop '* ERROR in StartFromEndPoints in rpath.f90: Incompatible atom numbers in start- and end-points'
    endif

    ! Initialize nimage-2 internal beads using method defined in pathinit.
    !
    call SetInternalCoords(rp, pathinit)

    return
  end subroutine StartFromEndPoints

  !
  !************************************************************************
  !> ReadPathFromFile
  !!
  !! Reads a file reaction path from a xyz file with nimage frames.
  !! First, we read the number of atoms and images, then create the path object.
  !! Finally, we read in the coordinates along the path.
  !!
  !! - rp: Reaction path object.
  !! - pathfile: filename for path.
  !!
  !************************************************************************
  !
  Subroutine ReadPathFromFile(rp, pathfile)
    type(rxp) :: rp
    integer :: i, na, ios, nimage
    character (len=*) :: pathfile
    logical :: there
    real(8), dimension(NAMAX) :: x, y, z
    character (len=2), dimension(NAMAX) :: label

    ! Check pathfile existence.
    !
    inquire( file = pathfile, exist = there )
    if (.not. there) stop '* ERROR in ReadPathFromFile in rpath.f90: specified input file does not exist'

    ! Open file.
    open(12, file = pathfile, status = 'unknown')

    ! First pass through...read number of atoms.
    !
    ios = 0
    nimage = 0
    do while (ios == 0)
      read(12, *, iostat=ios) na

      if (ios /= 0)exit

      read(12,*,iostat=ios)
      do i = 1, na
        read(12, *, iostat=ios) label(i), x(i), y(i), z(i)
      enddo
      nimage = nimage + 1
    enddo

    ! Close file.
    !
    close(12)

    ! Set number of images.
    !
    rp%nimage = nimage
    rp%na = na

    ! Allocate chemical structures.
    !
    allocate(rp%cx(rp%nimage))

    ! Open pathfile again.
    !
    open(12, file = pathfile, status = 'unknown')

    ! First pass through...read number of atoms.
    !
    ios = 0
    nimage = 0
    do while (ios == 0)
      read(12, *, iostat=ios) na

      if (ios /= 0) exit

      read(12, *, iostat=ios)
      if (ios /= 0) exit
      do i = 1, na
        read(12, *, iostat=ios) label(i), x(i), y(i), z(i)
        x(i) = x(i) * ang_to_bohr
        y(i) = y(i) * ang_to_bohr
        z(i) = z(i) * ang_to_bohr
      enddo
      nimage = nimage + 1

      ! Create chemical structure.
      !
      call CreateCXS(rp%cx(nimage), na, label, x, y, z)
    enddo
    close(12)
    ! Allocate space for Fourier coefficients.
    !
    allocate(rp%coeff(3, rp%na, rp%nimage))
    allocate(rp%pcoeff(3, rp%na, rp%nimage))
    allocate(rp%dcoeff(3, rp%na, rp%nimage))

  end Subroutine ReadPathFromFile


  !
  !************************************************************************
  !> SetPathConstraints
  !!
  !! Identify constrained atoms and DOFs in the path.
  !!
  !! - rp: Initialized path object.
  !! - NDOFconstr: Number of DOF constraints.
  !! - FixedDOF: Array containing integer number of constrained DOFs
  !! - Natomconstr: Number of atom constraints
  !! - Fixedatom: Array containing integer ids of fixed atoms.
  !!
  !************************************************************************
  !
  Subroutine SetPathConstraints(rp, NDOFconstr, FixedDOF, Natomconstr, Fixedatom)
    implicit none
    type(rxp) :: rp
    integer, intent(in) :: NDOFconstr, Natomconstr
    integer, intent(in), dimension(:) :: FixedDOF, FixedAtom
    integer :: i

    do i = 1, rp%nimage
      call SetCXSconstraints(rp%cx(i), NDOFconstr, FixedDOF, Natomconstr, FixedAtom)
    enddo

    return
  end Subroutine SetPathConstraints


  !
  !************************************************************************
  !> GetPathEnergy
  !!
  !! Calculates the energy for all structures along a reaction path.
  !! The pes.f90 module contains details of the calculation to perform.
  !!
  !! - rp: Path object.
  !!
  !************************************************************************
  !
  Subroutine GetPathEnergy(rp, success)
    implicit none
    type(rxp) :: rp
    integer :: i
    logical :: minimize, success

    ! Loop over images, calculating energy for each.
    !
    minimize = .false.
    do i = 1, rp%nimage
      if (.not. pesfull) then
        call GetMols(rp%cx(i))
      endif
      call AbInitio(rp%cx(i), 'ener', success)
    enddo

    return
  end Subroutine GetPathEnergy


  !
  !************************************************************************
  !> GetPathGradients
  !!
  !! Calculates the energy for all structures along a reaction path.
  !! The pes.f90 module contains details of the calculation to perform.
  !!
  !! - rp: Path object.
  !! - success: Result of Ab Initio calculation.
  !! - calc_all: Whether to force calculation of all cxs along rp.
  !!
  !************************************************************************
  !
  Subroutine GetPathGradients(rp, success, calc_all)
    implicit none
    type(rxp), intent(in) :: rp
    logical, intent(out) :: success
    logical, intent(in) :: calc_all
    integer :: i

    ! Loop over images, calculating energy for each.
    !
    do i = 1, rp%nimage
      ! Don't recalculate endpoints if they don't change.
      if ((.not. calc_all) .and. (i == 1 .or. i == rp%nimage) .and. &
          (.not. optendsbefore .or. .not. optendsduring)) cycle
      if (.not. pesfull) then
        call GetMols(rp%cx(i))
      endif
      call AbInitio(rp%cx(i), 'grad', success)
    enddo

    return
  end Subroutine GetPathGradients


  !
  !*************************************************************************
  !
  !> GetProjForces1
  !!
  !! Calculates the forces projected along the reaction path, as required in
  !! NEB or CINEB calculations. A previous PES calculation MUST have
  !! been performed for the reaction path.
  !!
  !! Note: By setting ci_flag = .true., this routine returns forces for
  !! CI-NEB.
  !!
  !! - rp: Reaction path object.
  !! - kon: Flag determining if the spring forces should be calculated or not
  !! - ci_flag: Logical flag indicating whether to return climbing-image
  !!            forces.
  !! - endflag: Logical flag indicatin whether or not the forces on
  !!            the end-point beads should be returned as the
  !!            "bare" (non-NEB) forces.
  !!
  !
  !*************************************************************************
  !
  Subroutine GetProjForces1(rp, kon, ci_flag, endflag)
    implicit none
    type (rxp) :: rp
    integer :: f, i, j, k, imax, natom, idof
    real(8) :: sum, kspring, dot, beta, betan, omega, delk, Ei, Emax, Eref, kmax
    real(8) :: x0, xm1, xp1, vmax, kvariability, kk(rp%nimage-1)
    double precision :: tp(rp%na*3), tm(rp%na*3), tt(rp%na*3), vp, vm, vv, vmx, vmn
    real(8), allocatable :: tangent(:,:,:), vec1(:,:), vec2(:,:)
    logical :: flag, ci_flag, endflag, kon

    Emax = -1d6
    do i = 1, rp%nimage
      if (rp%cx(i)%vcalc > Emax ) then
        Emax = rp%cx(i)%vcalc
        imax = i
      endif
    enddo
    if (kon) then
      kk(1:rp%nimage-1) = rp%ks(1:rp%nimage-1)
    else
      kk = 0.0
    endif

    ! Allocate space for tangents and vectors.
    !
    natom = rp%cx(1)%na
    allocate(tangent(3, natom, rp%nimage), vec1(3,natom), vec2(3,natom))
    tangent(:, :, :) = 0.d0
    vec1(:, :) = 0.d0
    vec2(:, :) = 0.d0

    ! Identify the intermediate image with the largest energy if we're running CI-NEB.
    !
    if (ci_flag) then
      vmax = -1d6
      do i = 2, rp%nimage-1
        if (rp%cx(i)%vcalc > vmax) then
          imax = i
          vmax = rp%cx(i)%vcalc
        endif
      enddo
    endif

    ! Calculate tangents.
    !
    do i = 1, rp%nimage

      if (i == 1) then
        sum = 0.0
        idof = 0
        do j = 1, natom
          if (.not. rp%cx(i)%fixedatom(j)) then
            do k = 1, 3
              idof = idof + 1
              if (.not. rp%cx(i)%fixeddof(idof)) then
                tangent(k, j, i) = rp%cx(2)%r(k, j) - rp%cx(1)%r(k, j)
                sum = sum + tangent(k, j, i) * tangent(k, j, i)
              endif
            enddo
          else
            idof = idof + 3
          endif
        enddo
        tangent(:, :, i) = tangent(:, :, i) / sqrt(sum)

      else if (i == rp%nimage) then

        sum = 0.0
        idof = 0
        do j = 1, natom
          if (.not. rp%cx(i)%fixedatom(j)) then
            do k = 1, 3
              idof = idof + 1
              if (.not. rp%cx(i)%fixeddof(idof)) then
                tangent(k, j, i) = rp%cx(rp%nimage)%r(k, j) - rp%cx(rp%nimage-1)%r(k, j)
                sum = sum + tangent(k, j, i) * tangent(k, j, i)
              endif
            enddo
          else
            idof = idof + 3
          endif
        enddo
        tangent(:, :, i) = tangent(:, :, i) / sqrt(sum)

      else  ! More accurate bisector vector !
        sum = 0.0
        idof = 0
        do j = 1, natom
          if (.not. rp%cx(i)%fixedatom(j)) then
            do k = 1, 3
              idof = idof + 1
              if (.not. rp%cx(i)%fixeddof(idof)) then
                vec1(k, j) = rp%cx(i+1)%r(k, j) - rp%cx(i)%r(k, j)
                sum = sum + vec1(k, j)*vec1(k, j)
              endif
            enddo
          else
            idof = idof + 3
          endif
        enddo
        vec1(:, :) = vec1(:, :) / sqrt(sum)

        sum = 0.0
        idof = 0
        do j = 1, natom
          if (.not. rp%cx(i)%fixedatom(j)) then
            do k = 1, 3
              idof = idof + 1
              if (.not.rp%cx(i)%fixeddof(idof)) then
                vec2(k, j) = rp%cx(i)%r(k, j) - rp%cx(i-1)%r(k, j)
                sum = sum + vec2(k, j)*vec2(k, j)
              endif
            enddo
          else
            idof = idof + 3
          endif
        enddo
        vec2(:, :) = vec2(:, :) / sqrt(sum)

        tangent(:, :, i) = vec1(:, :) + vec2(:, :)

        sum = 0.0
        idof = 0
        do j = 1, natom
          if (.not. rp%cx(i)%fixedatom(j)) then
            do k = 1, 3
              idof = idof + 1
              if (.not. rp%cx(i)%fixeddof(idof)) then
                sum = sum + tangent(k, j, i) * tangent(k, j, i)
              endif
            enddo
          else
            idof = idof + 3
          endif
        enddo
        tangent(:, :, i) = tangent(:, :, i) / sqrt(sum)

      endif

    enddo

    ! Loop over chain-states and calculate projected forces.
    !
    do i = 1, rp%nimage
      rp%cx(i)%force(:, :) = 0.0
    enddo

    do i = 1, rp%nimage

      ! Now calculate the spring forces on image i.
      !
      if (i > 1 .and. i < rp%nimage) then

        dot = 0.d0
        idof = 0
        do j = 1, natom
          if (.not. rp%cx(i)%fixedatom(j)) then
            do k = 1, 3
              idof = idof + 1
              if (.not. rp%cx(i)%fixeddof(idof)) then
                xp1 = rp%cx(i+1)%r(k, j)
                x0 = rp%cx(i)%r(k, j)
                xm1 = rp%cx(i-1)%r(k, j)

                ! dot = dot + kspring * ( ( xp1 - x0 ) - (x0 - xm1) ) * tangent(k,j,i)
                dot = dot + kk(i-1) * ( ( xp1 - x0 ) - (x0 - xm1) ) * tangent(k, j, i)
                ! rp%cx(i)%force(k,j) =   rp%cx(i)%force(k,j) + kk(i-1) * ( ( xp1 - x0 ) - (x0 - xm1) ) * tangent(k,j,i)

              endif
            enddo
          else
            idof = idof + 3
          endif
        enddo


      else if (i == 1) then

        dot = 0.d0
        idof = 0
        do j = 1, natom
          if (.not. rp%cx(i)%fixedatom(j)) then
            do k = 1, 3
              idof = idof + 1
              if (.not. rp%cx(i)%fixeddof(idof)) then
                xp1 = rp%cx(i+1)%r(k, j)
                x0 = rp%cx(i)%r(k, j)
                ! dot = dot + kspring* (  xp1 - x0 )* tangent(k,j,i)
                dot = dot + kk(i) * (xp1 - x0 )* tangent(k, j, i)
                ! rp%cx(i)%force(k,j) = rp%cx(i)%force(k,j) + kk(i) * (  xp1 - x0 )* tangent(k,j,i)
              endif
            enddo
          else
            idof = idof + 3
          endif
        enddo

      else if (i == rp%nimage) then

        dot = 0.d0
        idof = 0
        do j = 1, natom
          if (.not. rp%cx(i)%fixedatom(j)) then
            do k = 1, 3
              idof = idof + 1
              if (.not. rp%cx(i)%fixeddof(idof)) then
                x0 = rp%cx(i)%r(k, j)
                xm1 = rp%cx(i-1)%r(k, j)
                ! dot = dot - kspring*( x0-xm1 )* tangent(k,j,i)
                dot = dot - kk(i-1)*( x0-xm1 )* tangent(k, j, i)
                ! rp%cx(i)%force(k,j)=rp%cx(i)%force(k,j) - kk(i-1)*( x0-xm1 )* tangent(k,j,i)
              endif
            enddo
          else
            idof = idof + 3
          endif
        enddo

      endif

      ! Add the parallel spring contribution to the total forces.
      !
      idof = 0
      do j = 1, natom
        if (.not. rp%cx(i)%fixedatom(j)) then
          do k = 1, 3
            idof = idof + 1
            if (.not. rp%cx(i)%fixeddof(idof)) then
              rp%cx(i)%force(k, j) = dot * tangent(k, j, i)
            endif
          enddo
        else
          idof = idof + 3
        endif
      enddo

      ! Add the contribution of the perpendicular forces.
      !
      idof = 0
      dot = 0.d0
      do j = 1, natom
        if (.not. rp%cx(i)%fixedatom(j)) then
          do k = 1, 3
            idof = idof + 1
            if (.not. rp%cx(i)%fixeddof(idof)) then
              dot = dot + rp%cx(i)%dvdr(k, j) * tangent(k, j, i)
            endif
          enddo
        else
          idof = idof + 3
        endif
      enddo

      idof = 0
      do j = 1, natom
        if (.not. rp%cx(i)%fixedatom(j)) then
          do k = 1, 3
            idof = idof + 1
            if (.not. rp%cx(i)%fixeddof(idof)) then
              rp%cx(i)%force(k, j) = rp%cx(i)%force(k, j)-rp%cx(i)%dvdr(k, j) + dot * tangent(k, j, i)
            endif
          enddo
        else
          idof = idof + 3
        endif
      enddo

      ! Sort out CI-NEB force for imax (highest energy image).
      !
      if (ci_flag .and. i == imax) then

        idof = 0
        do j = 1, natom
          if (.not. rp%cx(i)%fixedatom(j)) then
            do k = 1, 3
              idof = idof + 1
              if (.not. rp%cx(i)%fixeddof(idof)) then
                rp%cx(i)%force(k, j) = -rp%cx(i)%dvdr(k, j)
              endif
            enddo
          else
            idof = idof + 3
          endif
        enddo

        dot = 0.d0
        idof = 0
        do j = 1, natom
          if (.not. rp%cx(i)%fixedatom(j)) then
            do k = 1, 3
              idof = idof + 1
              if (.not. rp%cx(i)%fixeddof(idof)) then
                dot = dot + rp%cx(i)%dvdr(k, j) * tangent(k, j, i)
              endif
            enddo
          else
            idof = idof + 3
          endif
        enddo

        idof = 0
        do j = 1, natom
          if (.not. rp%cx(i)%fixedatom(j)) then
            do k = 1, 3
              idof = idof + 1
              if (.not. rp%cx(i)%fixeddof(idof)) then
                rp%cx(i)%force(k, j) = rp%cx(i)%force(k, j)+ 2.d0 * dot * tangent(k, j, i)
              endif
            enddo
          else
            idof = idof + 3
          endif
        enddo
      endif
    enddo
    deallocate(tangent, vec1, vec2)

    ! If endflag is .TRUE., we just need to return the "bare" forces on the end-points.
    !
    if (endflag) then

      idof = 0
      do j = 1, natom
        if (.not. rp%cx(1)%fixedatom(j)) then
          do k = 1, 3
            idof = idof + 1
            if (.not. rp%cx(1)%fixeddof(idof)) then
              rp%cx(1)%force(k, j) = -rp%cx(1)%dvdr(k, j)
            endif
          enddo
        else
          idof = idof + 3
        endif
      enddo

      idof = 0
      do j = 1, natom
        if (.not. rp%cx(rp%nimage)%fixedatom(j)) then
          do k = 1, 3
            idof = idof + 1
            if (.not. rp%cx(rp%nimage)%fixeddof(idof)) then
              rp%cx(rp%nimage)%force(k, j) = -rp%cx(rp%nimage)%dvdr(k, j)
            endif
          enddo
        else
          idof = idof + 3
        endif
      enddo

    endif

    return
  end Subroutine GetProjForces1


  !
  !*************************************************************************
  !
  !> GetProjForces2
  !!
  !! Calculates the forces projected along the reaction path, as required in
  !! NEB or CINEB calculations. A previous PES calculation MUST have
  !! been performed for the reaction path. This routine uses a slightly
  !! improved version of tangents and forces as described by henkelman 00
  !!
  !! Note: By setting ci_flag = .true., this routine returns forces for
  !! CI-NEB.
  !!
  !! - rp: Reaction path object.
  !!   size is made.
  !! - ci_flag: Logical flag indicating whether to return climbing-image
  !!            forces.
  !! - endflag: Logical flag indicating whether or not the forces on
  !!            the end-point beads should be returned as the
  !!            "bare" (non-NEB) forces.
  !!
  !
  !*************************************************************************
  !
  Subroutine GetProjForces2(rp, kon, ci_flag, endflag)
    implicit none
    integer :: i, j, k, l, imax, idof, na, nd
    type (rxp) :: rp
    double precision :: kspring, Ei, kmax, delk, Emax, Eref, kvariability, kk(rp%nimage-1)
    double precision :: tp(rp%na*3), tm(rp%na*3), tt(rp%na*3), vp, vm, vv, vmx, vmn
    logical :: ci_flag, endflag, kon

    na = rp%na ; nd = rp%na*3
    Emax = -1d6
    do i = 1, rp%nimage
      if (rp%cx(i)%vcalc > Emax) then
        Emax = rp%cx(i)%vcalc
        imax = i
      endif
    enddo
    if (kon) then
      kk(1:rp%nimage-1) = rp%ks(1:rp%nimage-1)
    else
      kk = 0.0
    endif

    ! Calculate tangents and forces.
    !
    do i = 2, rp%nimage-1
      vv = rp%cx(i)%vcalc
      tp = reshape(rp%cx(i+1)%r - rp%cx(i)%r, (/nd/)) + epsil
      vp = rp%cx(i+1)%vcalc
      tm = reshape(rp%cx(i)%r - rp%cx(i-1)%r, (/nd/)) + epsil
      vm = rp%cx(i-1)%vcalc
      if(vp > vv .and. vv > vm) then
        tt = tp
      elseif(vp < vv .and. vv < vm)  then
        tt = tm
      elseif((vp > vv .and. vv < vm) .or. (vp < vv .and. vv > vm)) then
        vmx = max(abs(vp-vv), abs(vm-vv))
        vmn = min(abs(vp-vv), abs(vm-vv))
        if(vp > vm) then
          tt = tp*vmx+tm*vmn + epsil
        else
          tt = tp*vmn+tm*vmx + epsil
        endif
      else
        tt = tp/norm2(tp) + tm/norm2(tm) + epsil
      endif
                         ! if kk = 0 then make tt = 0
      tt = tt/norm2(tt) * (kk(i))/(kk(i)+epsil)
      if (ci_flag) then
        if(i /= imax) then
          ! Parallel forces to the band
          !
          rp%cx(i)%force(:, :) = reshape((kk(i)*norm2(tp)-kk(i-1)*norm2(tm))*tt, (/3, na/))
          ! Perpendicular forces of the system
          !
          rp%cx(i)%force(:, :) = rp%cx(i)%force(:, :) - rp%cx(i)%dvdr + &
              reshape(dot_product(reshape(rp%cx(i)%dvdr, (/nd/)), tt)*tt, (/3, na/))
        else
          ! climbing image force
          rp%cx(i)%force(:, :) = -rp%cx(i)%dvdr + 2.0d0 * &
              reshape(dot_product(reshape(rp%cx(i)%dvdr, (/nd/)), tt)*tt, (/3, na/))
        endif
      else
        ! Parallel forces to the band
        !
        rp%cx(i)%force(:, :) = reshape((kk(i)*norm2(tp)-kk(i-1)*norm2(tm))*tt, (/3, na/))
        ! Perpendicular forces of the system
        !
        rp%cx(i)%force(:, :) = rp%cx(i)%force(:, :) - rp%cx(i)%dvdr + &
            reshape(dot_product(reshape(rp%cx(i)%dvdr, (/nd/)), tt)*tt, (/3, na/))
      endif

    enddo

    ! first and last bead
    !
    if (endflag) then
      rp%cx(1)%force(:, :) = -rp%cx(1)%dvdr(:, :)
      rp%cx(rp%nimage)%force(:, :) = -rp%cx(rp%nimage)%dvdr(:, :)
    elseif(ci_flag) then
      if( 1 == imax) then
        tp = reshape(rp%cx(2)%r - rp%cx(1)%r, (/nd/))
        tp = tp/norm2(tp)
        rp%cx(1)%force(:, :) = -rp%cx(1)%dvdr + 2.0d0 * &
            reshape(dot_product(reshape(rp%cx(1)%dvdr, (/nd/)), tp)*tp, (/3, na/))
      elseif (rp%nimage == imax) then
        tm = reshape(rp%cx(rp%nimage)%r - rp%cx(rp%nimage-1)%r, (/nd/))
        tm = tm/norm2(tm)
        rp%cx(rp%nimage)%force(:, :) = -rp%cx(rp%nimage)%dvdr + 2.0d0 * &
            reshape(dot_product(reshape(rp%cx(rp%nimage)%dvdr, (/nd/)), tm)*tm, (/3, na/))
      endif
    endif

    ! set constraints...
    !
    do i = 1, rp%nimage
      idof = 0
      do j = 1, rp%na
        if (rp%cx(i)%fixedatom(j)) then
          rp%cx(i)%force(1:3, j) = 0.0d0
          idof = idof + 3
        else
          do k = 1, 3
            idof = idof + 1
            if (rp%cx(i)%fixeddof(idof)) then
                rp%cx(i)%force(k, j) = 0.0d0
            endif
          enddo
        endif
      enddo
    enddo

    !    print*,'WTF K: ',kk(1:rp%nimage-1)
    !    print*,'FORCES = ',rp%cx(5)%force

    return
  end Subroutine GetProjForces2


  Subroutine VariableSprings(rp)
    implicit none
    integer :: i, j, k, l, imax, idof, na, nd
    type (rxp) :: rp
    double precision :: kspring, Ei, kmax, delk, Emax, Eref, kvariability

    na = rp%na ; nd = rp%na*3
    Emax = -1d6
    do i = 1, rp%nimage
      if (rp%cx(i)%vcalc > Emax) then
        Emax = rp%cx(i)%vcalc
        imax = i
      endif
    enddo
    ! this is the extra strength or looseness wrt kspring that is made
    ! variable...:
    kvariability = 0.150d0
    ! determine spring strength
    !
    kmax = kspring + kspring * kvariability
    delk = 2.0d0 * kspring * kvariability
    Eref = max(rp%cx(1)%vcalc, rp%cx(rp%nimage)%vcalc)
    do i = 1, rp%nimage-1
      Ei = max(rp%cx(i)%vcalc, rp%cx(i+1)%vcalc)
      if (Ei > Eref) then
        rp%ks(i) = kmax - delk * (Emax - Ei)/(Emax - Eref)
      else
        rp%ks(i) = kmax - delk
      endif
    enddo

    return
  end Subroutine VariableSprings


  !
  !*************************************************************************
  !
  !> GetProjForces3
  !!
  !! Calculates the forces projected along the reaction path, as required in
  !! NEB or CINEB calculations.
  !! This routine uses improves the improvement in henkelmans two 00 papers
  !! as described in Kolsbjerg and Hammer 16
  !!
  !! Note: By setting ci_flag = .true., this routine returns forces for
  !! CI-NEB.
  !!
  !! - rp: Reaction path object.
  !!   size is made.
  !! - ci_flag: Logical flag indicating whether to return climbing-image
  !!            forces.
  !! - endflag: Logical flag indicatin whether or not the forces on
  !!            the end-point beads should be returned as the
  !!            "bare" (non-NEB) forces.
  !!
  !
  !*************************************************************************
  !
  Subroutine GetProjForces3(rp,kon,ci_flag,endflag)
    implicit none
    integer :: i, j, k, l, imax, idof, na, nd
    type (rxp) :: rp
    double precision :: kspring, Ei, kmax, delk, Emax, Eref, kvariability, kk(rp%nimage-1)
    double precision :: tp(rp%na*3), tm(rp%na*3), tt(rp%na*3), vp, vm, vv, vmx, vmn, Leq, tpn, tmn
    logical :: ci_flag, endflag, variablespring, kon

    na = rp%na ; nd = rp%na*3
    Emax = -1d6
    do i = 1, rp%nimage
      if (rp%cx(i)%vcalc > Emax ) then
        Emax = rp%cx(i)%vcalc
        imax = i
      endif
    enddo
    if (kon) then
      kk(1:rp%nimage-1) = rp%ks(1:rp%nimage-1)
    else
      kk = 0.0
    endif
    !do i = 2, rp%nimage
    ! Leq = Leq +norm2(reshape(rp%cx(i  )%r - rp%cx(i-1)%r,(/nd/)))
    !enddo
    !Leq = Leq/dble(rp%nimage-1)
    Leq = norm2(reshape(rp%cx(rp%nimage)%r, (/nd/)) - reshape(rp%cx(1)%r, (/nd/))) / dble(rp%nimage-1)
    ! Calculate tangents and forces.
    !
    do i = 2, rp%nimage-1
      rp%cx(i)%force(:, :) = 0.0d0
      vv = rp%cx(i)%vcalc
      tp = reshape(rp%cx(i+1)%r - rp%cx(i)%r, (/nd/))
      vp = rp%cx(i+1)%vcalc
      tm = reshape(rp%cx(i)%r - rp%cx(i-1)%r, (/nd/))
      vm = rp%cx(i-1)%vcalc
      tpn = norm2(tp) + epsil
      tmn = norm2(tm) + epsil
      if (ci_flag .and. i /= imax .or. .not. ci_flag) then
       ! Parallel forces to the band

        if (ci_flag .and. abs(i-imax) == 1) then
          vmx = max(abs(vp-vv), abs(vm-vv))
          vmn = min(abs(vp-vv), abs(vm-vv))
          rp%cx(i)%force(:, :) = reshape(((kk(i)*(tpn-Leq)*tp/tpn-kk(i-1)*(tmn-Leq)*tm/tmn))*(vmn/vmx), (/3, na/))
        else
          rp%cx(i)%force(:, :) = reshape((kk(i)*(tpn-Leq)*tp/tpn-kk(i-1)*(tmn-Leq)*tm/tmn), (/3, na/))
        endif
        !print*, 'FORCE NOW? = ', norm2(rp%cx(i)%force(:,:))

        tt = tp/tpn + tm/tmn + epsil
                              ! if kk = 0 then make tt = 0
        tt = (tt/norm2(tt))*(kk(i))/(kk(i)+epsil)
        ! Perpendicular forces of the system
        !
        !rp%cx(i)%force(:,:) = rp%cx(i)%force(:,:) - rp%cx(i)%dvdr
        rp%cx(i)%force(:, :) = rp%cx(i)%force(:, :) - rp%cx(i)%dvdr + &
            reshape(dot_product(reshape(rp%cx(i)%dvdr, (/nd/)), tt)*tt, (/3, na/))
        !print*, 'THIS BIT = ',norm2(reshape(dot_product(reshape(rp%cx(i)%dvdr,(/nd/)),tt)*tt,(/3,na/)))
      else
        ! climbing image force
        rp%cx(i)%force(:, :) = - rp%cx(i)%dvdr + 2.0d0 * &
            reshape(dot_product(reshape(rp%cx(i)%dvdr, (/nd/)), tt)*tt, (/3, na/))
      endif
    enddo

    ! first and last bead
    !
    if (endflag ) then
      rp%cx(1)%force(:, :) = -rp%cx(1)%dvdr(:, :)
      rp%cx(rp%nimage)%force(:, :) = -rp%cx(rp%nimage)%dvdr(:, :)
    elseif(ci_flag) then
      if(1 == imax) then
        tp = reshape(rp%cx(2)%r - rp%cx(1)%r, (/nd/))
        tp = tp/norm2(tp)
        rp%cx(1)%force(:, :) = -rp%cx(1)%dvdr + 2.0d0 * &
            reshape(dot_product(reshape(rp%cx(1)%dvdr, (/nd/)), tp)*tp, (/3, na/))
      elseif (rp%nimage == imax) then
        tm = reshape(rp%cx(rp%nimage)%r - rp%cx(rp%nimage-1)%r, (/nd/))
        tm = tm/norm2(tm)
        rp%cx(rp%nimage)%force(:, :) = -rp%cx(rp%nimage)%dvdr + 2.0d0 * &
            reshape(dot_product(reshape(rp%cx(rp%nimage)%dvdr, (/nd/)), tm)*tm, (/3, na/))
      endif
    endif

    ! set constraints...
    !
    do i = 1, rp%nimage
      idof = 0
      do j = 1, rp%na
        if (rp%cx(i)%fixedatom(j)) then
          rp%cx(i)%force(1:3, j) = 0.0d0
          idof = idof + 3
        else
          do k = 1, 3
            idof = idof + 1
            if (rp%cx(i)%fixeddof(idof)) then
              rp%cx(i)%force(k, j) = 0.0d0
            endif
          enddo
        endif
      enddo
    enddo

    return
  end Subroutine GetProjForces3


  !
  !*************************************************************************
  !
  !> ShimmyEndBeads
  !!
  !!
  !! - rp: Reaction path object.
  !! - kspring: NEB spring constant (in au), used as a reference for variable spring
  !!   size is made.
  !!
  !
  !*************************************************************************
  !
  Subroutine ShimmyEndBeads(rp, shimmied, beadcount, pend)
    implicit none
    integer :: i, j, k, l, imax, beadcount, pend, ne, ni
    type (rxp) :: rp
    double precision :: thresh, totlen, newlen, clen, olen, dell, rr(rp%nimage,3,rp%na)
    double precision :: dir(rp%na*3), tail(rp%na*3), front(rp%na*3), leng(rp%nimage), nlen, c
    logical :: success, flatbeads(rp%nimage), shimmied

    flatbeads = .false.
    !thresh = 0.0010d0
    thresh = 0.8E-004
    ! if gradients and derivatives are flat...
    !
    leng(1) = abs(norm2(reshape(rp%cx(1)%r, shape(dir)) - reshape(rp%cx(2)%r, shape(dir))))
    totlen = leng(1)
    do i = 2, rp%nimage-1
      !Call AbInitio(rp%cx(i), 'grad', success)
      leng(i) = abs(norm2(reshape(rp%cx(i)%r, shape(dir)) - reshape(rp%cx(i+1)%r, shape(dir))))
      if (abs((norm2(reshape(rp%cx(i-1)%dvdr, shape(dir))) - norm2(reshape(rp%cx(i)%dvdr, shape(dir)))*2.0d0 + &
            norm2(reshape(rp%cx(i+1)%dvdr, shape(dir))))/(leng(i-1) + leng(i))**2) <= thresh) then
        if (abs((rp%cx(i-1)%vcalc - rp%cx(i)%vcalc*2.0d0 + rp%cx(i+1)%vcalc)/(leng(i-1) + leng(i))**2) <= thresh*0.1) then
          flatbeads(i) = .true.
        endif
      endif
      totlen = totlen+leng(i)
      ! for debugging purposes..
      !print*, 'I = ', i, 'FLATBEADS = ', flatbeads(i)
      !print*, 'D = ', abs((norm2(reshape(rp%cx(i-1)%dvdr,shape(dir)))- norm2(reshape(rp%cx(i)%dvdr,shape(dir)))*2.0d0+ &
      !         norm2(reshape(rp%cx(i+1)%dvdr,shape(dir))))/(leng(i-1) + leng(i))**2), thresh
      !print*, 'V = ' ,abs((rp%cx(i-1)%vcalc- rp%cx(i)%vcalc*2.0d0+ rp%cx(i+1)%vcalc)/(leng(i-1)+ leng(i))**2), thresh*0.10
      !print*, rp%cx(i-1)%vcalc, rp%cx(i)%vcalc, rp%cx(i+1)%vcalc
    enddo

    pend = 0
    if (flatbeads(2)) pend = 1
    if (flatbeads(rp%nimage-1)) pend = pend + rp%nimage

    beadcount = 0
    ni = 1
    ! only shimmy up to half from the left...
    do i = 2, rp%nimage/2
      if (flatbeads(i)) then
        beadcount = beadcount + 1
        ni = i ! + 1
      else
        exit
      endif
    enddo
    ne = rp%nimage
    ! only shimmy up to  half from the right
    do i = rp%nimage-1, rp%nimage/2, -1
      if (flatbeads(i)) then
        beadcount = beadcount + 1
        ne = i !- 1
      else
        exit
      endif
    enddo

    !print*, 'NI = ', ni, 'Ne = ', ne, 'NIMAGE = ', rp%nimage
    shimmied = .false.
    if (beadcount >= 2) then
      if (pend > 0) shimmied = .true.
      do i = 1, rp%nimage
        rr(i, 1:3, 1:rp%na) = rp%cx(i)%r(1:3, 1:rp%na)
      enddo
      newlen = +0.00010d0
      do j = ni+1, ne
        newlen = newlen + leng(j-1)
      enddo
      dell = newlen/(rp%nimage-1)
      clen = 0.0d0
      nlen = 0.0d0
      !print*, 'TOTAL LENGTH = ', totlen, 'NEWLEN = ', newlen, 'DEL = ', dell
      do i = 1, rp%nimage
        clen = dell*(i-1)
        olen = 0.0d0
        !print*, 'CLEN = ', clen, 'I = ', i
        do j = ni+1, ne
          olen = olen + leng(j-1)
          if ( olen > clen ) then
            c = (clen-(olen-leng(j-1)))/leng(j-1)
            rp%cx(i)%r = rr(j-1, :, :) + (rr(j, :, :) - rr(j-1, :, :))*0.50d0*(1.0d0+tanh((c-1.0d0)*5.0d0))
            !rp%cx(i)%r = rr(j-1,:,:) + (rr(j,:,:) - rr(j-1,:,:))*(clen-(olen-leng(j-1)))/leng(j-1)
            ! for debugging purposes..
            !print*, 'FOUND AT OLEN = ', olen, 'J = ', j
            !print*, 'DIFF = ', clen-(olen-leng(j-1)), 'LENG = ', leng(j-1), 'RATIO = ', (clen-(olen-leng(j-1)))/leng(j-1)
            !if (i == 0) print*, 'LEN = 0', 'nlen =  0'
            !if (i  > 1) nlen = nlen +abs(norm2(rp%cx(i)%r-rp%cx(i-1)%r))
            !if (i  > 1) print*, 'LEN = ', abs(norm2(rp%cx(i)%r-rp%cx(i-1)%r)), 'NLEN = ', nlen
            exit
          endif
          !if (j .eq. ne .and. i .eq. rp%nimage) print*, 'LAST ADDEED', ne,rp%nimage
          if (j .eq. ne .and. i .eq. rp%nimage) &
            rp%cx(i)%r = rr(ne, :, :)
        enddo
        call AbInitio(rp%cx(i), 'grad', success)
        rp%cx(i)%p = 0.0d0
      enddo
      !print*, 'FINAL DIFFERENCE = ',rp%nimage, norm2(reshape(rp%cx(rp%nimage)%r-rr(ne,:,:),(/3*rp%na/)))
    endif
    return
  end subroutine ShimmyEndBeads


  !
  !*************************************************************************
  !
  !> GetFourierForces
  !!
  !! Calculates the Fourier expansion coefficients for the path rp.
  !!
  !! - rp: Input reaction-path object.
  !!
  !*************************************************************************
  !
  Subroutine GetFourierForces(rp)
    type(rxp) :: rp
    integer :: i, nb, k, j, m, na
    real(8) :: lambda

    nb = rp%nimage
    na = rp%na

    ! First, calculate the force on the end-points, including the contribution
    ! from the string correlation along the reaction-path.
    ! Note that the calculated derivatives on each image are passed in...
    !
    do i = 2, nb           ! Start at 2 to avoid overcounting i = 1 contribution.
      lambda = dble(i-1) / dble(nb-1)
      do j = 1, na
        do k = 1, 3
          rp%cx(1)%dvdr(k, j) = rp%cx(1)%dvdr(k, j) + (1.d0 - lambda) * rp%cx(i)%dvdr(k, j)
          rp%cx(nb)%dvdr(k, j) = rp%cx(nb)%dvdr(k, j) + lambda * rp%cx(i)%dvdr(k, j)
        enddo
      enddo
    enddo

    ! Divide through by total number of beads.
    !
    do i = 1, nb
      do j = 1, na
        do k = 1, 3
          rp%cx(i)%dvdr(k, j) = rp%cx(i)%dvdr(k, j) / dble(nb)
        enddo
      enddo
    enddo

    ! Now calculate the forces on the Fourier coefficients.
    !
    do m = 1, nb
      do j = 1, na
        do k = 1, 3
          rp%dcoeff(k, j, m) = 0.d0
          do i = 1, nb
            lambda = dble(i-1) / dble(nb-1)
            rp%dcoeff(k, j, m) = rp%dcoeff(k, j, m) + rp%cx(i)%dvdr(k, j) * sin(dble(m) * pivalue * lambda)
          enddo
        enddo
      enddo
    enddo

    return
  end Subroutine GetFourierForces


  !
  !*************************************************************************
  !
  !> FourierToPath
  !!
  !! Calculates the path positions using the path Fourier coefficients.
  !!
  !! - rp: Input reaction-path object.
  !!
  !*************************************************************************
  !
  Subroutine FourierToPath(rp)
    type(rxp) :: rp
    integer :: i, j, idof, nb, natom, k, m
    real(8) :: lambda, re, rs

    nb = rp%nimage
    natom = rp%na

    do i = 2, nb-1
      lambda = dble(i-1) / dble(nb-1)

      idof = 0
      do j = 1, natom
        if (.not. rp%cx(i)%fixedatom(j)) then
          do k = 1, 3
            idof = idof + 1
            if (.not. rp%cx(i)%fixeddof(idof)) then

              re = rp%cx(nb)%r(k, j)
              rs = rp%cx(1)%r(k, j)
              rp%cx(i)%r(k, j) = rp%cx(1)%r(k, j) + lambda * (re - rs)

              do m = 1, nb
                rp%cx(i)%r(k, j) = rp%cx(i)%r(k, j) + rp%coeff(k, j, m) * sin(dble(m) * pivalue * lambda )
              enddo

            endif
          enddo
        else
          idof = idof + 3
        endif
      enddo
    enddo

    return
  end Subroutine FourierToPath


  !
  !*************************************************************************
  !
  !> PathToFourier
  !!
  !! Calculates the Fourier coefficiencts from path positions.
  !!
  !! - rp: Input reaction-path object.
  !!
  !*************************************************************************
  !
  Subroutine PathToFourier(rp)
    type(rxp) :: rp
    integer :: i, j, idof, nb, natom, k, m
    real(8) :: lambda, re, rs, sum, dlambda, fac

    nb = rp%nimage
    natom = rp%na

    dlambda = 1.d0 / dble(nb-1)

    do j = 1, natom
      do k = 1, 3
        do m = 1, nb
          rp%coeff(k, j, m) = 0.d0

          sum = 0.d0

          do i = 1, nb
            lambda = dble(i-1) / dble(nb-1)
            fac = sin(m * pivalue * lambda)
            sum = sum + (rp%cx(i)%r(k, j) - rp%cx(1)%r(k, j) - lambda*(rp%cx(nb)%r(k, j) - rp%cx(1)%r(k, j)))*fac
          enddo
          sum = sum * 2.0 * dlambda
          rp%coeff(k, j, m) = sum
        enddo
      enddo
    enddo

    return
  end Subroutine PathToFourier


  !
  !*************************************************************************
  !
  !> FirstVVUpdate
  !!
  !! Evolves the Fourier coefficients and momenta through first half of a
  !! Velocity-Verlet time-step.
  !!
  !! - rp: Reaction path object
  !! - timestep: obvs...
  !! - mass: Fourier coefficient mass.
  !!
  !*************************************************************************
  !
  !
  Subroutine FirstVVUpdate(rp, timestep, mass)
    implicit none
    type(rxp) :: rp
    real(8) :: timestep, mass
    integer :: i, j, k, nb, idof

    ! Update the end-points
    !
    nb = rp%nimage
    idof = 0
    do j = 1, rp%na
      if (.not. rp%cx(1)%fixedatom(j)) then
        do k = 1, 3
          idof = idof + 1
          if (.not. rp%cx(1)%fixeddof(idof)) then
            rp%cx(1)%p(k, j) = rp%cx(1)%p(k, j) - 0.5d0 * timestep * rp%cx(1)%dvdr(k, j)
            rp%cx(nb)%p(k, j) = rp%cx(nb)%p(k, j) - 0.5d0 * timestep * rp%cx(nb)%dvdr(k, j)
            rp%cx(1)%r(k, j) = rp%cx(1)%r(k, j) + timestep * rp%cx(1)%p(k, j)/rp%cx(1)%mass(j)
            rp%cx(nb)%r(k, j) = rp%cx(nb)%r(k, j) + timestep * rp%cx(nb)%p(k, j)/rp%cx(nb)%mass(j)
          endif
        enddo
      else
        idof = idof + 3
      endif
    enddo

    ! Update coefficients.
    !
    do i = 1, rp%nimage
      idof = 0
      do j = 1, rp%na
        if (.not. rp%cx(i)%fixedatom(j)) then
          do k = 1, 3
            idof = idof + 1
            if (.not. rp%cx(i)%fixeddof(idof)) then
              rp%pcoeff(k, j, i) = rp%pcoeff(k, j, i) - 0.5 * timestep * rp%cx(i)%dvdr(k, j)
              rp%coeff(k, j, i) = rp%coeff(k, j, i) + rp%pcoeff(k, j, i) * timestep / mass
            endif
          enddo
        else
          idof = idof + 3
        endif
      enddo
    enddo

    return
  end Subroutine FirstVVUpdate


  !
  !*************************************************************************
  !
  !> SecondVVUpdate
  !!
  !! Evolves the Fourier coefficients and momenta through second half of a
  !! Velocity-Verlet time-step.
  !!
  !! - rp: Reaction path object
  !! - timestep: obvs...
  !! - mass: Fourier coefficient mass.
  !!
  !*************************************************************************
  !
  Subroutine SecondVVUpdate(rp, timestep, mass)
    implicit none
    type(rxp) :: rp
    real(8) :: timestep, mass
    integer :: i,j,k,nb,idof

    ! Update the end-points
    !
    nb = rp%nimage
    idof = 0
    do j = 1, rp%na
      if (.not. rp%cx(1)%fixedatom(j)) then
        do k = 1, 3
          idof = idof + 1
          if (.not. rp%cx(1)%fixeddof(idof)) then
            rp%cx(1)%p(k, j) = rp%cx(1)%p(k, j) - 0.5d0 * timestep * rp%cx(1)%dvdr(k, j)
            rp%cx(nb)%p(k, j) = rp%cx(nb)%p(k, j) - 0.5d0 * timestep * rp%cx(nb)%dvdr(k, j)
          endif
        enddo
      else
        idof = idof + 3
      endif
    enddo

    ! Update coefficients.
    !
    do i = 1, rp%nimage
      idof = 0
      do j = 1, rp%na
        if (.not. rp%cx(i)%fixedatom(j)) then
          do k = 1, 3
            idof = idof + 1
            if (.not. rp%cx(i)%fixeddof(idof)) then
              rp%pcoeff(k, j, i) = rp%pcoeff(k, j, i) - 0.5 * timestep * rp%cx(i)%dvdr(k, j)
            endif
          enddo
        else
          idof = idof + 3
        endif
      enddo
    enddo

    return
  end Subroutine SecondVVUpdate


  !
  !*************************************************************************
  !
  !> GetSpringForces
  !!
  !! Calculates the inter-bead spring potential and forces.
  !!
  !! - rp: Reaction path object.
  !! - kspring: Spring constant
  !!
  !*************************************************************************
  !
  Subroutine GetSpringForces(rp, kspring)
    implicit none
    type(rxp) :: rp
    integer :: i, j, k, nb, na
    real(8) :: dr, kspring

    na = rp%na
    nb = rp%nimage
    rp%vspring = 0.d0

    do i = 2, nb
      do j = 1, na
        do k = 1, 3

          dr = rp%cx(i)%r(k, j) - rp%cx(i-1)%r(k, j)
          rp%vspring = rp%vspring + kspring * dr * dr

          rp%cx(i)%dvdr(k, j) = rp%cx(i)%dvdr(k, j) + 2.d0 * kspring * dr
          rp%cx(i-1)%dvdr(k, j) = rp%cx(i-1)%dvdr(k, j) - 2.d0 * kspring * dr

        enddo
      enddo
    enddo

    return
  end Subroutine GetSpringForces


  !
  !*************************************************************************
  !
  !> StripInactiveFromPath
  !!
  !! Removes inactive molecules from the reaction path. Here, inactive means
  !! that the molecules are not involved in any bond changes.
  !!
  !! - rp: Reaction path object.
  !! - UpdateCons: flag determining if the FixedAtom and FixedDof should be updated to match the stripped CXS
  !! - stripfile: The filename where the stripped path is output.
  !!
  !*************************************************************************
  !
  Subroutine StripInactiveFromPath(rp, stripfile, FixedDOF, FixedAtom, NDOFconstr, Natomconstr, UpdateCons)
    implicit none
    type(rxp) :: rp
    integer :: i, j, k, nb, na, isum, l, icount
    integer :: tmpFixedDOF(3*NAMAX), tmpNDOFconstr, tmpFixedAtom(NAMAX), tmpNatomconstr
    integer :: FixedDOF(3*NAMAX), NDOFconstr, FixedAtom(NAMAX), Natomconstr, tmpatomidx(3)
    logical :: UpdateCons, isochk
    logical, dimension(:), allocatable :: keep_atom, keep_start, keep_end, found
    character (len=*) :: stripfile

    na = rp%na
    nb = rp%nimage
    allocate(keep_atom(na), keep_start(NMOLMAX), keep_end(NMOLMAX), found(na))

    ! First, get graphs and molecules for end-points.
    !
    call GetGraph(rp%cx(1))
    call GetGraph(rp%cx(nb))

    call GetMols(rp%cx(1))
    call GetMols(rp%cx(nb))

    ! Identify which atoms have changed bonding.
    !
    keep_atom(:) = .true.
    do k = 1, na
      isum = 0
      do i = 1, na
        isum = isum + abs(rp%cx(1)%graph(k, i) - rp%cx(nb)%graph(k, i))
      enddo
      if (isum == 0) keep_atom(k) = .false.
    enddo

    ! Identify which molecules need to be kept...
    !
    keep_start(:) = .false.
    keep_end(:) = .false.
    do k = 1, na
      if (keep_atom(k)) then
        do i = 1, rp%cx(1)%nmol
          do j = 1, rp%cx(1)%namol(i)
            if (rp%cx(1)%molid(i, j) == k) keep_start(i) = .true.
          enddo
        enddo
        do i = 1, rp%cx(nb)%nmol
          do j = 1, rp%cx(nb)%namol(i)
            if (rp%cx(nb)%molid(i, j) == k) keep_end(i) = .true.
          enddo
        enddo
      endif
    enddo

    ! Check pairs of molecules and see if they are the same on both sides, if
    ! they
    do i = 1, rp%cx(1)%nmol ; if (keep_start(i)) cycle
      do j = 1, rp%cx(nb)%nmol ; if (keep_end(j)) cycle
        if (rp%cx(1)%namol(i) .eq. rp%cx(nb)%namol(j)) then
          if (sum(abs(rp%cx(1)%molid(i, :rp%cx(1)%namol(i)) - rp%cx(nb)%molid(j, :rp%cx(nb)%namol(j)))) .eq. 0) then
!           if (.not.MoleculesAreSame(rp%cx(1),i,rp%cx(nb),j)) then
            keep_start(i) = .true.
            keep_end(j) = .true.
!           endif
          endif
        endif
      enddo
    enddo

    ! Set flags indicating that atoms have been allocated.
    !
    found(1:na) = .false.
    do i = 1, rp%cx(1)%nmol
      if (keep_start(i)) then
        do j = 1, rp%cx(1)%namol(i)
          l = rp%cx(1)%molid(i, j)
          found(l) = .true.
        enddo
      endif
    enddo
    do i = 1, rp%cx(nb)%nmol
      if (keep_end(i)) then
        do j = 1, rp%cx(nb)%namol(i)
          l = rp%cx(nb)%molid(i, j)
          if (.not.found(l)) then
            found(l) = .true.
          endif
        enddo
      endif
    enddo

    ! Count the number of atoms remaining.
    !
    icount = 0
    do i = 1, na
      if (found(i)) icount = icount + 1
    enddo

    if (UpdateCons) then
      !  Update the contraints matrix:
      !
      tmpFixedDOF = 0
      tmpNDOFconstr = 0
      l = 0
      do i = 1, na
        if (found(i)) then
          do k = 3*(i-1)+1, 3*(i-1)+3
            l = l + 1
            do j = 1, NDOFconstr
              if (FixedDOF(j) .eq. k) then
                tmpNDOFconstr = tmpNDOFconstr +1
                tmpFixedDOF(tmpNDOFconstr) = l
              endif
            enddo
          enddo
        endif
      enddo
      FixedDOF = tmpFixedDOF
      NDOFconstr = tmpNDOFconstr
      tmpFixedAtom = 0
      tmpNatomconstr = 0
      k = 0
      do i = 1, na
        if (found(i)) then
          k = k + 1
          do j = 1, Natomconstr
            if (FixedAtom(j) .eq. i) then
              tmpNatomconstr = tmpNatomconstr + 1
              tmpFixedAtom(tmpNatomconstr) = k
            endif
          enddo
        endif
      enddo
      FixedAtom = tmpFixedAtom
      Natomconstr = tmpNatomconstr
    endif
    if (sum(atomidx) .gt. 0) then
      ! update the atomidx NEB array...
      !
      k = 0
      tmpatomidx = 0
      do i = 1, na
        if (found(i)) then
          k = k + 1
          if (i .eq. atomidx(1)) then
            tmpatomidx(1) = k
          endif
          if (i .eq. atomidx(2)) then
            tmpatomidx(2) = k
          endif
          if (i .eq. atomidx(3)) then
            tmpatomidx(3) = k
          endif
        endif
      enddo
      if (tmpatomidx(1) .eq. 0 .or. tmpatomidx(2) .eq. 0 .or. tmpatomidx(3) .eq. 0) then
        write(logfile, '("Warning: could not find one of the aligned atoms indecies, switching to automatic")')
        atomidx = (/0, 0, 0/)
      else
        atomidx = tmpatomidx
      endif
    endif
    ! Create a new reaction path with the kept atoms and output this. If required, we
    ! can then create a new path by reading in the stripped path later.
    !
    open(31, file=stripfile, status='unknown')
    do i = 1, nb

      ! Write number of atoms and header
      !
      write(31, *) icount
      write(31, *)

      do j = 1, na
        if (found(j)) then
          write(31, *) rp%cx(i)%atomlabel(j), rp%cx(i)%r(1, j)*bohr_to_ang, &
              rp%cx(i)%r(2, j)*bohr_to_ang, &
              rp%cx(i)%r(3, j)*bohr_to_ang
        endif
      enddo
    enddo
    close(31)

    return
  end Subroutine StripInactiveFromPath


  !
  !*************************************************************************
  !
  !> errorr
  !!  calculates the absolute difference, elementwise between two matrices.
  !!
  !! - one, two: the matrices
  !! - na: ...
  !!
  !*************************************************************************
  !
  function errorr(one,two,na)
    integer :: i, j, k, l, m, n, nd, nb, na
    real(8) :: errorr
    real(8), dimension(3, na) :: one, two
    errorr = 0.0d0
    do i = 1, na
      do j = 1, 3
        errorr = errorr + abs(one(j, i)-two(j, i))
      enddo
    enddo

  end function


end Module rpath
