
!***************************************************************************************
!
!> @brief constants definition file, containing importance universal values used
!! throughout code.
!!
!! **You should avoid changing these values unless you are absolutely sure you know
!! what you're doing.**
!
!***************************************************************************************
!

Module constants

  implicit none

  !
  ! Useful constants.
  !
  integer, parameter :: NAMAX = 1000           !< Maximum number of atoms.
  integer, parameter :: NFMAX = 3 * NAMAX     !< Maximum number of DOFs.
  integer, parameter :: NBMAX = 20            !< Maximum number of images.
  integer, parameter :: NMOLMAX = 50         !< Maximum number of molecules.
  integer, parameter :: NDIRMAX = 50          !< Maximum number of directories to create.
  integer, parameter :: NTYPEMAX = 100        !< Maximum number of atom types.
  integer, parameter :: NTRYMAX = 1000        !< Maximum number of outer graph-move attempts.
  integer, parameter :: NTRYMAX2 = 5000       !< Maximum number of inner graph-move attempts.
  integer, parameter :: NMOVEMAX = 20        !< Maximum number of graph-moves
  integer, parameter :: NFORBIDMAX = 50        !< Maximum number of forbidden graphs.
  integer, parameter :: NAMOVEMAX = 6        !< Maximum number of atoms which can be involved in graph moves.
  integer, parameter :: MAXMOL = 4000 !< Maximum number of molecules stored in MolData datatype
  integer, parameter :: NMECHMAX = 250        !< Maximum number of stored mechanisms
  integer, parameter :: NRXNMAX = 25        !< Maximum number of stored mechanisms
  integer, parameter :: nvalmax = 12
  integer, parameter :: NELMAX = 10
  integer, parameter :: NBONDMAX = 100
  integer, parameter :: NANGLEMAX = 200
  integer, parameter :: NTORSMAX = 300

  !
  ! Useful numerical constants.
  !
  real(8), parameter :: bohr_to_ang   = 0.5291772108d0 !< Bohr to Angstrom conversion
  real(8), parameter :: ang_to_bohr   = 1.889726128d0  !< Angstrom to Bohr conversion
  real(8), parameter :: kcalmol_to_au = 0.001593601d0  !< Convert kcal/mol to au.
  real(8), parameter :: au_to_kjmol   = 2625.49963d0   !< convert au to kJ/mol
  real(8), parameter :: au_to_ev      = 27.211396d0
  real(8), parameter :: ev_to_au      = 0.036749322d0
  real(8), parameter :: ev_ang_to_au_bohr = 0.019446864d0
  real(8), parameter :: pivalue       = 3.141592654d0  !< PI!
  real(8), parameter :: Kboltz        = 3.166829d-6    !< Boltzmann constant in atomic units.
  real(8), parameter :: hbar          = 1.d0           !< hbar = 1 in atomic units
  real(8), parameter :: fs_to_au      = 41.341373d0    !< Convert femtosecond to atomic units
  real(8), parameter :: Eh_to_Kelvin  = 3.157733d+5    !< Convert Hartrees to Kelvin
  real(8), parameter :: Eh_to_wavenumber = 219474.63   !< Convert Hartrees to wavenumbers.
  real(8), parameter :: BIG           = 10000          !< A 'BIG' number which is occasionally useful.
  real(8), parameter :: SMALL         = 1d-5           !< A 'SMALL' number which is occasionally useful.
  real(8), parameter :: me_to_amu = 5.485799090d-4     !< Convert electron mass units to atomic mass units.
  real(8), parameter :: epsil = 10.0d0**(-20)
  real(8), parameter :: GCONV = 1d-7


  !> Set atomic masses. The MASS index can be referenced by atomic number;
  !! For example, MASS(1) = mass of hydrogen.
  !!
  real(8),parameter :: MASS(87) = (/1837.1527d0,  &      ! H
                                    0.d0,          &    ! He
                                    0.d0,          &    ! Li
                                     0.d0,       &    ! Be
                                    19852.3d0,      &    ! B
                                    21894.16747d0, &    ! C
                                    25533.1990d0,   &    ! N
                                    29156.9471d0, &      ! O
                                    34905.9013d0, &             ! F
                                    0.d0, &             ! Ne
                                    0.d0, &             ! Na
                                    0.d0, &             ! Mg
                                    0.d0, &             ! Al
                                    51594.6d0, &             ! Si
                                    56903.5d0, &             ! P
                                    0.d0, &             ! S
                                    0.d0, &             ! Cl
                                    0.d0, &             ! Ar
                                    0.d0, &             ! K
                                    0.d0, &             ! Ca
                                    0.d0, &             ! Sc
                                    0.d0, &             ! Ti
                                    0.d0, &             ! V
                                    0.d0, &             ! Cr
                                    0.d0, &             ! Mn
                                    102595.8d0, &             ! Fe
                                    107524.4367d0, &           ! Co
                                    0.d0, &             ! Ni
                                    0.d0, &             ! Cu
                                    0.d0, &             ! Zn
                                    0.d0, &             ! Ga
                                    0.d0, &             ! Ge
                                    0.d0, &             ! As
                                    0.d0, &             ! Se
                                    0.d0, &             ! Br
                                    0.d0, &             ! Kr
                                    0.d0, &             ! Rb
                                    0.d0, &             ! Sr
                                    0.d0, &             ! Y
                                    0.d0, &             ! Zr
                                    0.d0, &             ! Nb
                                    0.d0, &             ! Mo
                                    0.d0, &             ! Tc
                                    0.d0, &             ! Ru
                                    195509.8d0, &             ! Rh
                                    195509.8d0, &             ! Pd
                                    0.d0, &             ! Ag
                                    0.d0, &             ! Cd
                                    0.d0, &             ! In
                                    0.d0, &             ! Sn
                                    0.d0, &             ! Sb
                                    0.d0, &             ! Te
                                    0.d0, &             ! I
                                    0.d0, &             ! Xe
                                    0.d0, &             ! Cs
                                    0.d0, &             ! Ba
                                    0.d0, &             ! La
                                    0.d0, &             ! Ce
                                    0.d0, &             ! Pr
                                    0.d0, &             ! Nd
                                    0.d0, &             ! Pm
                                    0.d0, &             ! Sm
                                    0.d0, &             ! Eu
                                    0.d0, &             ! Gd
                                    0.d0, &             ! Tb
                                    0.d0, &             ! Dy
                                    0.d0, &             ! Ho
                                    0.d0, &             ! Er
                                    0.d0, &             ! Tm
                                    0.d0, &             ! Yb
                                    0.d0, &             ! Lu
                                    0.d0, &             ! Hf
                                    0.d0, &             ! Ta
                                    0.d0, &             ! W
                                    0.d0, &             ! Re
                                    0.d0, &             ! Os
                                    0.d0, &             ! Ir
                                    358399.0973d0, &             ! Pt
                                    0.d0, &             ! Au
                                    0.d0, &             ! Hg
                                    0.d0, &             ! Tl
                                    0.d0, &             ! Pb
                                    0.d0, &             ! Bi
                                    0.d0, &             ! Po
                                    0.d0, &             ! At
                                    0.d0, &             ! Rn
                                    50000.0d0/)        ! G1



  !> Set covalent radii. The CovRad(:) array is indexed by atomic number, so CovRad(1) is the
  !! covalent radius of hydrogen. Note that the values are given in Angstroms and
  !! converted into Bohr using the ang_to_bohr conversion.
  !
  real(8), parameter :: CovRad(87) = (/0.4d0 * ang_to_bohr, &     ! H !! 0.4?
                                      0.d0, &                      ! He
                                     0.d0, &                      ! Li
                                     0.d0, &                      ! Be
                                     0.84d0 * ang_to_bohr, &      ! B
                                     0.72d0 * ang_to_bohr, &      ! C      ! 0.72
                                     0.72d0 * ang_to_bohr, &      ! N
                                     0.62d0 * ang_to_bohr, &        ! O
                                    !! 0.72d0 * ang_to_bohr, &        ! O
                                     0.6d0 * ang_to_bohr, &     ! F
                                     0.d0, &             ! Ne
                                    0.d0, &             ! Na
                                    0.d0, &             ! Mg
                                    0.d0, &             ! Al
                                    1.15d0 * ang_to_bohr, &             ! Si
                                    1.07d0 * ang_to_bohr, &             ! P
                                    1.07d0*ang_to_bohr, &             ! S
                                    0.d0, &             ! Cl
                                    0.d0, &             ! Ar
                                    0.d0, &             ! K
                                    0.d0, &             ! Ca
                                    0.d0, &             ! Sc
                                    0.d0, &             ! Ti
                                    0.d0, &             ! V
                                    0.d0, &             ! Cr
                                    0.d0, &             ! Mn
                                    1.52d0 * ang_to_bohr, &             ! Fe
                                    1.52d0 * ang_to_bohr, &             ! Co - 1.52
                                    0.d0, &             ! Ni
                                    0.d0, &             ! Cu
                                    0.d0, &             ! Zn
                                    0.d0, &             ! Ga
                                    0.d0, &             ! Ge
                                    0.d0, &             ! As
                                    0.d0, &             ! Se
                                    0.d0, &             ! Br
                                    0.d0, &             ! Kr
                                    0.d0, &             ! Rb
                                    0.d0, &             ! Sr
                                    0.d0, &             ! Y
                                    0.d0, &             ! Zr
                                    0.d0, &             ! Nb
                                    0.d0, &             ! Mo
                                    0.d0, &             ! Tc
                                    0.d0, &             ! Ru
                                    1.42d0 * ang_to_bohr, &             ! Rh
                                    1.39d0 * ang_to_bohr, &             ! Pd
                                    0.d0, &             ! Ag
                                    0.d0, &             ! Cd
                                    0.d0, &             ! In
                                    0.d0, &             ! Sn
                                    0.d0, &             ! Sb
                                    0.d0, &             ! Te
                                    0.d0, &             ! I
                                    0.d0, &             ! Xe
                                    0.d0, &             ! Cs
                                    0.d0, &             ! Ba
                                    0.d0, &             ! La
                                    0.d0, &             ! Ce
                                    0.d0, &             ! Pr
                                    0.d0, &             ! Nd
                                    0.d0, &             ! Pm
                                    0.d0, &             ! Sm
                                    0.d0, &             ! Eu
                                    0.d0, &             ! Gd
                                    0.d0, &             ! Tb
                                    0.d0, &             ! Dy
                                    0.d0, &             ! Ho
                                    0.d0, &             ! Er
                                    0.d0, &             ! Tm
                                    0.d0, &             ! Yb
                                    0.d0, &             ! Lu
                                    0.d0, &             ! Hf
                                    0.d0, &             ! Ta
                                    0.d0, &             ! W
                                    0.d0, &             ! Re
                                    0.d0, &             ! Os
                                    0.d0, &             ! Ir
                                    1.46d0 * ang_to_bohr, &             ! Pt
                                    0.d0, &             ! Au
                                    0.d0, &             ! Hg
                                    0.d0, &             ! Tl
                                    0.d0, &             ! Pb
                                    0.d0, &             ! Bi
                                    0.d0, &             ! Po
                                    0.d0, &             ! At
                                    0.d0, &             ! Rn
                                    1.50d0/)              ! G1


  integer,parameter :: AValency(87) = (/1, &             ! H
                                       2, &             ! He
                                       1, &             ! Li
                                       2, &             ! Be
                                       3, &             ! B
                                       4, &             ! C
                                       5, &             ! N
                                       6, &             ! O
                                       7, &             ! F
                                       8, &             ! Ne
                                       1, &             ! Na
                                       2, &             ! Mg
                                       3, &             ! Al
                                       4, &             ! Si
                                       5, &             ! P
                                       6, &             ! S
                                       7, &             ! Cl
                                       8, &             ! Ar
                                       1, &             ! K
                                       2, &             ! Ca
                                       0, &             ! Sc
                                       0, &             ! Ti
                                       0, &             ! V
                                       0, &             ! Cr
                                       0, &             ! Mn
                                       0, &             ! Fe
                                       0, &             ! Co
                                       0, &             ! Ni
                                       0, &             ! Cu
                                       0, &             ! Zn
                                       0, &             ! Ga
                                       0, &             ! Ge
                                       0, &             ! As
                                       0, &             ! Se
                                       7, &             ! Br
                                       0, &             ! Kr
                                       0, &             ! Rb
                                       0, &             ! Sr
                                       0, &             ! Y
                                       0, &             ! Zr
                                       0, &             ! Nb
                                       0, &             ! Mo
                                       0, &             ! Tc
                                       0, &             ! Ru
                                       0, &             ! Rh
                                       0, &             ! Pd
                                       0, &             ! Ag
                                       0, &             ! Cd
                                       0, &             ! In
                                       0, &             ! Sn
                                       0, &             ! Sb
                                       0, &             ! Te
                                       7, &             ! I
                                       0, &             ! Xe
                                       0, &             ! Cs
                                       0, &             ! Ba
                                       0, &             ! La
                                       0, &             ! Ce
                                       0, &             ! Pr
                                       0, &             ! Nd
                                       0, &             ! Pm
                                       0, &             ! Sm
                                       0, &             ! Eu
                                       0, &             ! Gd
                                       0, &             ! Tb
                                       0, &             ! Dy
                                       0, &             ! Ho
                                       0, &             ! Er
                                       0, &             ! Tm
                                       0, &             ! Yb
                                       0, &             ! Lu
                                       0, &             ! Hf
                                       0, &             ! Ta
                                       0, &             ! W
                                       0, &             ! Re
                                       0, &             ! Os
                                       0, &             ! Ir
                                       0, &             ! Pt
                                       0, &             ! Au
                                       0, &             ! Hg
                                       0, &             ! Tl
                                       0, &             ! Pb
                                       0, &             ! Bi
                                       0, &             ! Po
                                       0, &             ! At
                                       0, &             ! Rn
                                       10 /)            ! G1



!  real(8), parameter :: bondingsf = 1.10d0       !< Scale factor applied to covalent radii to
  ! lammps optimises wrong C-H distance if made 1.1
  real(8), parameter :: bondingsf = 1.10d0       !< Scale factor applied to covalent radii to
                                                !! define bonding. Atoms are bonded if
                                                !! r(i,j) <= (covrad(i) + covrad(j)) * bondingsf

  real(8), parameter :: bondingrange1 = 0.25d0 * ang_to_bohr   !< Shift range over which atoms are restrained
                                                             !! in GDS simulations.

  real(8), parameter :: bondingrange2 = -0.10d0 * ang_to_bohr   !< Shift range over which atoms are restrained
                                                             !! in GDS simulations.
  !
  real(8), parameter :: RADIUS_MAX = 50.d0 * ang_to_bohr  ! 100?
  real(8), parameter :: RADIUS_MIN = 12.d0 * ang_to_bohr

  real(8), parameter :: LATTICESTEP = 4.d0 * ang_to_bohr

  real(8), parameter :: HIDPP = 2.750 ! au... its not clear..?
  !real(8), parameter :: HIDPP = 0.750 ! au... its not clear..?


  real(8) :: sbstrength(100,100) = 0.0d0
  real(8) :: dbstrength(10,10)  = 0.0d0
  real(8) :: tbstrength(10,10) = 0.0d0


  real(8), parameter :: GRPMINTHRESH = 1d-3      !< Force convergence target (in au) for minimization under GRP.
  real(8), parameter :: GRPMAXTHRESH = 1d-3      !< Force convergence target (in au) for minimization under GRP.


end Module constants
