!
!***************************************************************************************
!
!> @brief Contains global data definitions for CDE.
!
!***************************************************************************************
!

Module globaldata
  use constants
  implicit none

  ! Set the integer filenumber for the log file.
  !
  integer :: logfile = 11

  ! Filenames.
  character (len=50) :: inputfile           !< Input filename
  character (len=50) :: startfile, endfile  !< Filenames for start/end-points of new reaction path
  character (len=50) :: pathfile            !< Filename for full reaction-path for reading.
  character (len=50) :: movefile            !< Filename for graph moves.
  character (len=50) :: forbidfile          !< Filename for forbidden graphs.


  ! EVB parameters:
  integer :: evbtype, evbiter, nevbl, nevbdist,nevbqm
  logical :: evbvrep
  real(8) :: evbalpha1, evbalpha2, evbstep,evbmaxdl


  ! Global integer dummy:
  !
  integer :: idum

  ! System size definitions.
  !
  integer :: nimage                          !< Number of images in path
  integer :: f                               !< Number of DOFs
  integer :: na                              !< number of atoms

  ! Random numbers.
  !
  integer :: irun                            !< Random number seed


  ! Simulation parameters.
  !
  integer :: NEBiter               !< Total NEB/CINEB iterations
  integer :: NEBoutfreq            !< Output frequency during NEB.
  real(8) :: NEBconv               !< NEB convergence parameter for average force modulus.
  real(8) :: NEBMAXconv            !< NEB MAXIMUM convergence parameter for average force modulus.
  real(8) :: CIthresh              !< Force threshold to turn on climbing-image NEB.
  real(8) :: temperature, beta     !< Temperature and inverse temperature (1/kT)
  real(8) :: nebspring             !< Spring constant for NEB (au)
  real(8) :: NEBstep             !< Step-length for steepest-descent NEB
  character (len=8) :: calctype      !< Calculation type, NEB or GDS
  character (len=6) :: pathinit      !< Path initialization type, LINEAR
  character (len=6) :: pathoptmethod      !< Path optimization method, CINEB
  character (len=9) :: NEBmethod      !< (CI)NEB optimization method, CINEB
  logical :: startfrompath           !< Logical flag to signal whether or not we're starting from full path.
  logical :: stripinactive  !< Logical flag to signal that spectator molecules should be removed before optimization.
  logical :: optendsbefore !< Logical flag indicating whether to geometry-optimize path end-points before path optimization.
  logical :: optendsduring !< Logical flag indicating whether to optimize path end-points during path optimization.
  logical :: reconnect     !< Logical flag indicating whether to generate new linear path after end-point optimization.
  logical :: idpppath      !< Logical flag to turn the image-dependent pair potential NEB for a improved guess path
  logical :: fraginterpol  !< Logical flag to turn the fragment interpolation procedure - where fragments are first brought together before interpolating.
  logical :: optfragorient !< Logical flag to turn the fragment orientation optimization procedure - where fragments are oriented appropriately to enable the exchange of atoms

  ! Constrained DOFs and atoms.
  !
  integer :: NDOFconstr          !< Number of constrained DOFs
  integer :: FixedDOF(3*NAMAX)
  integer :: Natomconstr          !< Number of constrained DOFs
  integer :: FixedAtom(NAMAX)


  ! PES setup
  !
  character (len=6) :: PEStype             !< Identifies PES calculation type, 'ORCA'...
  character (len=25) :: PESfile            !< Identifies PES calculation template file.
  character (len=6) :: PESopttype          !< Identifies PES type for geometry optimisation.
  character (len=25) :: PESoptfile         !< Identifies geometry optimisation template file.
  character (len=100) :: PESexecutable      !< Identifies PES calculation executable.
  character (len=100) :: PESoptexecutable   !< Identifies geometry optimisation executable.
  logical :: PESfull                       !< Logical flag indicating whether we're running PES
                                           !! evaluations for the full system or independent molecules.


  ! Formats
  character (len=8), parameter :: fmt4 = '(I4.4)'


  ! GDS definitions:
  !
  integer :: nreactivetypes    !< Number of reactive atom-types.
  integer :: nfixtype          !< Number of types of fixed-bonds.
  character (len=2) :: reactivetype(NTYPEMAX)     !< Atom labels for reactive atoms.
  character(len = 2) :: fixedbondtype(NTYPEMAX,2) !< Atom labels for each type of fixed bonds
  logical :: reactive(NAMAX)                      !< Logical flag indicating whether or not atoms are reactive
  logical :: essentialmoves(NAMAX)                !< Logical flag indicating whether or not atoms are essential to moves
  logical :: essential(NAMAX)                     !< Logical flag indicating whether or not atoms are essential
  integer :: nessentialatoms                      !< Number of essential atoms labelled .TRUE. in essentialmoves(:)
  double precision ::intra_cutoff                 !< when performing a bond forming INTRA-molecular rxn move between atoms, set a cutoff distance

  integer :: nessentialatomsinmols
  integer :: essentialatomsinmols(NAMAX)      ! List of atoms which must be present in reacting molecules in
                                            ! a pthfinder calculation


                                                  !! part of molecules in a reaction
  logical :: fixedbonds(NAMAX,NAMAX)    !< Logical flag indicating fixed bonds.
  integer :: ngdsiter                    !< Number of GDS iterations
  real(8) :: gdsspring         !< Spring constant for inter-bead spring terms.
  real(8) :: gdsrestspring     !< Spring constant for bonded graph restraints
  real(8) :: nbstrength    !< Strength of non-bonded exponential graph-restraints
  real(8) :: nbrange       !< Range of non-bonded exponential graph-restraints
  real(8) :: kradius       !< Spring constant for intermolecular repulsion restraints.
  real(8) :: gdstemperature !< End-point temperature for GDS simulation.
  real(8) :: gdscoefftemp  !< Temperature for GDS coefficients
  real(8) :: gdscoeffmass  !< Mass for GDS coefficients
  real(8) :: timestep      !< GDS timestep
  real(8) :: VSthresh      !< Variable spring strength threshold (used suring CINEB)
  integer :: anebb         !< AutoNEB number of extra beads (also works as a logical switching Autoneb on)

  integer :: ngmove           !< Number of defined graph moves.
  integer :: namove(NMOVEMAX) !< Number of atoms in each graph move.
  integer :: gmstart(NMOVEMAX,NAMOVEMAX,NAMOVEMAX)     !< Start graph for each move
  integer :: gmend(NMOVEMAX,NAMOVEMAX,NAMOVEMAX)     !< End graph for each move
  character (len=4) :: movelabel(NMOVEMAX,NAMOVEMAX)  !< Labels for atoms involved in graph moves.
  character (len=7) :: movetype(NMOVEMAX)          !< For labelling as 'ionize' or 'eattach'
  real(8) :: moveprob(NMOVEMAX)     !< Probability of graph-move (not that these are normalized internally)
  integer :: gdsoutfreq !< GDS output frequency
  real(8) :: gdsthresh !< GDS graph-move attempt probability
  integer :: ngdsrelax !< Number of SD relaxation steps under graph restraints after graph moves
  real(8) :: gdsdtrelax !< Step size for SD relaxation after graph moves.
  logical :: optaftermove  !< Flag indicating whether to optimize on PES after each graph move.
  logical :: skiprepeats  !< Flag indicating whether to skip repeat graph moves.
  logical :: gatherreactivemol  !< Flag indicating whether to make the reactive molecles gather around into a "circle"
  logical :: simpleopt   !< optimises end images in NEB, all molecules in image at once...
  logical :: nebrestrend !< Flag indicating whether or not graph-restraint potenial is applied to end-points
                         !! during NEB optimization.
  logical :: removecollisions !< Flag indicating whether to refine path after move under purely repulsive
                              !! potential for non-bonded atoms.
  logical :: shimmybeads !< Flag determining if we want flat PES beads at the ends to shimmy (center the path )
  logical :: restartgds  !< Flag determining weather user wishes to read a ChemReacList file to restart
                         !! the calculation (xyz needed)
  integer :: nallowbonds   !< Number of defined bonding constraints for move-generated structures.
  character (len=2) :: allowbondsatom(NTYPEMAX,2) !< Atom labels for nallowbonds bonding constraints.
  integer :: allowbondsmax(NTYPEMAX)     !< Maximum number of bonds associated with nallowbonds constraints.
  logical :: forbidgraphs  !< Logical flag indicating whether to use the forbidden graphs file in allowing
  integer :: nforbid       !< Number of forbidden graph definitions.
  integer :: naforbid(NFORBIDMAX) !< Number of atoms in each forbidden graph.
  integer :: gforbid(NFORBIDMAX,NAMOVEMAX,NAMOVEMAX)     !< Forbidden graph definitions.
  character (len=4) :: forbidlabel(NFORBIDMAX,NAMOVEMAX)  !< Labels for atoms involved in forbidden graphs.
  integer :: maxmolcharge !< Maximum allowed charged on a single molecule.
  integer :: minmolcharge !< Minimum allowed charged on a single molecule.
  integer :: nchargemol  !< Maximum total charge at any reaction step.
  integer :: maxstepcharge   !< Maximum number of steps inolving charge changes.
  integer :: maxtotalcharge   !< Maximum number of steps inolving charge changes.
  integer :: projforcetype !< The type of projected Force used in NEB


  logical :: lmoldata !< logical flag to switch the instructions to print the molecular database on
  character (len=100) :: pathwaysfile   !< hold the file name for the pathway to the molecular database files (NEB+GDS DATA)

  ! Valence constraints:
  !
  integer :: nvalcon                         !< Number of atomic valence constraints for products
  character (len=2) :: valatom(NTYPEMAX)     !< Atom label corresponding to each set of valence constraints.
  integer :: valrange(NTYPEMAX,2)            !< Valence range for each valence constraint.
  logical :: valfz(NTYPEMAX)                 !< flag determining if this rule applies to frozen atoms bonding to non-frozen ones

  integer :: nrxval                          !< Number of atomic valence constraints for reactants.
  character (len=2) :: rxvalatom(NTYPEMAX,2)     !< Atom label corresponding to each set of reactant valence constraints.
  integer :: rxvalrange(NTYPEMAX,2)            !< Valence ranges for each reactant valence constraint.

  ! Molecular Valency contraints
  !
  logical              :: molecularvalency = .false.
  integer              :: nmv = 0                   !< the number of sets of molecular valence constraints
  character(2)         :: mvconel(30,20)            !< Elements in the molecular velency sets according to user
                                                    !< dimensions (setindex,atomic number)
  integer              :: mvconva(30,20)            !< Valency of element in in the molecular velency sets associated with mvconel
                                                    !< dimensions (setindex,valency assiciated with that element)
  integer              :: mvconne(30)               !< number of elements in each molecular valency set
  integer              :: mxvlel(30)                !< the element with the largers number of bonds allowed, used for
                                                    !< knowing how many times to loop for counting valences...see routine
  integer              :: ncharcon(30)              !< the maximun number of charged atoms (non-matched valency) this constrain allows for
                                                    !< ("charged" meaning radicals dont exist in this context, ie every atom has the full shell occupation)
  ! Molecular Valency contraints
  !
  logical              :: lewiscon = .false.
  integer              :: nlc = 0                   !< the number of sets of lewis constraints
  character(2)         :: lcconel(30,20)            !< Elements in the lewis constraint sets according to user
                                                    !< dimensions (setindex,atomic number)
  integer              :: lcconne(30)               !< number of elements in each molecular valency set
  integer              :: nlcchar(30)               !< the maximun number of charged atoms (non-matched valency) this constrain allows for
                                                    !< ("charged" meaning radicals dont exist in this context, ie every atom has the full shell occupation)
  logical              :: lcshift(30,20) = .false.  !< flag determinin if this element should be shifted out of the molecule before calculating constraints

  ! Molecular optimization parameters.
  !
  integer :: ngen   !< Number of GA generations.
  integer :: npop   !< Population for GA.
  integer :: ncross   !< Crossover events per generation.
  integer :: nmut  !< Number of mutation events per generation
  logical :: readcore
  integer :: nacore


  ! Checmial structure generation.
  !
  integer :: nheteromax, nheteromin  !< Max / min number of heteroatoms in generated structures.
  integer :: nelprob         !< Number of element probabilities defined in input
  character(len=2) :: elprob_label(NTYPEMAX)     !< Atomic label for element probability
  real(8) :: elprob(NTYPEMAX)       !< Appearance probability for each element in structure generation
  integer :: nmccxs        !< Maximum number of MC moves in structure generation
  real(8) :: mctemperature !< Temperature for MC structure generation.
  real(8) :: mcbondprob  !< Probability of bond creation in initial structure for MC optimization.
  real(8) :: nmolpenalty  !< Penalty factor for having more than 1 molecule in graph during MC generation.
  real(8) :: ringpenalty3     !< Penalty factor for having 3-rings in graph during MC generation.
  real(8) :: ringpenalty4     !< Penalty factor for having 4-rings in graph during MC generation.
  integer :: nheterolimit  !< Maximum number of heteroatoms in structures generated during GA.


  ! Pathfinder variables.
  !
  integer :: nrxn        !< Maximum number of reaction steps to locate in pathfinder calculation.
  integer :: nmcrxn      !< Number of MC annealing moves to attempt in initial graph search.
  integer :: nmcsearch   !< Number of MC moves during search phase.
  integer :: nmechmove   !< Maximum number of simultanerous changes in mechanism search.
  real(8) :: mcrxntemp   !< Temperature for MC reaction-string searching.
  integer :: igfunc      !< Graph function type.
  integer :: atomidx(3)  !< Atoms which are aligned before NEB.
  real(8) :: alphavbe    !< Scaling factor for bond-energy based contribution to graph-error function.

  !< Type definition for the parameteres used in the graphfinder algorithm
  type GraphPar
   integer              :: react(10), prod(10)  !< user supplied molecule ID from the Mol.dat file that will form the ends (products and reactants) of the
                                                                  !< reaction netowrk tree that they wish the code to calculate
   integer              :: nprod, nreact        !< the number of molecules of products and reactants the user has supplied
   integer              :: maxtreesolutions     !< the number of tree solutions the user requested to generate...
   integer              :: nni, nri             !< number of nodes and reactions the user want us to ignore during the graph generation algorithm...
   integer              :: maxcomb              !< the cutoff for the number of trees we keep from levels below the tree, to stop the combinatorial explosion
   integer              :: maxb                 !< the maximum number of nodes-with-more-than-one-reaction we are allowed to visit when traversing the tree network
   integer              :: mxevp                !< the maximum number of trees per eigenvector cluster the user wants to print out
   integer              :: treefav              !< -1: the tree should be unfavourable; 1: tree should be favourable;  0: it can be either (speaking about energetics)
   integer              :: maxisosteps          !< the maximum number of ismerization reactions per isomer group that can be found in the reaction trees
   double precision     :: rdiskstore           !< if a particuar branch has more than maxcomb*rdiskstore trees comming off it, those trees will be stored in disk, to avoid memory consumption
   double precision     :: mxndexpo             !< The order of magnitude of the maximum total number of nodes explored during the tree graph search prodedure.
   double precision     :: absthresh0         !< The initial sum of absolute energy of reactions threshold used to estimate the size of the network search... can be provided by the user
                                                !< if it has previously perfomed a calculation and already known a good number (printed in the logfile)
   integer, allocatable :: nignore(:)           !< an array of node ids to ignore when constructing the network
   integer, allocatable :: rignore(:)           !< an array of reaction ids to ignore when constructing the network
   logical              :: enforcefullrx        !< flag determining weather the user wants only trees that have all the reactants in the react array
   logical              :: readnds              !< reads the nds information from a restart binary file inputfilename_nds.bin if react, prod, nignore,rignore,Md are the same
   character(10)        :: threshmeth           !< weather to use absolute energies or lifetimes as the threshold of exploration of the network
  end type GraphPar
  type(GraphPar)        :: Gpar                 !< variable type of parameters used during graphfinder...

  ! matrix of molecules with atoms bonding/ breaking during a reaction
  ! if k1 = k2, id gives the id of atoms that molecule m1 shares to (if aid(id,1)) or takes from (if aid(id,2)) molecule m2 in the same image
  ! if k1 /= k2, id gives similar information but for inter-images
  ! ic gives information about what atoms of molecule m1 in image k are shared (are common, overlap) with molecule m2 of complementary image
  type molgm
     integer, allocatable :: id(:,:,:,:) ! the ID for naid and aid, given indecies id(k1,k2,m1,m2), where m1 = molecule id of image k1, m2 = molecule id of image k2
     integer, allocatable :: naid(:,:)   ! the number of atoms that these two molecules exchange naid(id,1:2)
     integer, allocatable :: aid(:,:,:)  ! the id of the atoms these two molecules exchange, now the last index aid(id,naid,1:2) =1 corresponds to the atom in molecule m1, =2 atom in molecule m2
     integer              :: nid         ! third dimension of id
     integer              :: na          ! fourth dimension of aid and dir
     integer, allocatable :: ic(:,:,:)   ! (k,m1,m2) first index: image,   second index, molecule m1 of image k, third, molecule m2 of complementary image to that of index k
     integer, allocatable :: naic(:)     ! k the image
     integer, allocatable :: aic(:,:)    ! k the image
     logical, allocatable :: arm(:,:,:)  ! what atoms are involved in molecular reaction for each image (true if those atoms form a bond)
  end type

  !> Type definition for the fingerprint of a reaction-path
  !
  type fingerprint
      sequence
      double precision, allocatable :: sprint(:,:,:)  !< sprint vector for beggining and end chemical structures (N elements)
      double precision, allocatable :: piv(:,:,:)     !< as above, but piv  (N*(N-1)/2 elements)
      double precision, allocatable :: menergy(:,:)   !< molecular energy per image, per active molecule i think
      double precision, allocatable :: mhash(:,:,:)   !< the molecular eigenvalues that characterises the graph
      integer                       :: graphmoveid    !< graph ID that was used on chemical structure IC
      integer                       :: ic             !< the start or end of the Path that was subject to the graph move
      integer                       :: nactmol(2)     !< number of active molecules in the start and end images
      character(100)                :: MakeBreakLab   !< label showing which elements are breaking/forming (symmetric)
      character(250), allocatable   :: Smiles(:,:)    !< an array of nactmol Smile strings in a unique order
      character(100), allocatable   :: mformula(:,:)  !< an array of nactmol molecular Formula strings in a unique order
      integer, allocatable          :: actnamol(:,:)  !< number of atoms per active molecules
      !integer, allocatable          :: mgraph(:,:,:)  !< vector of upper triangular graph matrix
      double precision              :: nbb            !< (number of bonds formed+broken)/2
      integer                       :: id             !< fingerprint ID
      type(fingerprint), pointer    :: next => null()
      type(fingerprint), pointer    :: prev => null()
  end type fingerprint

  integer                          :: moleculemax         !< Maximum number of fragments/molecules allowed to form ! seB
  ! Type definition for gds refreshing of templates
  !
  type refresh_template
   integer                          :: desirednmoves       !< Number of steps before reset to a template, 80% of the time !seb
   integer                          :: ntemplate           !< Number of defined templates seb]
   double precision                 :: DesiredNmoveProb    !< the probability that we reach the desirednmoves away from templates
   double precision , allocatable   :: templateprob(:)     !< the probability that
   character(len=50), allocatable   :: gdscxstemplate(:)   !< the filename of that template
  End type
  type(refresh_template) :: rtmpl
  type vtbashop
    logical                         :: MCaccepted      = .false.
    logical                         :: lbasinhopping = .false.
    double precision                :: T = 298.0d0
  end type
  type(vtbashop)                    :: bashop


end Module globaldata
