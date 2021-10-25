
!
!***************************************************************************************
!
!> @brief Functions
!
!> Library of Generic Function
!
!
!***************************************************************************************
!

Module functions
  use constants

contains

!***********************************************************************
!> is_numeric
!!
!! checks to see if input is number
!!
!********************************************************************
  function is_numeric(string)
    implicit none
    character(len=1), intent(in) :: string
    logical :: is_numeric
    is_numeric = .false.
    if (string.eq.'0') then
       is_numeric = .true.
    else if (string.eq.'1') then
       is_numeric = .true.
    else if (string.eq.'2') then
       is_numeric = .true.
    else if (string.eq.'3') then
       is_numeric = .true.
    else if (string.eq.'4') then
       is_numeric = .true.
    else if (string.eq.'5') then
       is_numeric = .true.
    else if (string.eq.'6') then
       is_numeric = .true.
    else if (string.eq.'7') then
       is_numeric = .true.
    else if (string.eq.'8') then
       is_numeric = .true.
    else if (string.eq.'9') then
       is_numeric = .true.
    endif
  end function is_numeric


!
!************************************************************************
!> xproduct
!!
!! cross product between 3d vectors, not for general space. returns vector
!!
!! - a,b: vectors
!!
!************************************************************************
!
  function xproduct(a, b)
    implicit none
    double precision :: xproduct(3), a(3), b(3)
    xproduct(1) = a(2)*b(3) - a(3)*b(2)
    xproduct(2) = a(3)*b(1) - a(1)*b(3)
    xproduct(3) = a(1)*b(2) - a(2)*b(1)
    return
  end function xproduct


  !
  !************************************************************************
  !> LabelToNumber
  !!
  !! Converts an input atomic label into a number.
  !!
  !! - label: the input atomic label.
  !!
  !************************************************************************
  !
  Function LabelToNumber(label)
    implicit none
    integer :: LabelToNumber
    character (len = 2) :: label

    select case (trim(label))

    case('H')
       LabelToNumber = 1

    case('He')
       LabelToNumber = 2

    case('Li')
       LabelToNumber = 3

    case('Be')
       LabelToNumber = 4

    case('B')
       LabelToNumber = 5

    case('C')
       LabelToNumber = 6

    case('N')
       LabelToNumber = 7

    case('O')
       LabelToNumber = 8

    case('F')
       LabelToNumber = 9

    case('Ne')
       LabelToNumber = 10

    case('Na')
       LabelToNumber = 11

    case('Mg')
       LabelToNumber = 12

    case('Al')
       LabelToNumber = 13

    case('Si')
       LabelToNumber = 14

    case('P')
       LabelToNumber = 15

    case('S')
       LabelToNumber = 16

    case('Cl')
       LabelToNumber = 17

    case('Ar')
       LabelToNumber = 18

    case('K')
       LabelToNumber = 19

    case('Ca')
       LabelToNumber = 20

    case('Sc')
       LabelToNumber = 21

    case('Ti')
       LabelToNumber = 22

    case('V')
       LabelToNumber = 23

    case('Cr')
       LabelToNumber = 24

    case('Mn')
       LabelToNumber = 25

    case('Fe')
       LabelToNumber = 26

    case('Co')
       LabelToNumber = 27

    case('Ni')
       LabelToNumber = 28

    case('Cu')
       LabelToNumber = 29

    case('Zn')
       LabelToNumber = 30

    case('Ga')
       LabelToNumber = 31

    case('Ge')
       LabelToNumber = 32

    case('As')
       LabelToNumber = 33

    case('Se')
       LabelToNumber = 34

    case('Br')
       LabelToNumber = 35

    case('Kr')
       LabelToNumber = 36

    case('Rb')
       LabelToNumber = 37

    case('Sr')
       LabelToNumber = 38

    case('Y')
       LabelToNumber = 39

    case('Zr')
       LabelToNumber = 40

    case('Nb')
       LabelToNumber = 41

    case('Mo')
       LabelToNumber = 42

    case('Tc')
       LabelToNumber = 43

    case('Ru')
       LabelToNumber = 44

    case('Rh')
       LabelToNumber = 45

    case('Pd')
       LabelToNumber = 46

    case('Ag')
       LabelToNumber = 47

    case('Cd')
       LabelToNumber = 48

    case('In')
       LabelToNumber = 49

    case('Sn')
       LabelToNumber = 50

    case('Sb')
       LabelToNumber = 51

    case('Te')
       LabelToNumber = 52

    case('I')
       LabelToNumber = 53

    case('Xe')
       LabelToNumber = 54

    case('Cs')
       LabelToNumber = 55

    case('Ba')
       LabelToNumber = 56

    case('La')
       LabelToNumber = 57

    case('Ce')
       LabelToNumber = 58

    case('Pr')
       LabelToNumber = 59

    case('Nd')
       LabelToNumber = 60

    case('Pm')
       LabelToNumber = 61

    case('Sm')
       LabelToNumber = 62

    case('Eu')
       LabelToNumber = 63

    case('Gd')
       LabelToNumber = 64

    case('Tb')
       LabelToNumber = 65

    case('Dy')
       LabelToNumber = 66

    case('Ho')
       LabelToNumber = 67

    case('Er')
       LabelToNumber = 68

    case('Tm')
       LabelToNumber = 69

    case('Yb')
       LabelToNumber = 70

    case('Lu')
       LabelToNumber = 71

    case('Hf')
       LabelToNumber = 72

    case('Ta')
       LabelToNumber = 73

    case('W')
       LabelToNumber = 74

    case('Re')
       LabelToNumber = 75

    case('Os')
       LabelToNumber = 76

    case('Ir')
       LabelToNumber = 77

    case('Pt')
       LabelToNumber = 78

    case('Au')
       LabelToNumber = 79

    case('Hg')
       LabelToNumber = 80

    case('Tl')
       LabelToNumber = 81

    case('Pb')
       LabelToNumber = 82

    case('Bi')
       LabelToNumber = 83

    case('Po')
       LabelToNumber = 84

    case('At')
       LabelToNumber = 85

    case('Rn')
       LabelToNumber = 86

    case('G1')
       LabelToNumber = 87

    case default
       stop '* Error in LabelToNumber in structure.f90'

    end select

    return
  end Function LabelToNumber

  !
  !************************************************************************
  !> NumberToLabel
  !!
  !! Converts an input atomic number into label.
  !!
  !! - num: the input atomic number.
  !!
  !************************************************************************
  !
  Function NumberToLabel(num)
    implicit none
    integer :: num
    character (len = 2) :: NumberToLabel

    select case (num)

    case(1)
       NumberToLabel = 'H'

    case(2)
       NumberToLabel = 'He'

    case(3)
       NumberToLabel = 'Li'

    case(4)
       NumberToLabel = 'Be'

    case(5)
       NumberToLabel = 'B'

    case(6)
       NumberToLabel = 'C'

    case(7)
       NumberToLabel = 'N'

    case(8)
       NumberToLabel = 'O'

    case(9)
       NumberToLabel = 'F'

    case(10)
       NumberToLabel = 'Ne'

    case(11)
       NumberToLabel = 'Na'

    case(12)
       NumberToLabel = 'Mg'

    case(13)
       NumberToLabel = 'Al'

    case(14)
       NumberToLabel = 'Si'

    case(15)
       NumberToLabel = 'P'

    case(16)
       NumberToLabel = 'S'

    case(17)
       NumberToLabel = 'Cl'

    case(18)
       NumberToLabel = 'Ar'

    case(19)
       NumberToLabel = 'K'

    case(20)
       NumberToLabel = 'Ca'

    case(21)
       NumberToLabel = 'Sc'

    case(22)
       NumberToLabel = 'Ti'

    case(23)
       NumberToLabel = 'V'

    case(24)
       NumberToLabel = 'Cr'

    case(25)
       NumberToLabel = 'Mn'

    case(26)
       NumberToLabel = 'Fe'

    case(27)
       NumberToLabel = 'Co'

    case(28)
       NumberToLabel = 'Ni'

    case(29)
       NumberToLabel = 'Cu'

    case(30)
       NumberToLabel = 'Zn'

    case(31)
       NumberToLabel = 'Ga'

    case(32)
       NumberToLabel = 'Ge'

    case(33)
       NumberToLabel = 'As'

    case(34)
       NumberToLabel = 'Se'

    case(35)
       NumberToLabel = 'Br'

    case(36)
       NumberToLabel = 'Kr'

    case(37)
       NumberToLabel = 'Rb'

    case(38)
       NumberToLabel = 'Sr'

    case(39)
       NumberToLabel = 'Y'

    case(40)
       NumberToLabel = 'Zr'

    case(41 )
       NumberToLabel = 'Nb'

    case(42)
       NumberToLabel = 'Mo'

    case(43)
       NumberToLabel = 'Tc'

    case(44)
       NumberToLabel = 'Ru'

    case(45)
       NumberToLabel = 'Rh'

    case(46)
       NumberToLabel = 'Pd'

    case(47)
       NumberToLabel = 'Ag'

    case(48)
       NumberToLabel = 'Cd'

    case(49)
       NumberToLabel = 'In'

    case(50)
       NumberToLabel = 'Sn'

    case(51)
       NumberToLabel = 'Sb'

    case(52)
       NumberToLabel = 'Te'

    case(53)
       NumberToLabel = 'I'

    case(54)
       NumberToLabel = 'Xe'

    case(55)
       NumberToLabel = 'Cs'

    case(56)
       NumberToLabel = 'Ba'

    case(57)
       NumberToLabel = 'La'

    case(58)
       NumberToLabel = 'Ce'

    case(59)
       NumberToLabel = 'Pr'

    case(60)
       NumberToLabel = 'Nd'

    case(61)
       NumberToLabel = 'Pm'

    case(62)
       NumberToLabel = 'Sm'

    case(63)
       NumberToLabel = 'Eu'

    case(64)
       NumberToLabel = 'Gd'

    case(65)
       NumberToLabel = 'Tb'

    case(66)
       NumberToLabel = 'Dy'

    case(67)
       NumberToLabel = 'Ho'

    case(68)
       NumberToLabel = 'Er'

    case(69)
       NumberToLabel = 'Tm'

    case(70)
       NumberToLabel = 'Yb'

    case(71)
       NumberToLabel = 'Lu'

    case(72)
       NumberToLabel = 'Hf'

    case(73)
       NumberToLabel = 'Ta'

    case(74)
       NumberToLabel = 'W'

    case(75)
       NumberToLabel = 'Re'

    case(76)
       NumberToLabel = 'Os'

    case(77)
       NumberToLabel = 'Ir'

    case(78)
       NumberToLabel = 'Pt'

    case(79)
       NumberToLabel = 'Au'

    case(80)
       NumberToLabel = 'Hg'

    case(81)
       NumberToLabel = 'Tl'

    case(82)
       NumberToLabel = 'Pb'

    case(83)
       NumberToLabel = 'Bi'

    case(84)
       NumberToLabel = 'Po'

    case(85)
       NumberToLabel = 'At'

    case(86)
       NumberToLabel = 'Rn'

    case(87)
       NumberToLabel = 'G1'

    case default
       stop '* Error in NumberToLabel in structure.f90'

    end select

    return
  end Function NumberToLabel


!
!==================================================================================
!
! Initialize random-number seed from user input.
!
! irun - returned Random number seed.
!
!==================================================================================
!

integer Function SetRanSeed( irun )
  implicit none
  integer :: n, irun, idum
  integer,allocatable :: seed(:)

  call random_seed(size=n)
  allocate(seed(n))
  seed = irun
  call random_seed(put=seed)
  deallocate(seed)
  SetRanSeed = 1

end Function SetRanSeed


  !
  !************************************************************************
  !> MassFromLabels
  !!
  !! gets the mass array from labels
  !!
  !! - labels : label array
  !!
  !************************************************************************
  !
  Function MassFromLabels(labels) result(mmass)
    implicit none
    integer                             :: na, i
    character(2)                        :: labels(:)
    double precision, allocatable       :: mmass(:)
    na = size(labels)
    allocate(mmass(na))
    do i = 1, na
      mmass(i) = MASS(LabelToNumber(labels(i)))
    enddo
    return
  end function MassFromLabels

End Module functions
