module infSurf_class

  use numPrecision
  use universalVariables
  use genericProcedures, only : fatalError
  use vector_class,      only : vector
  use dictionary_class,  only : dictionary
  use surface_inter,     only : surface, printSurfDef, surfaceShelf

  implicit none
  private

  !!
  !! Constants describing surface properties
  !!
  character(nameLen),parameter :: TYPE_NAME    = 'infSurf'

  !!
  !! Constructor
  !!
  interface infSurf
    module procedure infSurf_fromDict
  end interface

  !!
  !! Infinite surface for homogeneous regions
  !!
  type, public, extends (surface) :: infSurf
    private
  contains
    ! Initialisation & Identification procedures
    procedure :: init
    procedure :: type
    procedure :: getDef
    procedure :: boundingBox

    ! Runtime procedures
    procedure :: evaluate
    procedure :: distance
    procedure :: normalVector

  end type infSurf

contains

  !!
  !! Create an infinity surface object
  !!
  subroutine init(self, id)
    class(infSurf), intent(inout)   :: self
    integer(shortInt), intent(in)   :: id

    call self % setId(id)

    ! Hash and store surface definition
    call self % hashSurfDef()

  end subroutine init

  !!
  !! Returns an initialised instance of infSurf from dictionary and name
  !!
  function infSurf_fromDict(dict) result(new)
    class(dictionary), intent(in)  :: dict
    type(infSurf)                  :: new
    integer(shortInt)              :: id
    character(100), parameter :: Here ='infSurf_fromDict (infSurf_class.f90)'

    call dict % get(id, 'id')
    if(id < 1) call fatalError(Here,'Invalid surface id provided')

    call new % init(id)

  end function infSurf_fromDict

  !!
  !! Always return inside for an infinite surface
  !!
  elemental subroutine evaluate(self, res, r)
    class(infSurf), intent(in)              :: self
    real(defReal), intent(out)              :: res
    type(vector), intent(in)                :: r

    res = -ONE

  end subroutine evaluate

  !!
  !! Return parameter character containing TYPE NAME
  !!
  elemental function type(self)
    class(infSurf), intent(in) :: self
    character(nameLen)         :: type

    type = TYPE_NAME

  end function type

  !!
  !! Return string with definition of this surface
  !!
  pure subroutine getDef(self,string)
    class(infSurf), intent(in)             :: self
    character(:),allocatable,intent(inout) :: string

    string = printSurfDef(TYPE_NAME, [ZERO])

  end subroutine getDef

  !!
  !! Returns an axis alligned bouding box of surface -ve halfspace
  !!
  pure subroutine boundingBox(self,origin, halfwidth)
    class(infSurf), intent(in)              :: self
    real(defReal), dimension(3),intent(out) :: origin
    real(defReal), dimension(3),intent(out) :: halfwidth

    origin    = ZERO
    halfwidth = INFINITY

  end subroutine boundingBox

  !!
  !! Calculate distance to infinity's surface - always infinity
  !!
  elemental subroutine distance(self, dist, idx, r, u)
    class(infSurf), intent(in)     :: self
    real(defReal), intent(out)     :: dist
    integer(shortInt), intent(out) :: idx
    type(vector), intent(in)       :: r
    type(vector), intent(in)       :: u

    ! Set index
    idx = self % myIdx()

    dist = INFINITY

  end subroutine distance

  !!
  !! Supply the normal vector
  !!
  elemental function normalVector(self, r) result(normal)
    class(infSurf), intent(in)       :: self
    type(vector), intent(in)         :: r
    type(vector)                     :: normal

    ! Interpret infinate surface as a shpere centered at (0,0,0) with infinate radius
    ! Normal vector is just normalised current position
    normal = r / r % L2norm()

  end function normalVector

end module infSurf_class