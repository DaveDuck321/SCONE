module surfaceCurrentClerk_test

    use numPrecision
    use tallyResult_class,               only : tallyResult
    use surfaceCurrentClerk_class, only : surfaceCurrentClerk, SJResult
    use particle_class,                  only : particle, particleState, P_NEUTRON
    use particleDungeon_class,           only : particleDungeon
    use dictionary_class,                only : dictionary
    use scoreMemory_class,               only : scoreMemory
    use testNeutronDatabase_class,       only : testNeutronDatabase
    use outputFile_class,                only : outputFile
    use pFUnit_mod

    implicit none

  @testCase
    type, extends(TestCase) :: test_surfaceCurrentClerk
      private
      type(surfaceCurrentClerk) :: clerk
      integer(shortInt) :: N
    contains
      procedure :: setUp
      procedure :: tearDown
    end type test_surfaceCurrentClerk

  contains

    !!
    !! Test helper: turns a coordinate set into a multimap index
    !!
    integer(shortInt) function getMultiMapIndex(coord, N)
    implicit none
    integer(shortInt), dimension(3), intent(in) :: coord
    integer(shortInt), intent(in) :: N

    getMultiMapIndex = coord(1) + (coord(2) * N) + (coord(3) * (N**2))

    end function getMultiMapIndex

    !!
    !! Sets up test_surfaceCurrentClerk object we can use in a number of tests
    !!
    subroutine setUp(this)
      class(test_surfaceCurrentClerk), intent(inout) :: this
      type(dictionary)                                     :: dict
      type(dictionary)                                     :: multiMap, mapx, mapy, mapz
      character(nameLen)                                   :: name

      call mapx % init(6)
      call mapx % store('type', 'spaceMap')
      call mapx % store('axis', 'x')
      call mapx % store('grid', 'lin')
      call mapx % store('min', -1.0_defReal)
      call mapx % store('max', 1.0_defReal)
      call mapx % store('N', 2)

      call mapy % init(6)
      call mapy % store('type', 'spaceMap')
      call mapy % store('axis', 'y')
      call mapy % store('grid', 'lin')
      call mapy % store('min', -1.0_defReal)
      call mapy % store('max', 1.0_defReal)
      call mapy % store('N', 2)

      call mapz % init(6)
      call mapz % store('type', 'spaceMap')
      call mapz % store('axis', 'z')
      call mapz % store('grid', 'lin')
      call mapz % store('min', -1.0_defReal)
      call mapz % store('max', 1.0_defReal)
      call mapz % store('N', 2)

      call multiMap % init(5)
      call multiMap % store('type','multiMap')
      call multiMap % store('maps', ['mapx', 'mapy', 'mapz'])
      call multiMap % store('mapx', mapx)
      call multiMap % store('mapy', mapy)
      call multiMap % store('mapz', mapz)


      ! Build intput dictionary
      call dict % init(2)
      call dict % store('type','surfaceCurrentClerk')
      call dict % store('spacing', 1.0_defReal)
      call dict % store('spaceMap', multiMap)

      name = 'testClerk'
      call this % clerk % init(dict, name)
      this % N = 2

      call mapx % kill()
      call mapy % kill()
      call mapz % kill()
      call multiMap % kill()
      call dict % kill()
    end subroutine setUp

    !!
    !! Kills test_surfaceCurrentClerk object we can use in a number of tests
    !!
    subroutine tearDown(this)
      class(test_surfaceCurrentClerk), intent(inout) :: this

      call this % clerk % kill()

    end subroutine tearDown

  !!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
  !! PROPER TESTS BEGIN HERE
  !!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>


    !!
    !! Test correctness in a simple use case
    !!
  @Test
    subroutine testSimpleUseCase(this)
      class(test_surfaceCurrentClerk), intent(inout) :: this
      type(scoreMemory)                                    :: mem
      type(particle)                                       :: p
      type(particleState)                                  :: phase
      type(particleDungeon)                                :: pop
      type(testNeutronDatabase)                            :: xsData
      real(defReal)                                        :: val
      class(tallyResult), allocatable                      :: res
      integer(shortInt) :: i, idx
      real(defReal), parameter :: TOL = 1.0E-7


      ! Create score memory
      call mem % init(int(this % clerk % getSize(), longInt) , 1, batchSize = 1)
      call this % clerk % setMemAddress(1_longInt)

      ! Create test transport Nuclear Data
      call xsData % build(1.1_defReal, fissionXS = 1.1_defReal, nuFissionXS = 2.0_defReal)

      ! Create one particle that can be made to collide
      ! repeatedly in several materials

      ! Score some events
      p % type = P_NEUTRON
      p % w = 1.0

      ! Crosses boundary at 0.0 (with 0 velocity)
      p % preCollision % r = [-0.1_defReal, 0.1_defReal, 0.1_defReal]
      call p % coords % init([0.1_defReal, 0.1_defReal, 0.1_defReal], [0.0_defReal, 0.0_defReal, 0.0_defReal])

      call this % clerk % reportInColl(p, xsData, mem)

      ! Close cycle
      call mem % closeCycle(ONE)

      ! Verify results (TODO)
      ! print *, val
      ! do i=0, 3*2*2*2
      !   call mem % getResult(val, 1_longInt + i)
      !   print *, val
      ! end do
      ! @assertEqual(0.22580645161_defReal, val, TOL)

      ! Verify run-time result
      call this % clerk % getResult(res, mem)

      select type(res)
        class is (SJResult)
          @assertEqual(res % N, (this % N) ** 3)

          idx = 1 + getMultiMapIndex([1, 1, 1], this % N)
          ! print *, ""
          ! do i=1, 2*2*2
          !   print *, res % JM(1, i, 1)
          ! end do
          @assertEqual(1.0_defReal, res % JM(1, idx, 1), TOL)
        class default
          @assertEqual(0, 1)
      end select

      ! Clean
      call xsData % kill()
      call pop % kill()

    end subroutine testSimpleUseCase

  end module surfaceCurrentClerk_test
