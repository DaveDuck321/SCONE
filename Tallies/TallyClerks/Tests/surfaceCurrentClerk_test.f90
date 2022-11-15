module surfaceCurrentClerk_test

    use numPrecision
    use dictionary_class,                only : dictionary
    use multiMap_class,                  only : multiMap
    use outputFile_class,                only : outputFile
    use particle_class,                  only : particle, particleState, P_NEUTRON
    use particleDungeon_class,           only : particleDungeon
    use scoreMemory_class,               only : scoreMemory
    use surfaceCurrentClerk_class,       only : surfaceCurrentClerk, SJResult
    use tallyResult_class,               only : tallyResult
    use testNeutronDatabase_class,       only : testNeutronDatabase
    use pFUnit_mod

    implicit none

  @testCase
    type, extends(TestCase) :: test_surfaceCurrentClerk
      private
      type(surfaceCurrentClerk) :: clerk
      type(multiMap)            :: map
      integer(shortInt)         :: N
    contains
      procedure :: setUp
      procedure :: tearDown

      !! (Private) test helper function
      procedure :: printSurfaceCurrents
    end type test_surfaceCurrentClerk

  contains

    !!
    !! Test helper: visualize the surface current
    !!
    subroutine printSurfaceCurrents(this, currentMatrix, inState)
      class(test_surfaceCurrentClerk), intent(inout) :: this
      class(SJResult), intent(in)                    :: currentMatrix
      type(particleState), intent(in)                :: inState
      type(particleState)                            :: state
      integer(shortInt) :: range
      integer(shortInt) :: idx, x, y, z
      real(defReal)     :: scale

      state = inState
      scale = 1.0
      range = this % n / 2

      print *, ""
      print *, "Surface current matrix"
      do x = -range, range
        do y = -range, range
          do z = -range, range
            state % r = [0.5 + x * scale, 0.5 + y * scale, 0.5 + z * scale]
            idx = this % map % map(state)
            if (idx /= 0) then
              print *, "r: ", state % r, "mapping: ", idx, "  Current: ", currentMatrix % JM (1, idx, 1)
            end if
          end do
        end do
      end do

    end subroutine printSurfaceCurrents


    !!
    !! Sets up test_surfaceCurrentClerk object we can use in a number of tests
    !!
    subroutine setUp(this)
      class(test_surfaceCurrentClerk), intent(inout) :: this
      type(dictionary)                                     :: dict
      type(dictionary)                                     :: multiMapDict, mapx, mapy, mapz
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

      call multiMapDict % init(5)
      call multiMapDict % store('type','multiMap')
      call multiMapDict % store('maps', ['mapx', 'mapy', 'mapz'])
      call multiMapDict % store('mapx', mapx)
      call multiMapDict % store('mapy', mapy)
      call multiMapDict % store('mapz', mapz)


      ! Build intput dictionary
      call dict % init(2)
      call dict % store('type','surfaceCurrentClerk')
      call dict % store('spacing', [1.0_defReal, 1.0_defReal, 1.0_defReal])
      call dict % store('spaceMap', multiMapDict)

      name = 'testClerk'
      call this % clerk % init(dict, name)
      call this % map % init(multiMapDict)
      this % N = 2

      call mapx % kill()
      call mapy % kill()
      call mapz % kill()
      call multiMapDict % kill()
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
      type(particleState)                                  :: state
      type(particleDungeon)                                :: pop
      type(testNeutronDatabase)                            :: xsData
      real(defReal)                                        :: val
      integer(shortInt)                                    :: i, idx
      class(tallyResult), allocatable                      :: res
      real(defReal), dimension(3), parameter :: DEFAULT_VELOCITY = [0.0, 0.0, 0.0]
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

      ! Crosses boundary at x=0.0 in +ve direction with weight 1.0
      p % w = 1.0
      p % preCollision % r = [-0.1_defReal, 0.1_defReal, 0.1_defReal]
      call p % coords % init([0.1_defReal, 0.1_defReal, 0.1_defReal], DEFAULT_VELOCITY)
      call this % clerk % reportInColl(p, xsData, mem)

      ! Crosses boundary at x=0.0 in -ve direction with weight 0.25
      p % w = 0.25
      p % preCollision % r = [0.1_defReal, 0.1_defReal, 0.1_defReal]
      call p % coords % init([-0.1_defReal, 0.1_defReal, 0.1_defReal], DEFAULT_VELOCITY)
      call this % clerk % reportInColl(p, xsData, mem)

      ! Close cycle
      call mem % closeCycle(ONE)

      ! Verify run-time results
      call this % clerk % getResult(res, mem)

      select type(res)
        class is (SJResult)
          @assertEqual(res % N, (this % N) ** 3)

          ! Test current corresponding to the x=0 boundary
          state = p
          state % r = [-0.5, 0.5, 0.5]
          idx = this % map % map(state)
          @assertEqual(0.75_defReal, res % JM(1, idx, 1), TOL)

          ! Confirm no other currents have been recorded
          do i=1, (this % N ** 3)
            if (i /= idx) then
              @assertEqual(0.0_defReal, res % JM(1, i, 1), TOL)
            end if
          end do
        class default
          @assertEqual(0, 1)
      end select

      ! Clean
      call xsData % kill()
      call pop % kill()

    end subroutine testSimpleUseCase

end module surfaceCurrentClerk_test
