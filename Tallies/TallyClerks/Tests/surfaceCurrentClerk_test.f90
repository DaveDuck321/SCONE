module surfaceCurrentClerk_test

   use numPrecision
   use dictionary_class,                only : dictionary
   use outputFile_class,                only : outputFile
   use spaceMap_class,                  only : spaceMap
   use particle_class,                  only : particle, particleState, P_NEUTRON
   use particleDungeon_class,           only : particleDungeon
   use scoreMemory_class,               only : scoreMemory
   use surfaceCurrentClerk_class,       only : surfaceCurrentClerk, SJResult
   use tallyResult_class,               only : tallyResult
   use testNeutronDatabase_class,       only : testNeutronDatabase
   use pFUnit_mod

   use tallyMap_inter,             only : tallyMap
   use tallyMapFactory_func,       only : new_tallyMap

   implicit none

   @testCase
   type, extends(TestCase) :: test_surfaceCurrentClerk
      private
      type(surfaceCurrentClerk) :: clerk
      type(spaceMap)            :: mapX
      type(spaceMap)            :: mapY
      type(spaceMap)            :: mapZ
      integer(shortInt)         :: NSpace
      integer(shortInt)         :: NEnergy

   contains
      procedure :: setUp
      procedure :: tearDown

      !! (Private) test helper function
      procedure :: printSurfaceCurrents
      procedure :: getMapIndex
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
      range = this % NSpace / 2

      print *, ""
      print *, "Surface current matrix"
      do x = -range, range
         do y = -range, range
            do z = -range, range
               state % r = [0.5 + x * scale, 0.5 + y * scale, 0.5 + z * scale]
               idx = this % getMapIndex(state)
               if (idx /= 0) then
                  print *, "r: ", state % r, "mapping: ", idx, "  Current: ", currentMatrix % JM (1, 1, idx, 1)
               end if
            end do
         end do
      end do

   end subroutine printSurfaceCurrents


   function getMapIndex(this, state) result(index)
      class(test_surfaceCurrentClerk), intent(inout) :: this
      type(particleState), intent(in)                :: state
      integer(shortInt)                              :: index


      index =  this % mapX % map(state)
      index = index + (this % mapX % bins(0) + 1) * this % mapY % map(state)
      index = index + ((this % mapX % bins(0) + 1) * (this % mapY % bins(0) + 1)) * this % mapZ % map(state)
   end function getMapIndex
   !!
   !! Sets up test_surfaceCurrentClerk object we can use in a number of tests
   !!
   subroutine setUp(this)
      class(test_surfaceCurrentClerk), intent(inout) :: this
      type(dictionary)                                     :: dict
      type(dictionary)                                     :: mapx, mapy, mapz, energyMap
      character(nameLen)                                   :: name


      call energyMap % init(5)
      call energyMap % store('type', 'energyMap')
      call energyMap % store('grid', 'lin')
      call energyMap % store('min', 1.0E-9_defReal)
      call energyMap % store('max', 20.0_defReal)
      call energyMap % store('N', 2)

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

      ! Build intput dictionary
      call dict % init(4)
      call dict % store('type','surfaceCurrentClerk')
      call dict % store('spacing', [1.0_defReal, 1.0_defReal, 1.0_defReal])
      call dict % store('spaceMapX', mapx)
      call dict % store('spaceMapY', mapy)
      call dict % store('spaceMapZ', mapz)
      call dict % store('energyMap', energyMap)

      name = 'testClerk'
      call this % clerk % init(dict, name)
      call this % mapX % init(mapx)
      call this % mapY % init(mapy)
      call this % mapZ % init(mapz)
      this % NSpace = 3  ! 2 cells, 3 boundaries
      this % NEnergy = 2

      call mapx % kill()
      call mapy % kill()
      call mapz % kill()
      call energyMap % kill()
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
      p % E = 0.01_defReal
      p % isMG = .false.

      ! Crosses boundary at x=0.0 in +ve direction with weight 1.0
      p % w = 1.0
      p % preCollision % r = [-0.1_defReal, 0.1_defReal, 0.1_defReal]
      call p % coords % init([0.1_defReal, 0.1_defReal, 0.1_defReal], DEFAULT_VELOCITY)
      call this % clerk % reportTrans(p, xsData, mem)

      ! Crosses boundary at x=0.0 in -ve direction with weight 0.25
      p % w = 0.25
      p % preCollision % r = [0.1_defReal, 0.1_defReal, 0.1_defReal]
      call p % coords % init([-0.1_defReal, 0.1_defReal, 0.1_defReal], DEFAULT_VELOCITY)
      call this % clerk % reportTrans(p, xsData, mem)

      ! Close cycle
      call mem % closeCycle(ONE)

      ! Verify run-time results
      call this % clerk % getResult(res, mem)

      select type(res)
       class is (SJResult)
         @assertEqual(res % NSpace, (this % NSpace) ** 3)
         @assertEqual(res % NEnergy, (this % NEnergy))

         state = p

         ! Test current corresponding to the x=0 boundary
         state = p
         state % r = [-0.5, 0.5, 0.5]
         idx = this % getMapIndex(state)

         @assertEqual(0.75_defReal, res % JM(1, 1, idx + 1, 1), TOL)

         ! Confirm no other currents have been recorded
         do i=1, (this % NEnergy) * ((this % NSpace) ** 3)
            if (i /= idx + 1) then
               @assertEqual(0.0_defReal, res % JM(1, 1, i, 1), TOL)
            end if
         end do
       class default
         @assertEqual(0, 1)
      end select

      ! Clean
      call xsData % kill()
      call pop % kill()

   end subroutine testSimpleUseCase


   !!
   !! Test correctness for boundary conditions
   !!
   @Test
   subroutine testBoundaryConditions(this)
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
      p % E = 0.01_defReal
      p % isMG = .false.

      ! Crosses boundary at x=-1.0 in -ve direction with weight 1.0
      p % w = 1.0
      p % preCollision % r = [-0.9_defReal, 0.1_defReal, 0.1_defReal]
      call p % coords % init([-1.1_defReal, 0.1_defReal, 0.1_defReal], DEFAULT_VELOCITY)
      call this % clerk % reportTrans(p, xsData, mem)

      ! Crosses boundary at x=1.0 in +ve direction with weight 0.5
      p % w = 0.5
      p % preCollision % r = [0.9_defReal, 0.1_defReal, 0.1_defReal]
      call p % coords % init([1.1_defReal, 0.1_defReal, 0.1_defReal], DEFAULT_VELOCITY)
      call this % clerk % reportTrans(p, xsData, mem)

      ! Close cycle
      call mem % closeCycle(ONE)

      ! Verify run-time results
      call this % clerk % getResult(res, mem)

      select type(res)
       class is (SJResult)
         state = p

         ! Test current corresponding to the x=-1 boundary
         state % r = [-0.5, 0.5, 0.5]
         idx = this % getMapIndex(state) - 1  ! x-stride is 1... index-1 is the boundary
         @assertEqual(-1.0_defReal, res % JM(1, 1, idx + 1, 1), TOL)

         state % r = [0.5, 0.5, 0.5]
         idx = this % getMapIndex(state)  ! each cell contains the current of the rhs boundary
         @assertEqual(0.5_defReal, res % JM(1, 1, idx + 1, 1), TOL)
       class default
         @assertEqual(0, 1)
      end select

      ! Clean
      call xsData % kill()
      call pop % kill()

   end subroutine testBoundaryConditions

end module surfaceCurrentClerk_test
