module surfaceCurrentClerk_class

    use numPrecision
    use tallyCodes
    use genericProcedures,          only : fatalError
    use dictionary_class,           only : dictionary
    use particle_class,             only : particle, particleState
    use particleDungeon_class,      only : particleDungeon
    use outputFile_class,           only : outputFile

    ! Basic tally modules
    use scoreMemory_class,          only : scoreMemory
    use tallyClerk_inter,           only : tallyClerk, kill_super => kill
    use tallyResult_class,          only : tallyResult

    ! Nuclear Data
    use nuclearDatabase_inter,      only : nuclearDatabase
    use neutronMaterial_inter,      only : neutronMaterial, neutronMaterial_CptrCast

    ! Tally Maps
    use tallyMap_inter,             only : tallyMap
    use tallyMapFactory_func,       only : new_tallyMap

    implicit none
    private

    !!
    !! Surface current tally
    !!
    !! Private Members:
    !!   map      -> Map to segment the space
    !!   N        -> Number of Bins
    !!
    !! Interface:
    !!   tallyClerk Interface
    !!
    !! Sample dictionary input:
    !!
    !!  clerkName {
    !!      type surfaceCurrentClerk;
    !!      spaceMap {
    !!        type multiMap;
    !!        spacing 0.2;
    !!        maps (mapx mapy mapz);
    !!        mapx {
    !!          type spaceMap;
    !!          axis x;
    !!          grid lin;
    !!          min -10.0;
    !!          max 10.0;
    !!          N 10;
    !!        }
    !!        mapy {
    !!          type spaceMap;
    !!          axis y;
    !!          grid lin;
    !!          min -10.0;
    !!          max 10.0;
    !!          N 10;
    !!        }
    !!        mapz {
    !!          type spaceMap;
    !!          axis z;
    !!          grid lin;
    !!          min -10.0;
    !!          max 10.0;
    !!          N 10;
    !!        }
    !!      }
    !!    }
    !!  }
    !!
    type, public, extends(tallyClerk) :: surfaceCurrentClerk
      private
      !! Map defining the discretisation
      class(tallyMap), allocatable           :: map
      integer(shortInt)                      :: N = 0 !! Number of bins
      real(defReal)                          :: spacing
      !type(macroResponse)                    :: resp


    contains
      ! Procedures used during build
      procedure  :: init
      procedure  :: validReports
      procedure  :: getSize

      ! File reports and check status -> run-time procedures
      procedure  :: reportInColl
      procedure  :: reportCycleEnd

      ! Overwrite default run-time result procedure
      procedure  :: getResult

      ! Output procedures
      procedure  :: display
      procedure  :: print

      ! Deconstructor
      procedure  :: kill

      ! (Private) Utility
      procedure  :: reportCurrentInDirection
    end type surfaceCurrentClerk

    !!
    !! Surface current result class
    !!   Stored in column first order
    !!    dim1 -> 1 is J_{x+1/2, y, z}, 2 is J_{x, y+1/2, z}, 3 is J_{x, y, z+1/2}
    !!    dim2 -> 1 is values; 2 is STDs
    !!
    type,public, extends( tallyResult) :: SJResult
      integer(shortInt)                           :: N  = 0
      real(defReal), dimension(:,:,:),allocatable :: JM ! Current matrix
    end type SJResult

  contains

    !!
    !! Initialise clerk from dictionary and name
    !!
    !! See tallyClerk_inter for details
    !!
    subroutine init(self, dict, name)
      class(surfaceCurrentClerk), intent(inout) :: self
      class(dictionary), intent(in)                   :: dict
      character(nameLen), intent(in)                  :: name

      ! Assign name
      call self % setName(name)

      ! Read maps
      call new_tallyMap(self % map, dict % getDictPtr('spaceMap'))

      call dict % get(self % spacing, 'spacing')

      self % N = self % map % bins(0)

    end subroutine init

    !!
    !! Returns array of codes that represent diffrent reports
    !!
    !! See tallyClerk_inter for details
    !!
    function validReports(self) result(validCodes)
      class(surfaceCurrentClerk),intent(in) :: self
      integer(shortInt),dimension(:),allocatable  :: validCodes

      ! TODO: use trans_CODE here instead
      validCodes = [inColl_CODE, cycleEnd_Code]

    end function validReports

    !!
    !! Return memory size of the clerk
    !!
    !! See tallyClerk_inter for details
    !!
    elemental function getSize(self) result(S)
      class(surfaceCurrentClerk), intent(in) :: self
      integer(shortInt)                            :: S

      ! TODO: why not mult by 2 (value, std) and sizeof(Real)
      S = 3 * self % N

    end function getSize


    subroutine reportCurrentInDirection(self, mem, state_in, start, end, dir)
      class(surfaceCurrentClerk), intent(inout) :: self
      type(scoreMemory), intent(inout)          :: mem
      type(particleState), intent(in)           :: state_in
      real(defReal),dimension(3), intent(in)    :: start, end
      integer(shortInt), intent(in)             :: dir
      type(particleState)                       :: state
      integer(shortInt),dimension(3)            :: coordAtEnd, coordAtStart
      integer(longInt)                          :: addr, stepIdx
      real(defReal)                             :: currentContribution
      real(defReal),dimension(3)                :: step, diff
      integer(shortInt)                         :: i, offsetDueToSign, gridDirection

      ! Copy the state to allow mutation
      state = state_in

      ! Discretize onto the grid
      coordAtEnd = floor(end * (1.0/ (self % spacing)))
      coordAtStart = floor(start * (1.0/ (self % spacing)))

      diff = end - start
      currentContribution = (diff(dir)) / norm2(diff)

      if (diff(dir) > 0) then
        offsetDueToSign = 0
        gridDirection = 1
      else
        offsetDueToSign = -1
        gridDirection = -1
      end if

      ! Normalized so step(dir) == 1
      step = diff / diff(dir)

      ! Score the current across the final surface
      state % r (dir) = end(dir) + offsetDueToSign * step(dir) * (self % spacing)
      stepIdx = self % map % map(state)
      addr = self % getMemAddress() + ((dir - 1) * self % N) + stepIdx - 1
      call mem % score(currentContribution, addr)

      ! Score the current between all intermediate surfaces
      do i = 0, abs(coordAtEnd(dir) - coordAtStart(dir)) - 1, gridDirection
        state % r = start + (i + offsetDueToSign) * step * (self % spacing)
        stepIdx = self % map % map(state)

        addr = self % getMemAddress() + ((dir - 1) * self % N) + stepIdx - 1
        call mem % score(currentContribution, addr)
      end do

    end subroutine reportCurrentInDirection

    !!
    !! Process incoming collision report
    !!
    !! See tallyClerk_inter for details
    !!
    subroutine reportInColl(self, p, xsData, mem)
      class(surfaceCurrentClerk), intent(inout) :: self
      class(particle), intent(in)                     :: p
      class(nuclearDatabase),intent(inout)            :: xsData
      type(scoreMemory), intent(inout)                :: mem
      type(particleState)                             :: state
      integer(shortInt)                               :: sIdx, cIdx
      character(100), parameter :: Here = 'reportInColl (surfaceCurrentClerk_class.f90)'

      sIdx = self % map % map(p % preCollision)
      state = p
      cIdx = self % map % map(state)

      ! No boundaries crossed => no current
      if (sIdx == cIdx) return

      call self % reportCurrentInDirection(mem, state, p % preCollision % r, state % r, 1_shortInt)
      call self % reportCurrentInDirection(mem, state, p % preCollision % r, state % r, 2_shortInt)
      call self % reportCurrentInDirection(mem, state, p % preCollision % r, state % r, 3_shortInt)
    end subroutine reportInColl

    !!
    !! Process cycle end
    !!
    !! See tallyClerk_inter for details
    !!
    subroutine reportCycleEnd(self, end, mem)
      class(surfaceCurrentClerk), intent(inout) :: self
      class(particleDungeon), intent(in)              :: end
      type(scoreMemory), intent(inout)                :: mem


      if(mem % lastCycle()) then
        ! TODO: why normalize?
      end if

    end subroutine reportCycleEnd

    !!
    !! Return result from the clerk for interaction with Physics Package
    !! Returns SJResult defined in this module
    !! If res is already allocated to a JM of fitting size it reuses already allocated space
    !! This should improve performance when updating estimate of JM each cycle
    !!
    !! See tallyClerk_inter for details
    !!
    pure subroutine getResult(self, res, mem)
      class(surfaceCurrentClerk), intent(in)  :: self
      class(tallyResult),allocatable, intent(inout) :: res
      type(scoreMemory), intent(in)                 :: mem
      integer(shortInt)                             :: i, j
      integer(longInt)                              :: addr
      real(defReal)                                 :: val, STD

      ! Allocate result to SJResult
      ! Do not deallocate if already allocated to SJResult
      ! Its not too nice -> clean up
      if(allocated(res)) then
        select type(res)
          class is (SJResult)
            ! Do nothing

          class default ! Reallocate
            deallocate(res)
            allocate( SJResult :: res)
       end select

      else
        allocate( SJResult :: res)

      end if

      ! Load data into the JM
      select type(res)
        class is(SJResult)
          ! Check size and reallocate space if needed
          if (allocated(res % JM)) then
            if( any(shape(res % JM) /= [3, self % N, 2])) then
              deallocate(res % JM)
              allocate(res % JM(3, self % N, 2))
            end if
          else
            allocate(res % JM(3, self % N, 2))
          end if

          ! Set size of the JM
          res % N = self % N

          ! Load entries
          addr = self % getMemAddress() - 1
          do i = 1, 3
            do j = 1, self % N
              addr = addr + 1
              call mem % getResult(val, STD, addr)
              res % JM(i, j, 1) = val
              res % JM(i, j, 2) = STD
            end do
          end do

      end select
    end subroutine getResult

    !!
    !! Display convergence progress on the console
    !!
    !! See tallyClerk_inter for details
    !!
    subroutine display(self, mem)
      class(surfaceCurrentClerk), intent(in) :: self
      type(scoreMemory), intent(in)                :: mem

      print *, 'surfaceCurrentClerk does not support display yet'

    end subroutine display

    !!
    !! Write contents of the clerk to output file
    !!
    !! See tallyClerk_inter for details
    !!
    subroutine print(self, outFile, mem)
      class(surfaceCurrentClerk), intent(in) :: self
      class(outputFile), intent(inout)             :: outFile
      type(scoreMemory), intent(in)                :: mem
      integer(shortInt)                            :: i
      integer(longInt)                             :: addr
      real(defReal)                                :: val, std
      character(nameLen)                           :: name

      ! Begin block
      call outFile % startBlock(self % getName())

      ! Print map information
      call self % map % print(outFile)

      ! Print surface current matrix
      name = 'JM'
      addr = self % getMemAddress() - 1

      call outFile % startArray(name, [3, self % N])

      do i = 1, 3 * self % N
        addr = addr + 1
        call mem % getResult(val, std, addr)
        call outFile % addResult(val, std)
      end do
      call outFile % endArray()

      call outFile % endBlock()

    end subroutine print

    !!
    !! Returns to uninitialised state
    !!
    !! See tallyClerk_inter for details
    !!
    elemental subroutine kill(self)
      class(surfaceCurrentClerk), intent(inout) :: self

      ! Call superclass
      call kill_super(self)

      if(allocated(self % map)) deallocate(self % map)
      self % N = 0
    end subroutine kill

end module surfaceCurrentClerk_class
