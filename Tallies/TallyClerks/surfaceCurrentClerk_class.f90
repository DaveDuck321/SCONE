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
      integer(shortInt)                               :: sIdx, cIdx, stepIdx, stepPreviousIdx
      integer(longInt)                                :: addr
      integer(shortInt),dimension(3)                  :: coordAtStart, coordAtEnd
      real(defReal),dimension(3)                      :: step, diff, currentContribution
      integer(shortInt)                               :: i
      class(neutronMaterial), pointer                 :: mat
      character(100), parameter :: Here = 'reportInColl (surfaceCurrentClerk_class.f90)'

      ! Find starting index in the map
      ! It is important that preCollision is not changed by a collisionProcessor
      !before the particle is fed to the tally, otherwise results will be meaningless
      sIdx = self % map % map(p % preCollision)
      state = p
      cIdx = self % map % map(state)

      ! No boundaries crossed => no current
      if (sIdx == cIdx) return

      !! TODO clean below into three function calls
      coordAtEnd = floor((state % r) * (1.0/ (self % spacing)))
      coordAtStart = floor((p % preCollision % r) * (1.0/ (self % spacing)))
      diff = (state % r) - (p % preCollision % r)
      currentContribution = diff / norm2(diff)

      ! X
      stepIdx = sIdx
      step = diff / abs(diff(1))
      if (coordAtEnd(1) - coordAtStart(1) > 0) then
        do i = 1, abs(coordAtEnd(1) - coordAtStart(1)) - 1, 1
          state % r = (p % preCollision % r) + i * step * (self % spacing)

          stepPreviousIdx = stepIdx
          stepIdx = self % map % map(state)
          if (stepIdx == 0) then
            print *, sIdx, cIdx
            print *, state % r
            state = p
            print *, "goal: ", state % r
            print *, "start: ", p % preCollision % r
            print *, "step: ", i * step * (self % spacing)
            stop 37
          end if
          if (stepPreviousIdx == stepIdx) then
            print *, stepIdx, stepPreviousIdx
            stop 99
          end if

          addr = self % getMemAddress() + (0 * self % N) + stepIdx - 1
          call mem % score(currentContribution(1), addr)
        end do
        ! Score the remaining current
        addr = self % getMemAddress() + (0 * self % N) + cIdx - 1
        call mem % score(currentContribution(1), addr)
      end if


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
