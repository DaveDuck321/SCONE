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
    !!   NSpace   -> Number of spatial Bins
    !!   NEnergy  -> Number of energy Bins
    !!   spacing  -> Spacing of the spatial bins
    !!
    !! Interface:
    !!   tallyClerk Interface
    !!
    !! Sample dictionary input:
    !!
    !!  clerkName {
    !!      type surfaceCurrentClerk;
    !!      spacing (0.2, 0.2, 0.2);
    !!      energyMap { energyMap definition }
    !!      spaceMap { uniform spaceMap definition }
    !!
    type, public, extends(tallyClerk) :: surfaceCurrentClerk
      private
      !! Map defining the discretisation
      class(tallyMap), allocatable             :: spaceMap
      class(tallyMap), allocatable             :: energyMap
      integer(shortInt)                        :: NSpace = 0
      integer(shortInt)                        :: NEnergy = 0
      real(defReal), dimension(:), allocatable :: spacing


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
    !!    dim1 -> energy group
    !!    dim2 -> 1 is J_{x+1/2, y, z}, 2 is J_{x, y+1/2, z}, 3 is J_{x, y, z+1/2}
    !!    dim3 -> 1 is values; 2 is STDs
    !!
    type,public, extends(tallyResult) :: SJResult
      integer(shortInt)                           :: NEnergy  = 0
      integer(shortInt)                           :: NSpace  = 0
      real(defReal), dimension(:, :,:,:),allocatable :: JM ! Current matrix
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
      character(100), parameter :: Here = 'init (surfaceCurrentClerk_class.f90)'

      ! Assign name
      call self % setName(name)

      ! Read maps
      call new_tallyMap(self % energyMap, dict % getDictPtr('energyMap'))
      call new_tallyMap(self % spaceMap, dict % getDictPtr('spaceMap'))

      call dict % get(self % spacing, 'spacing')

      if (size(self % spacing) /= 3) then
        call fatalError(Here, "SurfaceCurrentClerk requires size('spacing') == 3")
      end if

      self % NSpace = self % spaceMap % bins(0)
      self % NEnergy = self % energyMap % bins(0)

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

      S = 3 * self % NEnergy * self % NSpace

    end function getSize


    subroutine reportCurrentInDirection(self, mem, weight, state_in, start, end, dir)
      class(surfaceCurrentClerk), intent(inout) :: self
      type(scoreMemory), intent(inout)          :: mem
      real(defReal), intent(in)                 :: weight
      type(particleState), intent(in)           :: state_in
      real(defReal),dimension(3), intent(in)    :: start, end
      integer(shortInt), intent(in)             :: dir
      type(particleState)                       :: state
      integer(shortInt)                         :: coordAtEnd, coordAtStart, i, offsetDueToSign, energyGroup
      integer(longInt)                          :: baseAddr, addr, stepIdx
      real(defReal)                             :: surfaceCrossSection, currentContribution
      real(defReal),dimension(3)                :: step, diff

      ! Copy the state to allow mutation
      state = state_in
      state % r = start

      ! Discretize onto the grid
      coordAtEnd = floor(end(dir) * (1.0/ (self % spacing(dir))))
      coordAtStart = floor(start(dir) * (1.0/ (self % spacing(dir))))

      ! No boundaries crossed => no current
      if (coordAtStart == coordAtEnd) return

      ! TODO: can the energy of a particle change during a transport?
      energyGroup = self % energyMap % map(state)
      if (energyGroup == 0) return  ! We're not interested in this energy group

      baseAddr = self % getMemAddress() + (3 * self % NSpace) * (energyGroup - 1)

      diff = end - start
      surfaceCrossSection = (self % spacing(1) * self % spacing(2) * self % spacing(3)) / (self % spacing(dir))
      currentContribution = (weight * diff(dir)) / (norm2(diff) * surfaceCrossSection)


      if (diff(dir) > 0) then
        offsetDueToSign = -1
      else
        offsetDueToSign = 0
      end if

      ! Normalized so step(dir) == spacing(dir)
      step = (diff / diff(dir)) * (self % spacing)

      ! Score the current between all intermediate surfaces
      do i = 1, abs(coordAtEnd - coordAtStart) - 1
        state % r = start + (i + offsetDueToSign) * step
        stepIdx = self % spaceMap % map(state)

        addr = baseAddr + ((dir - 1) * self % NSpace) + stepIdx - 1
        call mem % score(currentContribution, addr)
      end do

      ! Score the current into the last surface (without stepping past the end point)
      step = end - state % r

      state % r = end + offsetDueToSign * step
      stepIdx = self % spaceMap % map(state)
      addr = baseAddr + ((dir - 1) * self % NSpace) + stepIdx - 1
      call mem % score(currentContribution, addr)

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
      character(100), parameter :: Here = 'reportInColl (surfaceCurrentClerk_class.f90)'

      state = p
      call self % reportCurrentInDirection(mem, p % w, state, p % preCollision % r, state % r, 1)
      call self % reportCurrentInDirection(mem, p % w, state, p % preCollision % r, state % r, 2)
      call self % reportCurrentInDirection(mem, p % w, state, p % preCollision % r, state % r, 3)
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
      integer(shortInt)                             :: i, j, k
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
            if( any(shape(res % JM) /= [self % NEnergy, 3, self % NSpace, 2])) then
              deallocate(res % JM)
              allocate(res % JM(self % NEnergy, 3, self % NSpace, 2))
            end if
          else
            allocate(res % JM(self % NEnergy, 3, self % NSpace, 2))
          end if

          ! Set size of the JM
          res % NSpace = self % NSpace
          res % NEnergy = self % NEnergy

          ! Load entries
          addr = self % getMemAddress() - 1
          do i = 1, self % NEnergy
            do j = 1, 3
              do k = 1, self % NSpace
                addr = addr + 1
                call mem % getResult(val, STD, addr)
                res % JM(i, j, k, 1) = val
                res % JM(i, j, k, 2) = STD
              end do
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
      call self % spaceMap % print(outFile)

      ! Print surface current matrix
      name = 'JM'
      addr = self % getMemAddress() - 1

      call outFile % startArray(name, [self % NEnergy, 3, self % NSpace])

      do i = 1, self % getSize()
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

      if(allocated(self % spaceMap)) deallocate(self % spaceMap)
      self % NSpace = 0
      self % NEnergy = 0
    end subroutine kill

end module surfaceCurrentClerk_class
