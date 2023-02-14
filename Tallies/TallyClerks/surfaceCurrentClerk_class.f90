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

    type MapWrapper
      class(tallyMap), allocatable :: map
    end type

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
    !!      spaceMapX { uniform spaceMap definition }
    !!      spaceMapY { uniform spaceMap definition }
    !!      spaceMapZ { uniform spaceMap definition }
    !! }
    !!
    type, public, extends(tallyClerk) :: surfaceCurrentClerk
      private
      !! Map defining the discretisation
      class(MapWrapper), dimension(:), allocatable :: spaceMaps
      class(tallyMap), allocatable             :: energyMap

      ! Number of spacial memory bins (spaceMap + 1 in each dim)
      integer(shortInt)                        :: NSpace = 0
      integer(shortInt)                        :: NEnergy = 0
      real(defReal), dimension(:), allocatable :: spacing
    contains
      ! Procedures used during build
      procedure  :: init
      procedure  :: validReports
      procedure  :: getSize

      ! File reports and check status -> run-time procedures
      procedure  :: reportTrans
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

      allocate(self % spaceMaps(3))

      ! Read maps
      call new_tallyMap(self % energyMap, dict % getDictPtr('energyMap'))
      call new_tallyMap(self % spaceMaps(1) % map, dict % getDictPtr('spaceMapX'))
      call new_tallyMap(self % spaceMaps(2) % map, dict % getDictPtr('spaceMapY'))
      call new_tallyMap(self % spaceMaps(3) % map, dict % getDictPtr('spaceMapZ'))

      call dict % get(self % spacing, 'spacing')

      if (size(self % spacing) /= 3) then
        call fatalError(Here, "SurfaceCurrentClerk requires size('spacing') == 3")
      end if

      self % NEnergy = self % energyMap % bins(0)

      ! (xsize + 1) * (ysize + 1) * (zsize + 1)
      self % NSpace = 1
      self % NSpace = self % NSpace * (self % spaceMaps(1) % map % bins(0) + 1)
      self % NSpace = self % NSpace * (self % spaceMaps(2) % map % bins(0) + 1)
      self % NSpace = self % NSpace * (self % spaceMaps(3) % map % bins(0) + 1)

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
      validCodes = [trans_CODE, cycleEnd_Code]

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

    subroutine reportCurrentInDirection(self, mem, weight, state_in, startIn, endIn, dir)
      class(surfaceCurrentClerk), intent(inout) :: self
      type(scoreMemory), intent(inout)          :: mem
      real(defReal), intent(in)                 :: weight
      type(particleState), intent(in)           :: state_in
      real(defReal),dimension(3), intent(in)    :: startIn, endIn
      integer(shortInt), intent(in)             :: dir
      type(particleState)                       :: state
      integer(shortInt)                         :: coordAtEnd, coordAtStart, i, j, energyGroup, maxStepCount
      integer(longInt)                          :: baseAddr, addr, index, lastIndex, stride
      real(defReal)                             :: surfaceCrossSection, currentContribution
      real(defReal),dimension(3)                :: step, diff, startPos, endPos

      ! Copy the state to allow mutation
      state = state_in
      state % r = startIn

      ! Discretize onto the grid
      coordAtStart = self % spaceMaps(dir) % map % map(state)
      state % r = endIn
      coordAtEnd = self % spaceMaps(dir) % map % map(state)

      ! No boundaries crossed => no current
      if (coordAtStart == coordAtEnd) return

      if (coordAtStart == 0) then
        call fatalError("reportCurrentInDirection", "particle starting from illegal position")
      end if

      ! TODO: can the energy of a particle change during a transport?
      energyGroup = self % energyMap % map(state)
      if (energyGroup == 0) return  ! We're not interested in this energy group

      ! Index (x, y, z, dir, energy) (energy has the largest stride)
      baseAddr = self % getMemAddress() + (energyGroup - 1) * (3 * self % NSpace) + (dir - 1) * (self % NSpace)

      diff = endIn - startIn
      surfaceCrossSection = (self % spacing(1) * self % spacing(2) * self % spacing(3)) / (self % spacing(dir))
      currentContribution = (weight * diff(dir)) / (norm2(diff) * surfaceCrossSection)

      ! If we leave the mesh, we should only run for a single extra iteration
      if (coordAtEnd == 0) then
        if (diff(dir) > 0) then
          coordAtEnd = self % spaceMaps(dir) % map % bins(0) + 1
        else
          coordAtEnd = 0
        end if
      end if

      ! Ensure that all math is positive
      if (diff(dir) > 0) then
        startPos = startIn
        endPos = endIn
      else
        ! Swap the start and end coordinates
        startPos = endIn
        endPos = startIn
      end if

      ! Normalized so step(dir) == spacing(dir)
      step = (diff / diff(dir)) * (self % spacing)

      ! Score the current between all intermediate surfaces
      !    NOTE: we are always going from left to right
      if (step(dir) < 0) then
        call fatalError("reportCurrentInDirection", "negative step")
      end if

      ! Do first step to initialize problem
      stride = 1
      state % r = startPos
      lastIndex = 0
      do j = 1, size(self % spaceMaps)
        index = self % spaceMaps(j) % map % map(state)
        lastIndex = lastIndex + stride * index
        stride = stride * (self % spaceMaps(j) % map % bins(0) + 1)

        ! Indicies outside of mapping
        if (index == 0 .and. j /= dir) then
          lastIndex = -1
          exit
        end if
      end do
      if (stride /= self % NSpace) then
        call fatalError("reportCurrentInDirection", "bad index calculation")
      end if

      ! Step through mesh and record current
      do i = 1, abs(coordAtEnd - coordAtStart)
        state % r = startPos + i * step

        if (lastIndex /= -1) then
          addr = baseAddr + lastIndex
          call mem % score(currentContribution, addr)
        end if

        ! Map in space
        stride = 1
        lastIndex = 0
        do j = 1, size(self % spaceMaps)
          index = self % spaceMaps(j) % map % map(state)
          lastIndex = lastIndex + stride * index
          stride = stride * (self % spaceMaps(j) % map % bins(0) + 1)

          ! Indicies outside of mapping
          if (index == 0 .and. j /= dir) then
            lastIndex = -1
            exit
          end if
        end do
        if (stride /= self % NSpace) then
          call fatalError("reportCurrentInDirection", "bad index calculation")
        end if
      end do
    end subroutine reportCurrentInDirection

    !!
    !! Process incoming transport report
    !!
    !! See tallyClerk_inter for details
    !!
    subroutine reportTrans(self, p, xsData, mem)
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
    end subroutine reportTrans

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
      do i=1, size(self % spaceMaps)
        call self % spaceMaps(i) % map % print(outFile)
      end do

      ! Print surface current matrix
      name = 'JM'
      addr = self % getMemAddress() - 1

      call outFile % startArray(name, [self % NSpace, 3, self % NEnergy])

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

      if(allocated(self % spaceMaps)) deallocate(self % spaceMaps)
      self % NSpace = 0
      self % NEnergy = 0
    end subroutine kill

end module surfaceCurrentClerk_class
