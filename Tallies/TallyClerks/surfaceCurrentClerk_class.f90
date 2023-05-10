
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
    use geometry_inter, only : geometry

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
    !!
    !! Interface:
    !!   tallyClerk Interface
    !!
    !! Sample dictionary input:
    !!
    !!  clerkName {
    !!      type surfaceCurrentClerk;
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
      procedure :: reportCurrentDueToMovement
    end type surfaceCurrentClerk

    !!
    !! Surface current result class
    !!   Stored in column first order
    !!    dim1 -> energy group
    !!    dim2 -> 1 is J_{x+1/2, y, z}, 2 is J_{x, y+1/2, z}, 3 is J_{x, y, z+1/2}
    !!    dim3 -> Space coordinate
    !!
    type,public, extends(tallyResult) :: SJResult
      integer(shortInt)                           :: NEnergy  = 0
      integer(shortInt)                           :: NSpace  = 0
      real(defReal), dimension(:, :, :),allocatable :: JM ! Current matrix
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

      self % NEnergy = self % energyMap % bins(0)

      ! (xsize + 1) * (ysize + 1) * (zsize + 1)
      self % NSpace = 1
      self % NSpace = self % NSpace * (self % spaceMaps(1) % map % bins(0) + 1)
      self % NSpace = self % NSpace * (self % spaceMaps(2) % map % bins(0) + 1)
      self % NSpace = self % NSpace * (self % spaceMaps(3) % map % bins(0) + 1)

      call dict % get(self % spacing, 'spacing')

    end subroutine init

    !!
    !! Returns array of codes that represent diffrent reports
    !!
    !! See tallyClerk_inter for details
    !!
    function validReports(self) result(validCodes)
      class(surfaceCurrentClerk),intent(in) :: self
      integer(shortInt),dimension(:),allocatable  :: validCodes

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

    subroutine reportCurrentDueToMovement(self, mem, baseAddr, weight, dir, oldR, newR)
      class(surfaceCurrentClerk), intent(inout) :: self
      type(scoreMemory), intent(inout)          :: mem
      integer(longInt) :: baseAddr
      real(defReal) :: weight
      integer(shortInt), intent(in) :: dir
      real(defReal),dimension(3), intent(in) :: oldR, newR
      integer(shortInt) :: stride, index, oldDirIndex, newDirIndex
      type(particleState) :: stateNew, stateOld, negativeState, positiveState
      real(defReal),dimension(3) :: diff
      integer(shortInt) :: negBinIndex, posBinIndex, j
      real(defReal) :: currentContribution

      stateOld % isMG = .false.
      stateOld % r = oldR

      stateNew % isMG = .false.
      stateNew % r = newR

      oldDirIndex = self % spaceMaps(dir) % map % map(stateOld)
      newDirIndex = self % spaceMaps(dir) % map % map(stateNew)

      if (oldDirIndex == newDirIndex) then
        ! No boundaries crossed: no current
        return
      end if

      ! Ensure we travel from negative to positive
      if (oldR(dir) < newR(dir)) then
        negativeState = stateOld
        positiveState = stateNew
      else
        negativeState = stateNew
        positiveState = stateOld
      end if

      ! Lower boundary
      stride = 1
      negBinIndex = 0
      do j = 1, size(self % spaceMaps)
        index = self % spaceMaps(j) % map % map(negativeState)
        if (index == 0) then
          if (j == dir) then
          else
            negBinIndex = -1;
            exit;
          end if
        end if

        negBinIndex = negBinIndex + stride * index
        stride = stride * (self % spaceMaps(j) % map % bins(0) + 1)
      end do

      ! Upper boundary
      stride = 1
      posBinIndex = 0
      do j = 1, size(self % spaceMaps)
        index = self % spaceMaps(j) % map % map(positiveState)
        if (index == 0) then
          if (j == dir) then
            index = self % spaceMaps(j) % map % bins(0) + 1
          else
            posBinIndex = -1;
            exit;
          end if
        end if

        if (j == dir) then
          index = index - 1
        end if

        posBinIndex = posBinIndex + stride * index
        stride = stride * (self % spaceMaps(j) % map % bins(0) + 1)
      end do

      if (negBinIndex == -1 .and. posBinIndex == -1) then
        return
      end if
      if (negBinIndex == -1) then
        negBinIndex = posBinIndex
      end if
      if (posBinIndex == -1) then
        posBinIndex = negBinIndex
      end if

      ! Calculate the actual current
      diff = newR - oldR

      ! Negative contribution since we are backtracking
      currentContribution = -(weight * diff(dir)) / norm2(diff)

      if (posBinIndex == negBinIndex) then
        call mem % score(currentContribution, baseAddr + negBinIndex)
      else
        call mem % score(0.5 * currentContribution, baseAddr + negBinIndex)
        call mem % score(0.5 * currentContribution, baseAddr + posBinIndex)
      end if
      ! TODO: relax this assumption
      ! if (posBinIndex /= negBinIndex) then
      !   print *, "potential inaccuracy"
      !   print *, "dir: ", dir, "r ", negativeState % r
      !   print *, "pos: x: ", self % spaceMaps(1) % map % map(positiveState), "y: ", &
      !     self % spaceMaps(2) % map % map(positiveState), "z: ", self % spaceMaps(3) % map % map(positiveState)
      !   print *, "neg: x: ", self % spaceMaps(1) % map % map(negativeState), "y: ", &
      !     self % spaceMaps(2) % map % map(negativeState), "z: ", self % spaceMaps(3) % map % map(negativeState)
      ! end if

    end subroutine reportCurrentDueToMovement


    subroutine reportCurrentInDirection(self, mem, geom, endP, dir)
      class(surfaceCurrentClerk), intent(inout) :: self
      type(scoreMemory), intent(inout)          :: mem
      class(geometry), intent(in) :: geom
      type(particle), intent(in) :: endP
      integer(shortInt), intent(in) :: dir
      type(particleState) :: endPState, startPState, currentState
      integer(shortInt) :: energyGroup, eventOut
      integer(longInt) :: baseAddr
      type(particle) :: currentParticle
      real(defReal),dimension(3) :: oldR
      real(defReal) :: moveDist, distFromEnd

      endPState = endP
      startPState = endP % preTransition
      energyGroup = self % energyMap % map(startPState)
      if (energyGroup == 0) return  ! We're not interested in this energy group

      ! Index (x, y, z, dir, energy) (energy has the largest stride)
      baseAddr = self % getMemAddress() + (energyGroup - 1) * (3 * self % NSpace) + (dir - 1) * (self % NSpace)

      ! Traverse backward until we reach the starting position
      currentParticle = endP
      call currentParticle % build(endPState % r, -endPState % dir, endPState % E, endPState % wgt)

      oldR = endPState % r
      do
        ! Move partice in the geometry
        !   geom % moveGlobal (also points in correct direction)
        currentState = currentParticle
        distFromEnd = norm2(currentState % r - startPState % r)
        moveDist = min(self % spacing(1), self % spacing(2), self % spacing(3), distFromEnd)
        call geom % moveGlobal(currentParticle % coords, moveDist, eventOut)

        currentState = currentParticle
        call self % reportCurrentDueToMovement(mem, baseAddr, currentState % wgt, dir, oldR, currentState % r)
        oldR = currentState % r

        ! print *, self % spacing, distFromEnd, moveDist
        if (norm2(currentState % r - startPState % r) < 0.01_defReal) then
          ! Have we returned to the start?
          ! If so, transport must have finished
          exit
        end if
      end do

    end subroutine reportCurrentInDirection

    !!
    !! Process incoming transport report
    !!
    !! See tallyClerk_inter for details
    !!
    subroutine reportTrans(self, p, xsData, mem, geom)
      class(surfaceCurrentClerk), intent(inout) :: self
      class(particle), intent(in)                     :: p
      class(nuclearDatabase),intent(inout)            :: xsData
      type(scoreMemory), intent(inout)                :: mem
      class(geometry), intent(in) :: geom

      ! XXX: periodic boundaries: particle behind where its pointing
      call self % reportCurrentInDirection(mem, geom, p, 1)
      call self % reportCurrentInDirection(mem, geom, p, 2)
      call self % reportCurrentInDirection(mem, geom, p, 3)
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
            if(any(shape(res % JM) /= [self % NEnergy, 3, self % NSpace])) then
              deallocate(res % JM)
              allocate(res % JM(self % NEnergy, 3, self % NSpace))
            end if
          else
            allocate(res % JM(self % NEnergy, 3, self % NSpace))
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
                  res % JM(i, j, k) = val
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
      integer(shortInt)                            :: i, xsize, ysize, zsize
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

      xsize = self % spaceMaps(1) % map % bins(0) + 1
      ysize = self % spaceMaps(2) % map % bins(0) + 1
      zsize = self % spaceMaps(3) % map % bins(0) + 1
      call outFile % startArray(name, [xsize, ysize, zsize, 3, self % NEnergy])

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
