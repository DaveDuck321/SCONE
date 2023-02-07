
module transportXSClerk_class

    use numPrecision
    use tallyCodes
    use endfConstants
    use genericProcedures,          only : fatalError
    use dictionary_class,           only : dictionary
    use particle_class,             only : particle, particleState
    use particleDungeon_class,      only : particleDungeon
    use outputFile_class,           only : outputFile

    ! Nuclear Data interface
    use nuclearDatabase_inter,      only : nuclearDatabase
    use neutronXSPackages_class,    only : neutronMacroXSs
    use neutronMaterial_inter,      only : neutronMaterial,neutronMaterial_CptrCast

    ! Tally Filters
    use tallyFilter_inter,          only : tallyFilter
    use tallyFilterFactory_func,    only : new_tallyFilter

    ! Tally Maps
    use tallyMap_inter,             only : tallyMap
    use tallyMapFactory_func,       only : new_tallyMap

    ! Tally Responses
    use tallyResponseSlot_class,    only : tallyResponseSlot

    ! Tally Interfaces
    use tallyResult_class,          only : tallyResult, tallyResultEmpty
    use scoreMemory_class,          only : scoreMemory
    use tallyClerk_inter,           only : tallyClerk, kill_super => kill

    implicit none
    private

    !! Locations of diffrent bins wrt memory Address of the clerk
    integer(shortInt), parameter  :: MEM_SIZE = 7
    integer(longInt), parameter   :: FLX      = 1 ,&  ! Flux
                                     SCATT    = 2 ,&  ! Total scattering macroscopic reaction rate
                                     CAPT     = 3 ,&  ! Capture macroscopic reaction rate
                                     FISS     = 4 ,&  ! Fission macroscopic reaction rate
                                     NUFISS   = 5 ,&  ! NuFission macroscopic
                                     CHI      = 6 ,&  ! Fission neutron spectrum
                                     SCATT_EV = 7     ! Analog: number of scattering events


    !!
    !! Colision estimator of reaction rates
    !! Calculates flux weighted integral from collisions
    !!
    !! Private Members:
    !!   map      -> Space to store tally Map
    !!   response -> Array of responses
    !!   width    -> Number of responses (# of result bins for each map position)
    !!
    !! Interface
    !!   tallyClerk Interface
    !!
    !! SAMPLE DICTIOANRY INPUT:
    !!
    !! mytransportXSClerk {
    !!   type transportXSClerk;
    !!   energyMap { energyMap definition }
    !!   # spaceMap  { <other tallyMap definition> } #
    !! }
    !!
    type, public, extends(tallyClerk) :: transportXSClerk
      private
      ! Filter, Map & Vector of Responses
      class(tallyMap), allocatable      :: spaceMap
      class(tallyMap), allocatable      :: energyMap

      ! Useful data
      integer(shortInt) :: width   = 0
      integer(shortInt) :: energyN = 1
      integer(shortInt) :: matN    = 0

    contains
      ! Procedures used during build
      procedure  :: init
      procedure  :: kill
      procedure  :: validReports
      procedure  :: getSize

      ! File reports and check status -> run-time procedures
      procedure  :: reportInColl
      procedure  :: reportOutColl
      procedure  :: reportCycleEnd

      ! Output procedures
      procedure  :: print
      procedure  :: getResult
      procedure  :: processResult
      procedure  :: display

    end type transportXSClerk

    !!
    !! multi group cross section result class
    !!
    !! Public Members:
    !!   MGdata(:,:) -> the first dimension represents material/universe number
    !!                  the second dimension the following energy vectors:
    !!                  capture, fission, nu, chi, P0, P1, scattering matrix
    !!
    type,public, extends(tallyResult) :: xsMatrix
      real(defReal), dimension(:,:), allocatable  :: MGdata
    end type xsMatrix

  contains

    !!
    !! Initialise clerk from dictionary and name
    !!
    !! See tallyClerk_inter for details
    !!
    subroutine init(self, dict, name)
      class(transportXSClerk), intent(inout)   :: self
      class(dictionary), intent(in)     :: dict
      character(nameLen), intent(in)    :: name
      character(100), parameter :: Here =' init (transportXSClerk_class.f90)'

      ! Load energy map
      if( dict % isPresent('energyMap')) then
        call new_tallyMap(self % energyMap, dict % getDictPtr('energyMap'))
      else
        call fatalError(Here,'An energy grid is necessary for MG cross section generation')
      end if

      ! Get energy bin number
      self % energyN = self % energyMap % bins(0)

      ! Load space/material map and bin number
      if( dict % isPresent('spaceMap')) then
        call new_tallyMap(self % spaceMap, dict % getDictPtr('spaceMap'))
        self % matN = self % spaceMap % bins(0)
      else
        self % matN = 1
      end if

      ! Set width
      self % width = MEM_SIZE + 3 * self % energyN

    end subroutine init

    !!
    !! Return to uninitialised state
    !!
    elemental subroutine kill(self)
      class(transportXSClerk), intent(inout) :: self

      ! Superclass
      call kill_super(self)
      call self % energyMap % kill()

      ! Kill and deallocate map
      if(allocated(self % spaceMap)) then
        call self % spaceMap % kill()
        deallocate(self % spaceMap)
      end if

      self % width = 0

    end subroutine kill

    !!
    !! Returns array of codes that represent diffrent reports
    !!
    !! See tallyClerk_inter for details
    !!
    function validReports(self) result(validCodes)
      class(transportXSClerk),intent(in)                :: self
      integer(shortInt),dimension(:),allocatable :: validCodes

      validCodes = [inColl_CODE, outColl_CODE, cycleEnd_CODE]

    end function validReports

    !!
    !! Return memory size of the clerk
    !!
    !! See tallyClerk_inter for details
    !!
    elemental function getSize(self) result(S)
      class(transportXSClerk), intent(in) :: self
      integer(shortInt)            :: S

      S = self % width
      S = S * self % energyN * self % matN

    end function getSize

    !!
    !! Process incoming collision report
    !!
    !! See tallyClerk_inter for details
    !!
    subroutine reportInColl(self, p, xsData, mem)
      class(transportXSClerk), intent(inout)       :: self
      class(particle), intent(in)           :: p
      class(nuclearDatabase), intent(inout) :: xsData
      type(scoreMemory), intent(inout)      :: mem
      type(particleState)                   :: state
      type(neutronMacroXSs)                 :: xss
      class(neutronMaterial), pointer       :: mat
      real(defReal)                         :: totalXS, nuFissXS, captXS, fissXS, scattXS, flux
      integer(shortInt)                     :: enIdx, matIdx, binIdx
      integer(longInt)                      :: addr
      character(100), parameter :: Here =' reportInColl (collisionClerk_class.f90)'

      ! Get current particle state
      state = p

      ! Find bin indexes
      enIdx = self % energyMap % map(state)
      if(allocated(self % spaceMap)) then
        matIdx = self % spaceMap % map(state)
      else
        matIdx = 1
      end if

      ! Return if invalid bin index
      if (enIdx == 0 .or. matIdx == 0) return
      binIdx = self % energyN * (matIdx - 1) + enIdx

      ! Calculate bin address
      addr = self % getMemAddress() + self % width * (binIdx - 1) - 1

      mat => neutronMaterial_CptrCast(xsData % getMaterial( p % matIdx()))
      if(.not.associated(mat)) call fatalError(Here,'Unrecognised type of material was retrived from nuclearDatabase')
      call mat % getMacroXSs(xss, p)

      totalXS  = xss % total
      flux = p % w / totalXS

      ! Calculate reaction rates
      nuFissXS = xss % nuFission * flux
      captXS   = xss % capture * flux
      ! TODO: calculate fission h->g
      fissXS   = xss % fission * flux
      scattXS  = (xss % elasticScatter + xss % inelasticScatter) * flux

      ! Add scores to counters
      call mem % score(flux,     addr + FLX)
      call mem % score(nuFissXS, addr + NUFISS)
      call mem % score(captXS,   addr + CAPT)
      call mem % score(fissXS,   addr + FISS)
      call mem % score(scattXS,  addr + SCATT)

    end subroutine reportInColl

    !!
    !! Process outgoing collision report
    !!
    !! See tallyClerk_inter for details
    !!
    subroutine reportOutColl(self, p, MT, muL, xsData, mem)
      class(transportXSClerk), intent(inout) :: self
      class(particle), intent(in)             :: p
      integer(shortInt), intent(in)           :: MT
      real(defReal), intent(in)               :: muL
      class(nuclearDatabase),intent(inout)    :: xsData
      type(scoreMemory), intent(inout)        :: mem
      type(particleState)                     :: preColl, postColl
      real(defReal)                           :: score, prod, mu
      integer(shortInt)                       :: enIdx, matIdx, binIdx, binEnOut
      integer(longInt)                        :: addr

      ! Get current particle state
      preColl = p % preCollision
      postColl = p

      ! Considers weight for scattering multiplicity
      select case(MT)
        case(N_2N, N_2Na, N_2Nd, N_2Nf, N_2Np, N_2N2a, N_2Nl(1):N_2Nl(16))
          score = 2.0_defReal
        case(N_3N, N_3Na, N_3Nf, N_3Np)
          score = 3.0_defReal
        case(N_4N)
          score = 4.0_defReal
        case default
          score = ONE
      end select

      ! Score in case of scattering events
      select case(MT)
      case(N_N_ELASTIC, N_N_INELASTIC, N_Na, N_Np, N_Nd, N_Nt, N_2N, N_2Na, N_2Nd, N_2Nf, &
          N_2Np, N_2N2a, N_2Nl(1):N_2Nl(16), N_3N, N_3Na, N_3Nf, N_3Np, N_4N, N_Nl(1):N_Nl(40), N_Ncont)

          ! Find bin indexes
          enIdx = self % energyMap % map(preColl)
          if(allocated(self % spaceMap)) then
            matIdx = self % spaceMap % map(preColl)
          else
            matIdx = 1
          end if

          ! Return if invalid bin index
          if (enIdx == 0 .or. matIdx == 0) return
          binIdx = self % energyN * (matIdx - 1) + enIdx

          addr = self % getMemAddress() + self % width * (binIdx - 1) - 1

          ! Score one scattering event from group g
          call mem % score(preColl % wgt, addr + SCATT_EV)

          ! Get bin of outgoing energy
          binEnOut = self % energyMap % map(postColl)
          ! Return if invalid bin index
          if (binEnOut == 0) return

          ! Score scattering event from group g to g'
          call mem % score(preColl % wgt, addr + SCATT_EV + binEnOut)

          ! Score outgoing scattering angle for P1 matrix
          mu = muL * preColl % wgt
          call mem % score(mu, addr + SCATT_EV + self % energyN + binEnOut)

          ! Score multiplicity matrix
          prod = score * preColl % wgt
          call mem % score(prod, addr + SCATT_EV + 2*self % energyN + binEnOut)

        case default
          ! Do nothing
      end select

    end subroutine reportOutColl

    !!
    !! Process end of the cycle
    !!
    !! See tallyClerk_inter for details
    !!
    subroutine reportCycleEnd(self, end, mem)
      class(transportXSClerk), intent(inout)     :: self
      class(particleDungeon), intent(in)  :: end
      type(scoreMemory), intent(inout)    :: mem
      integer(longInt)                    :: addr, N, i, binIdx, enIdx, matIdx

      ! Score fission spectrum
      N = end % popSize()
      do i = 1,N

        ! Find bin indexes
        enIdx = self % energyMap % map(end % prisoners(i))
        if(allocated(self % spaceMap)) then
          matIdx = self % spaceMap % map(end % prisoners(i))
        else
          matIdx = 1
        end if

        ! Return if invalid bin index
        if (enIdx == 0 .or. matIdx == 0) cycle
        binIdx = self % energyN * (matIdx - 1) + enIdx

        addr = self % getMemAddress() + self % width * (binIdx - 1) - 1
        call mem % score(ONE,  addr + CHI)
      end do

    end subroutine reportCycleEnd

    !!
    !! Calculates the multi-group constants
    !!
    pure subroutine processResult(self, mem, sigmaF_res, sigmaC_res, transpFL_res, transpOS_res, &
                                  nuBar_res, chiTot_res, P0_res, P1_res, prod_res, fluxG, &
                                  scatteringProductionXS, fissionProductionXS)
      class(transportXSClerk), intent(in)     :: self
      type(scoreMemory), intent(in)    :: mem
      real(defReal), dimension(:), allocatable, intent(out) :: sigmaF_res
      real(defReal), dimension(:), allocatable, intent(out) :: nuBar_res
      real(defReal), dimension(:), allocatable, intent(out) :: sigmaC_res
      real(defReal), dimension(:), allocatable, intent(out) :: transpFL_res
      real(defReal), dimension(:), allocatable, intent(out) :: transpOS_res
      real(defReal), dimension(:), allocatable, intent(out) :: chiTot_res
      real(defReal), dimension(:), allocatable, intent(out) :: P0_res
      real(defReal), dimension(:), allocatable, intent(out) :: P1_res
      real(defReal), dimension(:), allocatable, intent(out) :: prod_res
      real(defReal), dimension(:), allocatable, intent(out) :: fluxG
      real(defReal), dimension(:), allocatable, intent(out) :: scatteringProductionXS
      real(defReal), dimension(:), allocatable, intent(out) :: fissionProductionXS
      real(defReal), dimension(:,:), allocatable  :: delta
      real(defReal), dimension(:), allocatable    :: total
      integer(longInt) :: addr, N, Nm, i, j, k, m, materialEnergyId
      real(defReal)    :: nuF, cap, fis, flux, chiT, P0, P1, prod, sumChi, &
                          scat, scattEv, nuFstd, capStd, fisStd, fluxStd, chiTstd, &
                          P0std, P1std, prodStd, scatStd, scattEvStd, scatXS, scatXSstd, &
                          fissionXs, chiSum

      N = self % energyN
      Nm = self % matN
      ! Allocate arrays for MG xss
      allocate(sigmaF_res(N*Nm), sigmaC_res(N*Nm), nuBar_res(N*Nm), chiTot_res(N*Nm), &
               P0_res(N*N*Nm), P1_res(N*N*Nm), prod_res(N*N*Nm), transpFL_res(N*Nm), transpOS_res(N*Nm), &
               total(N*Nm) , fluxG(N*Nm), delta(Nm, N), scatteringProductionXS(N*N*Nm), fissionProductionXS(N*N*Nm))

      ! Loop over all materials and energies
      sumChi = 0
      m = 1
      delta = ZERO
      do i = 1, N*Nm
        ! Retrieve results from memory
        addr = self % getMemAddress() + self % width * (i - 1) - 1
        call mem % getResult(fis,  fisStd,  addr + FISS)
        call mem % getResult(cap,  capStd,  addr + CAPT)
        call mem % getResult(scat, scatStd, addr + SCATT)
        call mem % getResult(flux, fluxStd, addr + FLX)
        call mem % getResult(nuF,  nuFstd,  addr + NUFISS)
        call mem % getResult(chiT, chiTstd, addr + CHI)  ! Probably mapping fission in cell -> energy level at given index
        call mem % getResult(scattEv, scattEvStd, addr + SCATT_EV)
        ! Calculate MG constants
        ! Account for the risk of NaNs by division by zero
        if (flux == ZERO) then
          sigmaC_res(i) = ZERO
          scatXS    = ZERO
          scatXSstd = ZERO
        else
          sigmaC_res(i) = cap/flux
          scatXS    = scat/flux
          scatXSstd = scatXS * sqrt((scatStd/scat)**2 + (fluxStd/flux)**2)
        end if
        ! Calculate fission production term
        if (fis == ZERO) then
          sigmaF_res(i) = ZERO
          nuBar_res(i) = ZERO
        else
          sigmaF_res(i) = fis/flux
          nuBar_res(i)  = nuF/fis
        end if
        ! Calculate and normalise fission spectrum
        chiTot_res(i) = chiT
        sumChi = sumChi + chiT
        if (mod(i,N) == 0) then
          if (sumChi /= ZERO) chiTot_res((i+1-N):i) = chiTot_res((i+1-N):i)/sumChi
          sumChi = 0
        end if
        ! Store total xs and flux
        total(i) = sigmaC_res(i) + scatXS + sigmaF_res(i)
        fluxG(i) = flux
        ! Loop over outgoing energies
        do j = 1, N
          ! Get results
          call mem % getResult(P0, P0std, addr + SCATT_EV + j)
          call mem % getResult(P1, P1std, addr + SCATT_EV + N + j)
          call mem % getResult(prod, prodStd, addr + SCATT_EV + 2*N + j)
          ! Calculate scattering matrices
          if (P0 == ZERO) then
            P0_res(N*(i-1)+j) = ZERO
            P1_res(N*(i-1)+j) = ZERO
            prod_res(N*(i-1)+j) = ONE
          else
            P0_res(N*(i-1)+j) = P0/scattEv*scatXS
            P1_res(N*(i-1)+j) = P1/scattEv*scatXS
            prod_res(N*(i-1)+j) = prod/P0
            delta(m,j) = delta(m,j) + P1/scattEv*scat
          end if

        end do

        if (mod(i, N) == 0) m = m + 1

        transpOS_res(i) = total(i) - sum(P1_res((N*(i-1)+1):(N*(i-1)+N)))

      end do

      ! XXX: P0 * prod_res is scattering production cross section
      ! XXX: transpFL and transpOS are the same quantity
      scatteringProductionXS = P0_res * prod_res

      ! Calculate fission production rate
      do i = 1, Nm
        do j = 1, N
          materialEnergyId = N * (i - 1) + j
          fissionXs = sigmaF_res(materialEnergyId)
          do k = 1, N
            fissionProductionXS(N * (materialEnergyId - 1) + k) = fissionXs * chiTot_res( N * (i - 1) + k)
          end do
        end do
      end do

      ! Calculate transport cross section with in-scatter approximation
      do i = 1, Nm
        do j = 1, N
          if (fluxG(N*(i-1)+j) == ZERO) then
            transpFL_res(N*(i-1)+j) = ZERO
          else
            transpFL_res(N*(i-1)+j) = total(N*(i-1)+j) - delta(i,j)/fluxG(N*(i-1)+j)
          end if
        end do
      end do

    end subroutine processResult

    !!
    !! Outputs results
    !!
    !! See tallyClerk_inter for details
    !!
    pure subroutine getResult(self, res, mem)
      class(transportXSClerk), intent(in)                     :: self
      class(tallyResult), allocatable, intent(inout)   :: res
      type(scoreMemory), intent(in)                    :: mem
      integer(shortInt)                                :: N, i, j, k
      real(defReal), dimension(:), allocatable :: sigmaF, sigmaC, nuBar, chiT, P0, P1, prod, transFL, &
                                                  transOS, fissionProductionXS, scatteringProductionXS
      real(defReal), dimension(:), allocatable :: flux
      real(defReal), dimension(:,:), allocatable :: matrix

      call self % processResult(mem, sigmaF, sigmaC, transFL, transOS, nuBar, chiT, P0, P1, &
                                prod, flux, scatteringProductionXS, fissionProductionXS)

      N = self % energyN
      allocate(matrix(self % matN, N))

      ! Loop over spatial regions
      do i = 1, self % matN
        ! Loop over energies
        do j = 1, N
          ! TODO: HACK!! I only record the transport cross section
          matrix(i, j) = transFL(N-(j-1) + N*(i-1))
        end do
      end do

      allocate(res, source = xsMatrix(matrix))

    end subroutine getResult

    !!
    !! Write contents of the clerk to output file
    !!
    !! See tallyClerk_inter for details
    !!
    subroutine print(self, outFile, mem)
      class(transportXSClerk), intent(in)               :: self
      class(outputFile), intent(inout)           :: outFile
      type(scoreMemory), intent(in)              :: mem
      integer(shortInt)                          :: i
      integer(shortInt),dimension(:),allocatable :: resArrayShape
      character(nameLen)                         :: name
      real(defReal), dimension(:), allocatable :: sigmaF, sigmaC, nuBar, chiT, &
                                                    P0, P1, prod, transFL, transOS, &
                                                    scatteringProductionXS, flux, fissionProductionXS

      call self % processResult(mem, sigmaF, sigmaC, transFL, transOS, nuBar, chiT, P0, P1,&
                                prod, flux, scatteringProductionXS, fissionProductionXS)

      ! Begin block
      name = 'TransportXS'
      call outFile % startBlock(name)

      ! Allocate space for resultShape array`
      if (allocated(self % spaceMap)) then
        allocate(resArrayShape(self % spaceMap % dimensions() + 1))
      else
        allocate(resArrayShape(1))
      end if

      ! Print energy map information
      call self % energyMap % print(outFile)
      resArrayShape(1) = self % energyN

      ! If TransportXS clerk has a space map print map information
      if (allocated(self % spaceMap)) then
        call self % spaceMap % print(outFile)
        resArrayShape(2:(self % spaceMap % dimensions() + 1)) = self % spaceMap % binArrayShape()
      end if

      ! Print results
      name = 'transportXS'
      call outFile % startArray(name, resArrayShape)
      do i=1,product(resArrayShape)
        call outFile % addResult(transFL(i), ONE)
      end do
      call outFile % endArray()

      deallocate(resArrayShape)
      allocate(resArrayShape(self % spaceMap % dimensions() + 2))

      ! Scattering production matrix shap
      resArrayShape(1) = self % energyN
      resArrayShape(2) = self % energyN
      if (allocated(self % spaceMap)) then
        resArrayShape(3:(self % spaceMap % dimensions() + 2)) = self % spaceMap % binArrayShape()
      end if

      name = 'scatteringProductionXS'
      call outFile % startArray(name, resArrayShape)
      do i=1,product(resArrayShape)
        call outFile % addResult(scatteringProductionXS(i), ONE)
      end do
      call outFile % endArray()

      name = 'fissionProductionXS'
      call outFile % startArray(name, resArrayShape)
      do i=1,product(resArrayShape)
        call outFile % addResult(fissionProductionXS(i), ONE)
      end do
      call outFile % endArray()
      call outFile % endBlock()

    end subroutine print

    !!
    !! Display convergance progress on the console
    !!
    !! See tallyClerk_inter for details
    !!
    subroutine display(self, mem)
      class(transportXSClerk), intent(in)  :: self
      type(scoreMemory), intent(in) :: mem

      print *, 'transportXSClerk does not support display yet'

    end subroutine display

  end module transportXSClerk_class
