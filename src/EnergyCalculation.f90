module EnergyCalculation
    use StoreInput
    use PairsTripletsQuadruplets
    implicit none
    save
    private
    public :: EnergyFunc

    real, parameter :: kCC = 310.0, kCH = 340.0     !force constants
    !real, parameter :: rCH = 1.090, rCC = 1.526     !equilibrium bond length
    real, parameter :: KthetaCCC = 40.0, KthetaCCH = 50.0, KthetaHCH = 35.0     !Ktheta constants per bond combination
    real, parameter :: thetaCCC = 114.0, thetaCCH = 109.5, thetaHCH = 109.5     !theta constants per bond combination
    real, parameter :: VABCD = 1.40, gamma = 0.0, ntors = 2.0   !constants for torsional energy
    real, parameter :: eC = 0.1094, eH = 0.0157, RC = 1.9080, RH = 1.4590   !constants for Van der Waals term
    real, parameter :: qC = -0.145, qH = 0.145   !partial charge constants 
    real, parameter :: pi = 3.14159


contains
    real function EnergyFunc(DataArray, numtriplets, numdihedrals)
        type(atom), allocatable, intent(in) :: DataArray(:)
        integer, intent(in)                 :: numtriplets, numdihedrals
        integer, allocatable                :: PairMatrix(:,:), TripletMatrix(:,:), DihedralMatrix(:,:)
        
        call FindPairs(DataArray, PairMatrix)
        call FindTriplets(PairMatrix, numtriplets, TripletMatrix)
        call FindQuadruplets(PairMatrix, numdihedrals, DihedralMatrix)

        EnergyFunc = StretchEnergy(PairMatrix, DataArray) + BendEnergy(TripletMatrix, DataArray) +&
            TorsEnergy(DihedralMatrix, DataArray) + NonBondEnergy(PairMatrix, DataArray)
        print *, 'Total energy term is:', EnergyFunc
    end function EnergyFunc



    real function StretchEnergy(PairMatrix, DataArray)
        integer, intent(in), allocatable        :: PairMatrix(:,:)
        type(atom), intent(in), allocatable     :: DataArray(:)
        integer                                 :: looplength, i, j, countr=0
        real                                  :: rAB
        real, allocatable                       :: stretchenergies(:)

        looplength = size(PairMatrix(1,:))
        allocate(stretchenergies(looplength-1))

        do i = 1, looplength
            do j = i+1, looplength
                if (PairMatrix(i,j) == 1) then
                    rAB = sqrt((DataArray(i)%x-DataArray(j)%x)**2+(DataArray(i)%y-DataArray(j)%y)**2&
                    +(DataArray(i)%z-DataArray(j)%z)**2)
                    countr = countr + 1

                    if (DataArray(i)%element == DataArray(j)%element) then
                        stretchenergies(countr) = kCC*((rAB - rCC)**2)
                    elseif (DataArray(i)%element /= DataArray(j)%element) then
                        stretchenergies(countr) = kCH*((rAB - rCH)**2)
                    endif 
                else 
                    cycle
                endif
            enddo
        enddo
        StretchEnergy = sum(stretchenergies)
        print *, 'Stretching energy term is:', StretchEnergy
    end function StretchEnergy


    real function BendEnergy(TripletMatrix, DataArray)
        integer, intent(in), allocatable        :: TripletMatrix(:,:)
        type(atom), intent(in), allocatable     :: DataArray(:)
        integer                                 :: looplength, i, A, B, C, countH, countC
        real, allocatable                       :: bendenergies(:)
        real                                  :: theta
        real                                  :: xBA, yBA, zBA, xBC, yBC, zBC, vBA(3), vBC(3)

        looplength = size(TripletMatrix(1,:))
        allocate(bendenergies(looplength))

        do i = 1, looplength
            A = TripletMatrix(1,i)
            B = TripletMatrix(2,i)
            C = TripletMatrix(3,i)
            xBA = DataArray(A)%x-DataArray(B)%x
            yBA = DataArray(A)%y-DataArray(B)%y
            zBA = DataArray(A)%z-DataArray(B)%z          
            vBA = (/xBA,yBA,zBA/)
            xBC = DataArray(C)%x-DataArray(B)%x
            yBC = DataArray(C)%y-DataArray(B)%y
            zBC = DataArray(C)%z-DataArray(B)%z 
            vBC = (/xBC,yBC,zBC/)

            theta = acosd(dot_product(vBA,vBC)/((sqrt(vBA(1)**2+vBA(2)**2+vBA(3)**2))*(sqrt(vBC(1)**2+vBC(2)**2+vBC(3)**2))))

            countH = 0
            countC = 0
            if (DataArray(A)%element == 'C') then
                countC = countC + 1
            else
                countH = countH + 1
            endif
            if (DataArray(B)%element == 'C') then
                countC = countC + 1
            else
                countH = countH + 1
            endif
            if (DataArray(C)%element == 'C') then
                countC = countC + 1
            else
                countH = countH + 1
            endif

            if (countC == 3) then
                bendenergies(i) = kthetaCCC*((theta-thetaCCC)**2)
            elseif (countC == 2 .and. countH == 1) then
                bendenergies(i) = KthetaCCH*((theta-thetaCCH)**2)
            elseif (countC == 1 .and. countH == 2) then
                bendenergies(i) = KthetaHCH*((theta-thetaHCH)**2)
            endif
        enddo
        BendEnergy = sum(bendenergies)
        print *, 'Bending energy term is:', BendEnergy
    end function BendEnergy


    real function TorsEnergy(DihedralMatrix, DataArray)
        integer, intent(in), allocatable        :: DihedralMatrix(:,:)
        type(atom), intent(in), allocatable     :: DataArray(:)
        integer                                 :: looplength, i, A, B, C, D
        real, allocatable                       :: torsenergies(:)
        real                                  :: theta, xBA, yBA, zBA, xBC, yBC, zBC, xCD, yCD, zCD
        real                                  :: vBA(3), vBC(3), vCB(3), vCD(3), orthABC(3), orthDCB(3)

        looplength = size(DihedralMatrix(1,:))
        allocate(torsenergies(looplength))

        do i = 1, looplength
            A = DihedralMatrix(1,i)
            B = DihedralMatrix(2,i)
            C = DihedralMatrix(3,i)
            D = DihedralMatrix(4,i)
            xBA = DataArray(A)%x-DataArray(B)%x
            yBA = DataArray(A)%y-DataArray(B)%y
            zBA = DataArray(A)%z-DataArray(B)%z          
            vBA = (/xBA,yBA,zBA/)
            xBC = DataArray(C)%x-DataArray(B)%x
            yBC = DataArray(C)%y-DataArray(B)%y
            zBC = DataArray(C)%z-DataArray(B)%z 
            vBC = (/xBC,yBC,zBC/)
            vCB = -vBC
            xCD = DataArray(D)%x-DataArray(C)%x
            yCD = DataArray(D)%y-DataArray(C)%y
            zCD = DataArray(D)%z-DataArray(C)%z 
            vCD = (/xCD,yCD,zCD/)

            orthABC(1) = vBA(2)*vBC(3)-vBA(3)*vBC(2)
            orthABC(2) = vBA(3)*vBC(1)-vBA(1)*vBC(3)
            orthABC(3) = vBA(1)*vBC(2)-vBA(2)*vBC(1)

            orthDCB(1) = vCB(2)*vCD(3)-vCB(3)*vCD(2)
            orthDCB(2) = vCB(3)*vCD(1)-vCB(1)*vCD(3)
            orthDCB(3) = vCB(1)*vCD(2)-vCB(2)*vCD(1)

            theta = acosd(dot_product(orthABC,orthDCB)/((sqrt(orthABC(1)**2+orthABC(2)**2+orthABC(3)**2))*&
                (sqrt(orthDCB(1)**2+orthDCB(2)**2+orthDCB(3)**2))))

            torsenergies(i) = VABCD*abs(1+cosd(ntors*theta-gamma))
        enddo

        TorsEnergy = sum(torsenergies)
        print *, 'Torsional energy term is:', TorsEnergy
    end function TorsEnergy

    real function NonBondEnergy(PairMatrix, DataArray)
        integer, allocatable, intent(in)    :: PairMatrix(:,:)
        type(atom), allocatable, intent(in) :: DataArray(:)
        real, allocatable                   :: vdwenergies(:), nbenergies(:)
        integer                             :: looplength, nonbondedpairs, i, j, countr=0
        real                              :: Rij2, Rij, Rij6, Rij12
        real                                :: eij, RijStar, Aij, Bij, qi, qj
        
        looplength = size(PairMatrix(1,:))
        nonbondedpairs = (looplength**2-looplength)/2-(looplength-1)
        allocate(vdwenergies(nonbondedpairs), nbenergies(nonbondedpairs))

        do i = 1, looplength
            do j = i+1, looplength
                if (PairMatrix(i,j) == 0) then
                    countr = countr+1
                    
                    Rij2 = (DataArray(i)%x-DataArray(j)%x)**2+(DataArray(i)%y-DataArray(j)%y)**2&
                    +(DataArray(i)%z-DataArray(j)%z)**2
                    Rij = sqrt(Rij2)
                    Rij6 = Rij2**3
                    Rij12 = Rij6**2
                    
                    if (DataArray(i)%element == 'C' .and. DataArray(j)%element == 'C') then
                        eij = eC
                        RijStar = 2*RC
                        qi = qC
                        qj = qC
                    elseif (DataArray(i)%element == 'H' .and. DataArray(j)%element == 'H') then
                        eij = eH
                        RijStar = 2*RH
                        qi = qH
                        qj = qH
                    else    !one is C and one is H
                        eij = sqrt(eC*eH)
                        RijStar = RC + RH
                        qi = qC     !the order does not really matter here
                        qj = qH
                    endif

                    Aij = eij*RijStar**12
                    Bij = 2*eij*RijStar**6

                    vdwenergies(countr) = abs(Aij/Rij12 - Bij/Rij6)
                    nbenergies(countr)  = (qi*qj)/(4*pi*eij*Rij)
                endif
            enddo
        enddo
        
        NonBondEnergy = sum(vdwenergies) + sum(nbenergies)
        print *, 'Nonbonded energy term is:', NonBondEnergy
    end function NonBondEnergy

end module EnergyCalculation