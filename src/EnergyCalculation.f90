module EnergyCalculation
    use StoreInput
    implicit none
    save
    private
    public :: EnergyFunc, StretchEnergy, BendEnergy, TorsEnergy

    real, parameter :: kCC = 310.0, kCH = 340.0     !force constants
    real, parameter :: rCH = 1.090, rCC = 1.526     !equilibrium bond length
    real, parameter :: KthetaCCC = 40.0, KthetaCCH = 50.0, KthetaHCH = 35.0     !Ktheta constants per bond combination
    real, parameter :: thetaCCC = 114.0, thetaCCH = 109.5, thetaHCH = 109.5     !theta constants per bond combination
    real, parameter :: VABCD = 1.40, gamma = 0.0, ntors = 2.0   !constants for torsional energy


contains
    real function EnergyFunc()

        EnergyFunc = 1
    end function EnergyFunc



    real function StretchEnergy(PairMatrix, DataArray)
        integer, intent(in), allocatable        :: PairMatrix(:,:)
        type(atom), intent(in), allocatable     :: DataArray(:)
        integer                                 :: looplength, i, j, countr=0
        real*8                                  :: rAB
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
    end function StretchEnergy


    real function BendEnergy(TripletMatrix, DataArray)
        integer, intent(in), allocatable        :: TripletMatrix(:,:)
        type(atom), intent(in), allocatable     :: DataArray(:)
        integer                                 :: looplength, i, A, B, C, countH, countC
        real, allocatable                       :: bendenergies(:)
        real*8                                  :: theta
        real*8                                  :: xBA, yBA, zBA, xBC, yBC, zBC, vBA(3), vBC(3)

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
    end function BendEnergy


    real function TorsEnergy(DihedralMatrix, DataArray)
        integer, intent(in), allocatable        :: DihedralMatrix(:,:)
        type(atom), intent(in), allocatable     :: DataArray(:)
        integer                                 :: looplength, i, A, B, C, D
        real, allocatable                       :: torsenergies(:)
        real*8                                  :: theta, xBA, yBA, zBA, xBC, yBC, zBC, xCD, yCD, zCD
        real*8                                  :: vBA(3), vBC(3), vCB(3), vCD(3), orthABC(3), orthDCB(3)

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
    end function TorsEnergy

    real function NonBondEnergy()
        NonBondEnergy = 1
    end function NonBondEnergy

end module EnergyCalculation