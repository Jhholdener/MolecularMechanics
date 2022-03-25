module EnergyCalculation
    use StoreInput
    implicit none
    save
    private
    public :: EnergyFunc, StretchEnergy

    real, parameter :: kCC = 310.0, kCH = 340.0     !force constants
    real, parameter :: rCH = 1.090, rCC = 1.526     !equilibrium bond length
    real, parameter :: KthetaCCC = 40.0, KthetaCCH = 50.0, KthetaHCH = 35.0     !Ktheta constants per bond combination
    real, parameter :: thetaCCC = 114.0, thetaCCH = 109.5, thetaHCH = 109.5     !theta constants per bond combination


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
                    rAB = sqrt((DataArray(i)%x_coord-DataArray(j)%x_coord)**2+(DataArray(i)%y_coord-DataArray(j)%y_coord)**2&
                    +(DataArray(i)%z_coord-DataArray(j)%z_coord)**2)
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


    real function BendEnergy(PairMatrix, DataArray)
        integer, intent(in), allocatable        :: PairMatrix(:,:)
        type(atom), intent(in), allocatable     :: DataArray(:)
        integer                                 :: looplength, i, j, k, countr=0
        real, allocatable                       :: bendenergies(:)
        real                                    :: theta

        looplength = size(PairMatrix(1,:))
        allocate(bendenergies(looplength-1))

        do i = 1, looplength
            do j = i+1, looplength
                
                
                if (PairMatrix(i,j) == 1) then
                    do k = j+1, looplength

                    enddo
                else
                    cycle
                endif


            enddo
        enddo



        BendEnergy = sum(bendenergies)
    end function BendEnergy

    real function TorsEnergy()

        TorsEnergy = 1
    end function TorsEnergy

    real function NonBondEnergy()
        NonBondEnergy = 1
    end function NonBondEnergy

end module EnergyCalculation