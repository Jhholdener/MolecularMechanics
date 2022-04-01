module MetropolisAlgorithm
    use StoreInput
    use EnergyCalculation
    use PairsTripletsQuadruplets
    implicit none
    save
    private
    public :: Metropolis

    real, parameter :: kb = 1.38064852e-23, T = 298

contains

    subroutine Metropolis(filename, NewCoordinates)
        character(len=*), intent(in)            :: filename
        type(atom), allocatable                 :: DataArray(:), DataArrayNew(:), DataArrayRandom(:)
        type(atom), allocatable, intent(out)    :: NewCoordinates(:)
        integer, allocatable                    :: PairMatrix(:,:), TripletMatrix(:,:), DihedralMatrix(:,:)
        integer                                 :: numtriplets, numdihedrals, i, nloops, countr
        real                                    :: Energy, EnergyInitial, Energynew, DeltaEnergy
        real                                    :: factorR, randomx, randomy, randomz, pA, randomq

        call StoreData(filename, DataArray, numtriplets, numdihedrals)
        
        call FindPairs(DataArray, PairMatrix)
        call FindTriplets(PairMatrix, numtriplets, TripletMatrix)
        call FindQuadruplets(PairMatrix, numdihedrals, DihedralMatrix)

        Energy = EnergyFunc(DataArray, PairMatrix, TripletMatrix, DihedralMatrix)
        EnergyInitial = Energy
        

        factorR = 0.0001
        nloops = 100
        countr = 0
        
        do while (countr < 100)
            DataArrayRandom = DataArray
            DataArrayNew = DataArray
            
            call random_number(DataArrayRandom%x)
            call random_number(DataArrayRandom%y)
            call random_number(DataArrayRandom%z)
        
            do i = 1, size(DataArray)
                call random_number(randomx)
                if (randomx > 0.5) then           
                    DataArrayNew(i)%x = DataArray(i)%x + factorR*DataArrayRandom(i)%x
                else
                    DataArrayNew(i)%x = DataArray(i)%x - factorR*DataArrayRandom(i)%x
                endif

                call random_number(randomy)
                if (randomy > 0.5) then           
                    DataArrayNew(i)%y = DataArray(i)%y + factorR*DataArrayRandom(i)%y
                else
                    DataArrayNew(i)%y = DataArray(i)%y - factorR*DataArrayRandom(i)%y
                endif

                call random_number(randomz)
                if (randomz > 0.5) then           
                    DataArrayNew(i)%z = DataArray(i)%z + factorR*DataArrayRandom(i)%z
                else
                    DataArrayNew(i)%z = DataArray(i)%z - factorR*DataArrayRandom(i)%z
                endif
            enddo

            Energynew = EnergyFunc(DataArrayNew, PairMatrix, TripletMatrix, DihedralMatrix)
            DeltaEnergy = EnergyNew - Energy
            print *, 'New total Energy is', Energynew
            
            pA = exp(-DeltaEnergy/(kb*T))

            if (DeltaEnergy < 0) then
                DataArray = DataArrayNew
                Energy = Energynew
                countr = 0
                print *, 'Change accepted'
            else
                call random_number(randomq)
                if (randomq < pA) then
                    DataArray = DataArrayNew
                    Energy = Energynew
                    countr = 0
                    print *, 'Change accepted'
                else
                    countr = countr +1
                    print *, 'Change declined'
                endif
            endif
        enddo
        print *, 'The initial total energy was:           ', EnergyInitial
        print *, 'The final total energy of the system is:', Energy
        NewCoordinates = DataArray
    end subroutine Metropolis
end module MetropolisAlgorithm