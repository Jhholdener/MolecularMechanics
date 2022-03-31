module MetropolisAlgorithm
    use StoreInput
    use EnergyCalculation
    use PairsTripletsQuadruplets
    implicit none
    save
    private
    public :: Metropolis

contains

    subroutine Metropolis(filename)
        character(len=*), intent(in)    :: filename
        type(atom), allocatable         :: DataArray(:), DataArrayNew(:), DataArrayRandom(:)
        integer, allocatable            :: PairMatrix(:,:), TripletMatrix(:,:), DihedralMatrix(:,:)
        
        integer                         :: numtriplets, numdihedrals, i, n, nloops
        real                            :: Energy, Energynew, DeltaEnergy, factorR, randomx, randomy, randomz

        call StoreData(filename, DataArray, numtriplets, numdihedrals)
        
        call FindPairs(DataArray, PairMatrix)
        call FindTriplets(PairMatrix, numtriplets, TripletMatrix)
        call FindQuadruplets(PairMatrix, numdihedrals, DihedralMatrix)

        Energy = EnergyFunc(DataArray, PairMatrix, TripletMatrix, DihedralMatrix)
        

        factorR = 0.001
        nloops = 100
        
        do n = 1, nloops
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
            !print*, Energynew
            !DeltaEnergy = Energy - EnergyNew
        enddo
        print *, DataArrayNew(1)%x, DataArray(1)%x
    end subroutine Metropolis
end module MetropolisAlgorithm