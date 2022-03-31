module MetropolisAlgorithm
    use StoreInput
    use EnergyCalculation
    implicit none
    save
    private

contains

    subroutine Metropolis(filename)
        character(len=*), intent(in)    :: filename
        type(atom), allocatable         :: DataArray(:)
        integer                         :: numtriplets, numdihedrals
        real                            :: Energy, Energynew

        call StoreData(filename, DataArray, numtriplets, numdihedrals)
        Energy = EnergyFunc(DataArray, numtriplets, numdihedrals)


    end subroutine Metropolis
end module MetropolisAlgorithm