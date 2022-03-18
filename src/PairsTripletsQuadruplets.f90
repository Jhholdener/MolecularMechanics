module PairsTripletsQuadruplets
    use StoreInput
    implicit none


contains
    subroutine FindPairs(filename, BondMatrix)
        character(len=*), intent(in)        :: filename
        integer, allocatable, intent(out)   :: BondMatrix(:,:)
        type (atom), allocatable            :: DataArray(:)
        INTEGER                             :: i

        call StoreData(filename, DataArray) 
        !make a matrix having all atom pairs as 1 and the rest as 0!
        do i = 1, size(DataArray)
            print *, DataArray(i)
        enddo


    end subroutine FindPairs
end module PairsTripletsQuadruplets