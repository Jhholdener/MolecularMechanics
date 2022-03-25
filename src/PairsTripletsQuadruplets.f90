module PairsTripletsQuadruplets
    use StoreInput
    implicit none
    real, parameter :: rCH = 1.090, rCC = 1.526, threshold = 0.10

contains
    subroutine FindPairs(filename, PairMatrix)
        character(len=*), intent(in)        :: filename
        integer, allocatable, intent(out)   :: PairMatrix(:,:)
        type(atom), allocatable             :: DataArray(:)
        integer                             :: i, j
        real*8                              :: bondlength

        call StoreData(filename, DataArray) 
        ! make a matrix having all atom pairs as 1 and the rest as 0!
        allocate(PairMatrix(size(DataArray),size(DataArray)))
        do i = 1, size(DataArray)
            do j = 1, size(DataArray)
                bondlength = sqrt((DataArray(i)%x_coord-DataArray(j)%x_coord)**2+(DataArray(i)%y_coord-DataArray(j)%y_coord)**2&
                +(DataArray(i)%z_coord-DataArray(j)%z_coord)**2)
                if (DataArray(i)%element == DataArray(j)%element) then
                    if (abs(bondlength-rCC) > threshold) then
                        PairMatrix(i,j) = 0
                    elseif ((abs(bondlength-rCC) <= threshold)) then
                        PairMatrix(i,j) = 1
                    endif
                elseif (DataArray(i)%element /= DataArray(j)%element) then
                    if (abs(bondlength-rCH) > threshold) then
                        PairMatrix(i,j) = 0
                    elseif ((abs(bondlength-rCH) <= threshold)) then
                        PairMatrix(i,j) = 1
                    endif
                endif 
            enddo
        enddo
        print "(14i3)", PairMatrix
    end subroutine FindPairs

    subroutine FindTriplets

    end subroutine FindTriplets
end module PairsTripletsQuadruplets