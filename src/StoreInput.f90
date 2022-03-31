module StoreInput
    implicit none
    private
    public :: Atom, StoreData

    type Atom
        character(len=2)    :: element
        real                :: x, y, z
    end type

contains

    subroutine StoreData(filename, DataArray, numtriplets, numdihedrals) 
        character(len=*)                        :: filename
        type (Atom), allocatable, intent(out)   :: DataArray(:)
        integer, intent(out)                    :: numtriplets, numdihedrals
        integer                                 :: numatoms, i, n, countC=0

        open(unit=20, file=filename, action='read')
        read(20, *) numatoms
        allocate(DataArray(numatoms))
        do i = 1, numatoms
            read(20, *) DataArray(i)%element, DataArray(i)%x, DataArray(i)%y, DataArray(i)%z
        enddo
        close(20)

        do n = 1, numatoms   
            if (DataArray(n)%element == 'C') then
                countC = countC + 1
            endif
        enddo     
        numtriplets = countC*6
        numdihedrals = (countC - 1) * 9
    end subroutine StoreData

end module StoreInput