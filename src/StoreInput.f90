module StoreInput
    implicit none
    private
    public :: Atom, StoreData

    type Atom
        character(len=2)    :: element
        real*8              :: x_coord, y_coord, z_coord
    end type

contains

    subroutine StoreData(filename, DataArray) 
        character(len=*)                        :: filename
        type(Atom), allocatable, intent(out)    :: DataArray(:)
        integer                                 :: numatoms, i 

        open(unit=20, file=filename, action='read')
        read(20, *) numatoms
        allocate(DataArray(numatoms))
        do i = 1,numatoms
            read(20, *) DataArray(i)%element, DataArray(i)%x_coord, DataArray(i)%y_coord, DataArray(i)%z_coord
        enddo
        close(20)
    end subroutine StoreData

end module StoreInput