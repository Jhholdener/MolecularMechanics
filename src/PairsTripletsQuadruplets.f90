module PairsTripletsQuadruplets
    use StoreInput
    implicit none
    real, parameter :: rCH = 1.090, rCC = 1.526, threshold = 0.10

contains
    subroutine FindPairs(DataArray, PairMatrix)
        integer, allocatable, intent(out)   :: PairMatrix(:,:)
        type(atom), allocatable, intent(in) :: DataArray(:)
        integer                             :: i, j
        real*8                              :: bondlength
        
        ! make a matrix having all atom pairs as 1 and the rest as 0!
        allocate(PairMatrix(size(DataArray),size(DataArray)))
        do i = 1, size(DataArray)
            do j = 1, size(DataArray)
                bondlength = sqrt((DataArray(i)%x-DataArray(j)%x)**2+(DataArray(i)%y-DataArray(j)%y)**2&
                +(DataArray(i)%z-DataArray(j)%z)**2)
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
        !print "(14i3)", PairMatrix
    end subroutine FindPairs

    subroutine FindTriplets(PairMatrix, numtriplets, TripletMatrix)
        integer, allocatable, intent(in)        :: PairMatrix(:,:)
        integer, intent(in)                     :: numtriplets
        integer, allocatable, intent(out)       :: TripletMatrix(:,:) 
        integer                                 :: i, j, k, n, counta=0, countr=0, looplength, newcomb(3)
        integer, allocatable                    :: seencomb(:,:)
        !integer, allocatable                    :: temp(:,:)
        logical                                 :: check

        allocate(TripletMatrix(3,numtriplets))
        allocate(seencomb(3, numtriplets*2))
        looplength = size(PairMatrix(1,:))

        do i = 1, looplength
            do j = 1, looplength
                if (PairMatrix(i,j) == 1) then
                    do k = 1, looplength
                        if (PairMatrix(k,j) == 1 .and. k /= i) then
                            newcomb(1) = i
                            newcomb(2) = j
                            newcomb(3) = k
                            
                            counta = counta + 1
                            seencomb(1,counta) = i
                            seencomb(2,counta) = j
                            seencomb(3,counta) = k
                                                        
                        ! I tried to make it better with allocation, but the deallocation did not seem to work...
                            !if (allocated(seencomb) .eqv. .true.) then
                            !    print *, 'size',ubound(seencomb,2)
                            !    allocate(temp(3, ubound(seencomb,2)))
                            !    temp = seencomb
                            !    deallocate(seencomb)
                            !    
                            !    allocate(seencomb(3, size(temp(1,:))+1))
                            !    print *, 'size',ubound(seencomb,2)
                            !    seencomb(:, 1:size(temp)) = temp
                            !    
                            !
                            !    seencomb(1, ubound(seencomb,2)) = i
                            !    seencomb(2, ubound(seencomb,2)) = j
                            !    seencomb(3, ubound(seencomb,2)) = k
                            !    
                            !    print *,'seencomb =',seencomb
                            !    
                            !    deallocate(temp)

                            !elseif (allocated(seencomb) .eqv. .false.) then
                            !    allocate(seencomb(3,1))
                            !    seencomb(1,1) = i
                            !    seencomb(2,1) = j
                            !    seencomb(3,1) = k
                            !    print *, 'seencomb =',seencomb
                            !endif
                            !print *, 'seencomb =',seencomb

                        !checking if reverse combination already exists
                            check = .true.
                            do n = 1, size(seencomb(1,:))
                                if (newcomb(1) == seencomb(3,n) .and. newcomb(2) == seencomb(2,n)&
                                .and. newcomb(3) == seencomb(1,n)) then
                                    check = .false.
                                endif
                            enddo
                            
                            if (check .eqv. .true.) then
                                countr = countr + 1    
                                TripletMatrix(1,countr) = newcomb(1)
                                TripletMatrix(2,countr) = newcomb(2)
                                TripletMatrix(3,countr) = newcomb(3)
                            endif

                        endif
                    enddo

                endif
            enddo
        enddo
        !print "(3i3)", TripletMatrix
    end subroutine FindTriplets

    subroutine FindQuadruplets(PairMatrix, numdihedrals, DihedralMatrix)
        integer, allocatable, intent(in)        :: PairMatrix(:,:)
        integer, intent(in)                     :: numdihedrals
        integer, allocatable, intent(out)       :: DihedralMatrix(:,:) 
        integer                                 :: i, j, k, l, n, counta=0, countr=0, looplength, newcomb(4)
        integer, allocatable                    :: seencomb(:,:)
        logical                                 :: check

        allocate(DihedralMatrix(4,numdihedrals))
        allocate(seencomb(4, numdihedrals*2))
        looplength = size(PairMatrix(1,:))

        do i = 1, looplength
            do j = 1, looplength
                if (PairMatrix(i,j) == 1) then
                    do k = 1, looplength
                        if (PairMatrix(k,j) == 1 .and. k /= i) then
                            do l = 1, looplength
                                if (PairMatrix(k,l) == 1 .and. l /= j) then
                                    newcomb(1) = i
                                    newcomb(2) = j
                                    newcomb(3) = k
                                    newcomb(4) = l
                            
                                    counta = counta + 1
                                    seencomb(1,counta) = i
                                    seencomb(2,counta) = j
                                    seencomb(3,counta) = k
                                    seencomb(4,counta) = l
                                                        
                                !checking if reverse combination already exists
                                    check = .true.
                                    do n = 1, size(seencomb(1,:))
                                        if (newcomb(1) == seencomb(4,n) .and. newcomb(2) == seencomb(3,n)&
                                        .and. newcomb(3) == seencomb(2,n) .and. newcomb(4) == seencomb(1,n)) then
                                            check = .false.
                                        endif
                                    enddo
                            
                                    if (check .eqv. .true.) then
                                        countr = countr + 1    
                                        DihedralMatrix(1,countr) = newcomb(1)
                                        DihedralMatrix(2,countr) = newcomb(2)
                                        DihedralMatrix(3,countr) = newcomb(3)
                                        DihedralMatrix(4,countr) = newcomb(4)
                                    endif
                                endif
                            enddo                         
                        endif
                    enddo
                endif
            enddo
        enddo
        !print "(4i3)", DihedralMatrix
    end subroutine FindQuadruplets
end module PairsTripletsQuadruplets