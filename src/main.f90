program prog1
   use StoreInput
   Use PairsTripletsQuadruplets

   implicit none
   type(atom), allocatable :: test(:)
   integer, allocatable :: test2(:,:)
   call FindPairs("ch4.xyz", test2)
   !print *, test
end program