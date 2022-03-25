program prog1
   use StoreInput
   use PairsTripletsQuadruplets
   use EnergyCalculation

   implicit none
   
   integer              :: test(3,3)
   integer, allocatable :: PairMatrix(:,:)
   type(atom), allocatable :: DataArray(:)
   real                 :: dummy

   test(1, :) = (/1,2,3/)
   test(2, :) = (/4,5,6/)
   test(3, :) = (/7,8,9/)

   print "(3i3)", test
   print *, test(2,:)
   
   call StoreData("c4h10.xyz", DataArray)
   call FindPairs("c4h10.xyz", PairMatrix)
   dummy = StretchEnergy(PairMatrix, DataArray)
   print *, dummy
end program