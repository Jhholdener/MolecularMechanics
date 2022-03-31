program prog1
   use StoreInput
   use EnergyCalculation

   implicit none
   
   integer              ::  numtriplets, numdihedrals
   type(atom), allocatable :: DataArray(:)
   real                 :: dummy5

   
   call StoreData("c4h10.xyz", DataArray, numtriplets, numdihedrals)

   dummy5 = EnergyFunc(DataArray, numtriplets, numdihedrals)
   
end program