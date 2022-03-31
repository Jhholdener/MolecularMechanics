program prog1
   use StoreInput
   use EnergyCalculation
   use MetropolisAlgorithm

   implicit none
   
   integer              ::  numtriplets, numdihedrals
   type(atom), allocatable :: DataArray(:)
   real                 :: dummy5

   call Metropolis("c4h10.xyz")
   !call StoreData("c4h10.xyz", DataArray, numtriplets, numdihedrals)

   !dummy5 = EnergyFunc(DataArray, numtriplets, numdihedrals)
   
end program