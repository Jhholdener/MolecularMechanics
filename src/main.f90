program prog1
   use StoreInput
   use PairsTripletsQuadruplets
   use EnergyCalculation

   implicit none
   
   integer              ::  numtriplets, numdihedrals
   integer, allocatable :: PairMatrix(:,:), TripletMatrix(:,:), DihedralMatrix(:,:)
   type(atom), allocatable :: DataArray(:)
   real                 :: dummy1, dummy2, dummy3

   
   call StoreData("c4h10.xyz", DataArray, numtriplets, numdihedrals)
   call FindPairs(DataArray, PairMatrix)
   call FindTriplets(PairMatrix,numtriplets,TripletMatrix)
   call FindQuadruplets(PairMatrix, numdihedrals,DihedralMatrix)
   dummy1 = StretchEnergy(PairMatrix, DataArray)
   dummy2 = BendEnergy(TripletMatrix, DataArray)
   dummy3 = TorsEnergy(DihedralMatrix, DataArray)
   print *, 'stretchenergy:',dummy1
   print *, 'bendenergy:',dummy2
   print *, 'torsenergy:',dummy3

end program