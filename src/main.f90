program prog1
   use StoreInput
   use MetropolisAlgorithm

   implicit none
   
   type(atom), allocatable :: NewCoordinates(:)

   call Metropolis("c4h10.xyz", NewCoordinates)

   
end program