# MolecularMechanics
This program will try to find the minimum energy configuration of single bonded alkanes.
Just fill in your molecule.xyz file into the metropolis subroutine call and the total 
enegry of the system will be calculated. Then the positions of the atoms will be changed
slightly and the energy will be recalculated. If this energy is lower than the previous 
energy, the new configuration will be accepted and is used for the next cycle. The cycle
is infinite until 100 position shifts in a row are not accepted (due to higher energy).

How the program works:
    1. The main program calls the metropolis algortithm with the filename
    2. The metropolis algorithm module stores the data in an atom type function and calls 
        3 modules that determines the Pairs, Triplets and Quadruplets(Dihedrals) that are
        present in the system.
    3. Then the energy function is called to calculate the energy of this system. Within
        this energy function the Stretch-, Bending-, Torsional-, and Nonbonded-energy 
        terms are calculated
    4. The metropolis algorithm then generates random arrays with numbers between -1 and 
        1 to shift the positions of the atoms. These random numbers are reduced to
        decimals with a factor (factorR) and added to the previous accepted positions.
    5. The energy function is called again and the metropolis algorithm module checks if 
        the energy is lower than the previous energy. If so, it is accepted as the new
        position and everything repeats until the new position is declined 100 times in a
        row.
