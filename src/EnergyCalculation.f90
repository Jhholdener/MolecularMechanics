module EnergyCalculation
    implicit none
    save
    private

contains
    real function EnergyFunc()

        EnergyFunc = 1
    end function EnergyFunc



    real function StretchEnergy()

        StretchEnergy = 1
    end function StretchEnergy


    real function BendEnergy()

        BendEnergy = 1
    end function BendEnergy

    real function TorsEnergy()

        TorsEnergy = 1
    end function TorsEnergy

    real function NonBondEnergy()
        NonBondEnergy = 1
    end function NonBondEnergy

end module EnergyCalculation