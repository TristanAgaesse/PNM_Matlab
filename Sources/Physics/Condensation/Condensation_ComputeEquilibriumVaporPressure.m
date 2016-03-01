function equilibriumVaporPressure = Condensation_ComputeEquilibriumVaporPressure(temperature)
    %http://fr.wikipedia.org/wiki/Pression_de_vapeur_saturante

    %Rankine formula
    equilibriumVaporPressure = 1e5*exp(13.7-5120./temperature) ;

end
