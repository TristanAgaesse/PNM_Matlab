function equilibriumVaporPressure = Condensation_ComputeEquilibriumVaporPressure(temperature)

    %Antoine's equation used in B.Straubhaar's PhD
    Patm = 101325;
    A = 5.40221;
    B = 1838.675;
    C = -31.737;
    equilibriumVaporPressure = Patm*10.^(A-B./(temperature+C+273)) ;
    
    
    %Rankine formula : http://fr.wikipedia.org/wiki/Pression_de_vapeur_saturante
    %equilibriumVaporPressure = 1e5*exp(13.7-5120./temperature) ;
    
end
