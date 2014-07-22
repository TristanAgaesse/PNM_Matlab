function capillaryPressureCurve = ComputeCapillaryPressureCurve(network,numBoundaryInlet,wettability)
%COMPUTECAPILLARYPRESSURECURVE Calcule la courbe de pression capillaire
% input : - network
%           - numBoundaryInlet
%           - wettability : 'currentWettability', 'hydrophobic', 'hydrophilic' or 'random'.
% output : capillaryPressureCurve (saturation en fonction de la pression)
    
    
    [cluster,criticalPressures,breakthroughPressure,invasionPressureList]=ComputeInvasionPercolationReservoir(network,numBoundaryInlet,[],wettability);
    
    
    totalSaturation=ComputeSaturation(cluster,network,'totalSaturation');

    




end

