
network=CreateGDLNetwork(dimension,nPore,xLength,yLength,zLength)


%Invasion initiale du reseau
clusterOptions.ThroatPressure = 'LaplaceCylinder' ;
clusterOptions.Coalescence = 'numberOfInvadedNeighbours';
clusterOptions.SurfaceTension = 60e-3; %Water/air surface tension at 80°C

inletLink = network.GetL
outletLink = 

[initialCluster,breakthroughPressure,invasionPressureList] = ComputeInvasionPercolation(network,inletLink,outletLink,'currentWettability',clusterOptions);

%Degradation
options.nIterations = 200 ;
options.ClusterGrowth = false  ;
options.MechanismeDegradation='uniforme','uniformeDansEau','sommeVitesses'
options.RechercheNextInvadedLink='exactTimeForLaplaceLaw' ;

clusterPressure=0.9*invasionPressureList(end);

outputInformation=ComputeHydrophobicityLoss(network,initialCluster,clusterPressure,inletLink,outletLink,options);

infoDegradation=postTraitementDegradation(network,outputInformation);
save('DegradationUniformeEau','infosUniformeEau')



