

options.RechercheNextInvadedLink='localLinearisationOfCapillaryPressure';
clusterOptions.Coalescence = 'numberOfInvadedNeighbours';

options.nIterations=100;
options.MechanismeDegradation='uniforme';
options.RechercheNextInvadedLink='localLinearisationOfCapillaryPressure';
floodingStepInformationUniforme=ComputeHydrophobicityLoss(network,inletLink,outletLink,options,clusterOptions);
infosUniforme=postTraitementDegradation(network,floodingStepInformationUniforme);

save('DegradationUniforme','infosUniforme')
clear('floodingStepInformationUniforme','infosUniforme')

save('network','network')


network.RemoveLinkData('ContactAngle')
network.RemoveLinkData('InitialContactAngle')
options.nIterations=100;
options.MechanismeDegradation='uniformeDansEau';

floodingStepInformationUniformeEau=ComputeHydrophobicityLoss(network,inletLink,outletLink,options,clusterOptions);
infosUniformeEau=postTraitementDegradation(network,floodingStepInformationUniformeEau);

save('DegradationUniformeEau','infosUniformeEau')
clear('floodingStepInformationUniformeEau','infosUniformeEau')




network.RemoveLinkData('ContactAngle')
network.RemoveLinkData('InitialContactAngle')
options.nIterations=400;
options.MechanismeDegradation='sommeVitesses';

floodingStepInformationVitesse=ComputeHydrophobicityLoss(network,inletLink,outletLink,options,clusterOptions);
infosVitesse=postTraitementDegradation(network,floodingStepInformationVitesse);

save('DegradationVitesse','infosVitesse')
clear('floodingStepInformationVitesse','infosVitesse')



options.RechercheNextInvadedLink='linearDecreaseOfCapillaryPressure';

network.RemoveLinkData('ContactAngle')
network.RemoveLinkData('InitialContactAngle')
options.nIterations=400;
options.MechanismeDegradation='sommeVitesses';

floodingStepInformationVitesse=ComputeHydrophobicityLoss(network,inletLink,outletLink,options,clusterOptions);
infosVitesse=postTraitementDegradation(network,floodingStepInformationVitesse);

save('LinearDecreaseDegradationVitesse','infosVitesse')
clear('floodingStepInformationVitesse','infosVitesse')




clusterOptions.Coalescence = 'none';

options.nIterations=100;
options.MechanismeDegradation='uniforme';
options.RechercheNextInvadedLink='localLinearisationOfCapillaryPressure';
floodingStepInformationUniforme=ComputeHydrophobicityLoss(network,inletLink,outletLink,options,clusterOptions);
infosUniforme=postTraitementDegradation(network,floodingStepInformationUniforme);

save('NoCoalescenceDegradationUniforme','infosUniforme')
clear('floodingStepInformationUniforme','infosUniforme')




network.RemoveLinkData('ContactAngle')
network.RemoveLinkData('InitialContactAngle')
options.nIterations=100;
options.MechanismeDegradation='uniformeDansEau';

floodingStepInformationUniformeEau=ComputeHydrophobicityLoss(network,inletLink,outletLink,options,clusterOptions);
infosUniformeEau=postTraitementDegradation(network,floodingStepInformationUniformeEau);

save('NoCoalescenceDegradationUniformeEau','infosUniformeEau')
clear('floodingStepInformationUniformeEau','infosUniformeEau')




network.RemoveLinkData('ContactAngle')
network.RemoveLinkData('InitialContactAngle')
options.nIterations=400;
options.MechanismeDegradation='sommeVitesses';

floodingStepInformationVitesse=ComputeHydrophobicityLoss(network,inletLink,outletLink,options,clusterOptions);
infosVitesse=postTraitementDegradation(network,floodingStepInformationVitesse);

save('NoCoalescenceDegradationVitesse','infosVitesse')
clear('floodingStepInformationVitesse','infosVitesse')


%save