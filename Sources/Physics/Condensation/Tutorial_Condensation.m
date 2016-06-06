[ network,viewer ]=CreateNetwork('GDL_2D');

contactAngle=pi*110/180*ones(network.GetNumberOfLinks,1);
network.AddNewLinkData(contactAngle,'ContactAngle')


options.TemperatureInletLinks = network.GetLinksFrontiere([1 2 3]); % GDL/MPL interface
options.TemperatureOutletLinks = network.GetLinksFrontiere(5);      % Rib
options.TemperatureInlet = 90+273;  % GDL/MPL interface temperature
options.TemperatureOutlet = 70+273; % rib temperature

options.LiquidWaterOutletLinks = network.GetLinksFrontiere([4 6]); %channel

options.AirPressure = 1.5e5 ;

options.VaporInletLinks = network.GetLinksFrontiere([1 2 3]);  % GDL/MPL interface
options.VaporOutletLinks = network.GetLinksFrontiere([4 6]);   % Channel
options.RelativeHumidityInlet = 0.95;
options.RelativeHumidityOutlet = 0.80;

clusterOptions.Coalescence = 'none' ;
clusterOptions.CapillaryPressureLaw = 'LaplaceCylinder';
clusterOptions.SurfaceTension = 60e-3;
options.ClusterOptions = clusterOptions;


%[cluster,outputInformation] = ComputeCondensation( network, options ) ;
Condensation_main

figure
viewer.View('PoreField',outputInformation.TemperatureField)
title('Temperature field')

figure
viewer.View('PoreField',outputInformation.PartialVaporPressure{1})
title('Partial vapor pressure')

figure
viewer.View('PoreList',cluster.GetInvadedPores)
title('Invaded pores')
