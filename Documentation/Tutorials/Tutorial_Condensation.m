[ network,viewer ]=CreateNetwork('GDL_2D');

contactAngle=pi*110/180*ones(network.GetNumberOfLinks);
network.AddNewLinkData(contactAngle,'ContactAngle')


options.TemperatureInletLinks = network.GetLinksFrontiere([1 2 3]); % GDL/MPL interface
options.TemperatureOutletLinks = network.GetLinksFrontiere(5);      % Rib
options.TemperatureInlet = 273+80;  % GDL/MPL interface temperature
options.TemperatureOutlet = 273+75; % rib temperature

options.LiquidWaterOutletLinks = network.GetLinksFrontiere([4 6]); %channel

options.AirPressure = 1.5e5 ;

options.VaporInletLinks = network.GetLinksFrontiere([1 2 3]);  % GDL/MPL interface
options.VaporOutletLinks = network.GetLinksFrontiere([4 6]);   % Channel
options.RelativeHumidityInlet = 0.95;
options.RelativeHumidityOutlet = 0.8;

clusterOptions.Coalescence = 'none' ;
clusterOptions.CapillaryPressureLaw = 'LaplaceCylinder';
clusterOptions.SurfaceTension = 60e-3;
options.ClusterOptions = clusterOptions;


[cluster,outputInformation] = ComputeCondensation( network, options ) ;


viewer.View('PoreField',outputInformation.TemperatureField)

viewer.View('PoreField',outputInformation.PartialVaporPressure{1})

viewer.View('PoreList',cluster.GetInvadedPores)
