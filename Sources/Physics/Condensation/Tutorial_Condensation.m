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
minTemp=60;
maxTemp=100;
temperature = outputInformation.TemperatureField-273;
viewer.ViewPoreData(temperature,[minTemp,maxTemp]);
%viewer.View('PoreField',outputInformation.TemperatureField)
title('Temperature field')

figure
minVP=1e4;
maxVP=1e5;
equilibriumVP = outputInformation.EquilibriumVaporPressure;
viewer.ViewPoreData(equilibriumVP,[minVP,maxVP]);
title('Equilibrium Vapor Pressure')

figure
VPbegin = outputInformation.NucleationInfos.PartialVaporPressure{1};
viewer.ViewPoreData(VPbegin,[minVP,maxVP])
title('Partial vapor pressure begining')

figure
minRH=0;
maxRH=1.5;
RHFieldBegin=outputInformation.NucleationInfos.PartialVaporPressure{1}./outputInformation.EquilibriumVaporPressure;
viewer.ViewPoreData(RHFieldBegin,[minRH,maxRH])
title('Relative humidity begining')

figure
RHFieldEnd=outputInformation.NucleationInfos.PartialVaporPressure{end}./outputInformation.EquilibriumVaporPressure;
viewer.ViewPoreData(RHFieldEnd,[minRH,maxRH])
title('Relative humidity nucleation end')

figure
maxRHEvolution=outputInformation.NucleationInfos.MaxRH;
plot(maxRHEvolution);
xlabel('Nucleation Step');
ylabel('max RH in the network');
title('Max RH evolution during nucleation')
% figure
% viewer.View('PoreList',cluster.GetInvadedPores)
% title('Invaded pores')
