%% Tutorial for water condensation in a gas diffusion layer

%% Set-up and run simulation 

% Set-up the network
folderName=fileparts(which('Tutorial_Condensation.m'));
geometryFile = strcat(folderName,'/GDL_2D');
[ network,viewer ]=CreateNetwork(geometryFile);

poreVolume = network.ComputeAllPoreVolume;
network.AddNewPoreData(poreVolume,'Volume');

contactAngle=pi*110/180*ones(network.GetNumberOfLinks,1);
network.AddNewLinkData(contactAngle,'ContactAngle')

% Set-up the options for condensation algorithm
nPore=network.GetNumberOfPores;

options.TemperatureTransportPores = 1:network.GetNumberOfPores ;
options.TemperatureInletLinks = network.GetLinksFrontiere([1 2 3]); % GDL/MPL interface
options.TemperatureOutletLinks = network.GetLinksFrontiere(5);      % Rib
options.TemperatureInlet = 90+273;                                  % GDL/MPL interface temperature
options.TemperatureOutlet = 50+273;                                 % rib temperature
options.TemperaturePoreHeatConductivity=1*ones(nPore,1);

options.LiquidWaterOutletLinks = network.GetLinksFrontiere([4 6]);  % channel

options.AirPressure = 1.5e5 ;

options.VaporTransportPores = 1:network.GetNumberOfPores;
options.VaporInletLinks = network.GetLinksFrontiere([1 2 3]);       % GDL/MPL interface
options.VaporOutletLinks = network.GetLinksFrontiere([4 6]);        % Channel
options.RelativeHumidityInlet = 0.95;
options.RelativeHumidityOutlet = 0.90;
options.VaporPoreDiffusivity = 2e-5*ones(nPore,1); % 02 in N2 at ambiant conditions

clusterOptions.Coalescence = 'none' ;
clusterOptions.CapillaryPressureLaw = 'LaplaceCylinder';
clusterOptions.SurfaceTension = 60e-3;
options.ClusterOptions = clusterOptions;

% Run condensation simulation
[condensationInfos,condensationClusters]= ComputeCondensation( network, options ) ;


%% Display results

% Plots Nucleation results

figure
minTemp=60;
maxTemp=100;
temperature = condensationInfos.TemperatureField-273;
viewer.ViewPoreData(temperature,[minTemp,maxTemp]);
title('Temperature field')

figure
minVP=1e4;
maxVP=1e5;
equilibriumVP = condensationInfos.EquilibriumVaporPressure;
viewer.ViewPoreData(equilibriumVP,[minVP,maxVP]);
title('Equilibrium Vapor Pressure')

figure
VPbegin = condensationInfos.NucleationInfos.PartialVaporPressure{1};
viewer.ViewPoreData(VPbegin,[minVP,maxVP])
title('Partial vapor pressure begining')

figure
minRH=0;
maxRH=1.5;
RHFieldBegin=condensationInfos.NucleationInfos.PartialVaporPressure{1}./condensationInfos.EquilibriumVaporPressure;
viewer.ViewPoreData(RHFieldBegin,[minRH,maxRH])
title('Relative humidity begining')

figure
RHFieldEnd=condensationInfos.NucleationInfos.PartialVaporPressure{end}./condensationInfos.EquilibriumVaporPressure;
viewer.ViewPoreData(RHFieldEnd,[minRH,maxRH])
title('Relative humidity nucleation end')

figure
maxRHEvolution=condensationInfos.NucleationInfos.MaxRH;
plot(maxRHEvolution);
xlabel('Nucleation Step');
ylabel('max RH in the network');
title('Max RH evolution during nucleation')

figure
effDiff=condensationInfos.NucleationInfos.EffectiveDiffusion;
plot(effDiff);
xlabel('Nucleation Step');
ylabel('Effective diffusion conductivity');
title('Effective diffusion conductivity evolution during nucleation')


% Plots DiffusionControledCondensation results

figure
time=condensationInfos.ClusterGrowthInfos.InvasionTime;
plot(time);
xlabel('Nucleation Step');
ylabel('Invasion Time (s)');
title('Time between two invasions during capillarity controled condensation')


figure
nPore = network.GetNumberOfPores;
invadedPores = zeros(1,nPore);
nInvaded = length(condensationInfos.ClusterGrowthInfos.InvadedPore);
invadedPores(condensationInfos.ClusterGrowthInfos.InvadedPore)=2+(1:nInvaded)/nInvaded;
invadedPores(condensationInfos.NucleationInfos.InvadedPore)=1;
viewer.ViewPoreData(invadedPores,[0,3])
title('Invaded Pores - color is function of invasion order')



