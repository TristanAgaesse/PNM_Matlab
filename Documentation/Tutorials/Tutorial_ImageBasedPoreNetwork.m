%Tutorial : using Image based pore networks

%First load the mat file which contains the results from python pore
%network extraction
%load(extractionResultFile)

%Create the pore network object
InputScriptImageBasedPNM
network=CreateImageBasedPoreNetwork(inputContainerMap);

%Computes conductances for diffusion

diffusivity = 2e-5; % 02 in N2 at ambiant conditions
parameters.PoreBulkProp = diffusivity;
parameters.GeometricModel.Pore = 'Cylinder';
diameterForDiffusion=2*network.GetLinkData('CapillaryRadius');     %CapillaryRadius=inscribed circle
network.AddNewLinkData(diameterForDiffusion,'Diameter')
parameters.LinkBulkProp = 0;
parameters.GeometricModel.Link = 'None';
conductancesDiffusion = LocalScaleComputeConductancesDiffusion(network,parameters);


%Perform a diffusion simulation
transportPores = 1:network.GetNumberOfPores ;

boundaryConditions.inletLink = network.GetLinksFrontiere(1);
boundaryConditions.outletLink = network.GetLinksFrontiere(2);
boundaryConditions.inletType = 'Dirichlet' ;
boundaryConditions.outletType = 'Dirichlet' ;
nInletLink = length(boundaryConditions.inletLink);
nOutletLink = length(boundaryConditions.outletLink);
boundaryConditions.inletValue = 1*ones(1,nInletLink);
boundaryConditions.outletValue = 0.1*ones(1,nOutletLink);

[ concentrations, ~, diffusionCoefficient ] = ComputeLinearTransport( ...
            network,transportPores, ...
            conductancesDiffusion, ...
            boundaryConditions  ...
            );

