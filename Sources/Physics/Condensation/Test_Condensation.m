function tests = Test_Condensation
%TEST_CONDENSATION Test Condensation
    
    tests = functiontests(localfunctions);
    
end


function test_ComputeTemperatureField(testCase)
    
    %Vary network : cube, canal/rib, extracted from image...
    network = CreateNetwork('GDL_2D');
    
    temperatureInlet = 1;
    temperatureOutlet = 1;
    temperatureInletLinks = network.GetLinksFrontiere([1 2 3]); % GDL/MPL interface
    temperatureOutletLinks = network.GetLinksFrontiere(5);      % Rib 
    [temperature,heatTransferCoefficient] = Condensation_ComputeTemperatureField(network,...
        temperatureInlet,temperatureOutlet,...
        temperatureInletLinks,temperatureOutletLinks);
    
    % Test effective conductivity
    refSolution = 1;
    verifyEqual(testCase,heatTransferCoefficient,refSolution)
    
    
    
    % Test temperature at inlet and outlet
    
end


function test_ComputeEquilibriumVaporPressure(testCase)
    
    % Test with uniform temperature field
    nPore = 10;
    temperature = 1*ones(nPore,1);
    
    equilibriumVaporPressure = Condensation_ComputeEquilibriumVaporPressure(temperature);
    
    % Test field with canal/rib configuration
    
end


function test_Nucleation(testCase)


end

