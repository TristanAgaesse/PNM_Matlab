function Condensation_CheckInputs(network, options)

    assert(isa(network,'PoreNetwork'))
    
    assert(isfield(options,'TemperatureInletLinks'));
    assert(isfield(options,'TemperatureOutletLinks'));
    assert(isfield(options,'TemperatureInlet'));
    assert(isfield(options,'TemperatureOutlet'));
    assert(isfield(options,'LiquidWaterOutletLinks'));
    assert(isfield(options,'AirPressure'));
    assert(isfield(options,'VaporInletLinks'));
    assert(isfield(options,'VaporOutletLinks'));
    assert(isfield(options,'RelativeHumidityInlet'));
    assert(isfield(options,'RelativeHumidityOutlet'));
    assert(isfield(options,'ClusterOptions'));

end