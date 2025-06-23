function config = loadConfig()
    % LOADCONFIG Initialize configuration parameters for diode model fitting
    %   config = LOADCONFIG() returns a structure containing physical constants
    %   and fitting parameters used in the diode model simulation.
    %
    %   The returned structure contains:
    %   - physics: Physical constants and device parameters
    %   - fitting: Voltage thresholds for different fitting regions
    
    % Physical constants and device parameters
    config.physics = struct(...
        'q', 1.602e-19, ...    % Elementary charge (C)
        'kb', 1.38e-23, ...    % Boltzmann constant (J/K)
        'T', 300, ...          % Temperature (K)
        'n', 1.4, ...          % Ideality factor - increase to better fit negative region
        'm', 2.4 ...           % Exponent for recombination current component
    );
    
    % Calculate thermal voltage factor
    config.physics.A = config.physics.q / (config.physics.kb * config.physics.T);
    
    % Voltage region thresholds for fitting
    config.fitting = struct(...
        'neg_voltage_threshold', -0.2, ... % Negative voltage region threshold (V)
        'pos_voltage_threshold', 0.1 ...   % Positive voltage region threshold (V)
    );
end