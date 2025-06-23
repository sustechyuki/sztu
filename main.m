function main()
    % Diode Model Fitting Tool
    %
    % This program fits experimental current-voltage (I-V) data to a modified diode model
    % that includes series resistance, shunt resistance, and non-ohmic behavior.
    % 
    % Features:
    %   - Load and validate experimental I-V data
    %   - Initialize parameters using Lambert W function or historical parameters
    %   - Perform non-linear fitting to optimize model parameters
    %   - Interactive parameter adjustment for fine-tuning results
    %   - Visualization of fitting results and error analysis
    %   - Export results to various formats (MAT, CSV, TXT, PNG)
    %
    % Model parameters:
    %   - J0: Saturation current density [A]
    %   - Rs: Series resistance [Ohm]
    %   - Rsh: Shunt resistance [Ohm]
    %   - k: Non-ohmic coefficient
    %
    % Author: h.tang;y.zheng;k.luo;
    % Date: 2025-03
    %
    % Load configuration and data
    config = loadConfig();
    [data_V, data_JD] = loadData();
    
    % Validate input data
    validateInputData(data_V, data_JD);
    
    % Ask whether to load initial parameters from history file
    use_history = input('Use historical parameter file as initial parameters? (y/n): ', 's');
    
    if strcmpi(use_history, 'y')
        % 使用历史参数文件
        params = loadHistoricalParameters(data_V, data_JD, config);
    else
        % Initialize parameters - pass data for Lambert W function estimation
        params = initializeParameters(data_V, data_JD, config);
    end
    
    % Perform fitting
    [optimized_params, fit_results] = performFitting(data_V, data_JD, params, config);
    
    % Calculate current components
    currents = calculateCurrents(data_V, optimized_params, config);
    
    % Plot results
    plotResults(data_V, data_JD, fit_results, currents);
    
    % Save fitting results and images
    save_results = input('Save fitting results? (y/n): ', 's');
    if strcmpi(save_results, 'y')
        saveResults(data_V, data_JD, optimized_params, fit_results, currents);
    end
    
    % Output results
    displayResults(optimized_params);
    
    % Ask user if satisfied with fitting results, if not enter interactive adjustment mode
    interactive_adjust = input('Enter interactive parameter adjustment mode? (y/n): ', 's');
    if strcmpi(interactive_adjust, 'y')
        [refined_params, refined_fit] = interactiveParameterAdjustment(data_V, data_JD, optimized_params, config);
        
        % Calculate new current components
        refined_currents = calculateCurrents(data_V, refined_params, config);
        
        % Plot new results
        figure;
        plotResults(data_V, data_JD, refined_fit, refined_currents);
        
        % Output adjusted results
        displayResults(refined_params);
        
        % Save adjusted parameters
        saveAdjustedParameters(refined_params);
    end
end

