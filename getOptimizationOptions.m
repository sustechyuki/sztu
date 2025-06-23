function options = getOptimizationOptions(algorithm, detail_level)
    % GETOPTIMIZATIONOPTIONS Get optimization options for specified algorithm
    %   OPTIONS = GETOPTIMIZATIONOPTIONS(ALGORITHM, DETAIL_LEVEL) returns optimization options
    %   for the specified algorithm.
    %
    %   Inputs:
    %     ALGORITHM - Optimization algorithm name, options: 'levenberg-marquardt', 'trust-region-reflective'
    %     DETAIL_LEVEL - Display detail level, options: 'none', 'iter', 'iter-detailed'
    %
    %   Outputs:
    %     OPTIONS - Optimization options structure
    
    if nargin < 2
        detail_level = 'iter-detailed';
    end
    
    % Basic options settings
    common_options = {
        'FunctionTolerance', 1e-14, ...    % Increased precision from 1e-12
        'OptimalityTolerance', 1e-14, ...  % Increased precision from 1e-12
        'StepTolerance', 1e-10, ...        % Increased precision from 1e-8
        'MaxFunctionEvaluations', 15000, ... % Increased from 12000
        'MaxIterations', 5000, ...         % Increased from 4000
        'Display', detail_level
    };
    
    % Set specific options based on algorithm
    switch lower(algorithm)
        case 'levenberg-marquardt'
            options = optimoptions('lsqnonlin', ...
                common_options{:}, ...
                'Algorithm', 'levenberg-marquardt', ...
                'FiniteDifferenceType', 'central', ...
                'FiniteDifferenceStepSize', 1e-8);  % Increased precision from 1e-6
            
        case 'trust-region-reflective'
            options = optimoptions('lsqnonlin', ...
                common_options{:}, ...
                'Algorithm', 'trust-region-reflective', ...
                'FiniteDifferenceType', 'central', ...
                'FiniteDifferenceStepSize', 1e-8, ... % Increased precision from 1e-6
                'DiffMaxChange', 1e-1, ...
                'DiffMinChange', 1e-10, ...          % Increased precision from 1e-8
                'SubproblemAlgorithm', 'factorization'); % Added for better numerical stability
            
        otherwise
            error('Unsupported algorithm: %s', algorithm);
    end
end