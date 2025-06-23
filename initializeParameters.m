function params = initializeParameters(data_V, data_JD, config)
    % Initialize parameters for diode model fitting
    %
    % This function estimates initial parameters for the diode model using the Lambert W method.
    % It analyzes different regions of the I-V curve to estimate each parameter:
    % - J0 from moderate forward bias region using Lambert W function
    % - Rs from high forward voltage region slope
    % - Rsh from negative voltage region slope
    % - k is set to a default value for non-ohmic behavior
    %
    % Inputs:
    %   data_V - Measured voltage data
    %   data_JD - Measured current density data
    %   config - Configuration structure with physical constants
    %
    % Output:
    %   params - Structure containing initial parameters, bounds and scale factors
    
    % Initial parameters [J0, Rs, Rsh, k]
    % J0: Saturation current, mainly affects forward characteristics
    % Rs: Series resistance, affects linearity in high current region
    % Rsh: Shunt resistance, mainly affects leakage current in low voltage region
    % k: Non-ohmic coefficient, especially important for fitting negative voltage region
    
    % 获取Lambert W辅助函数
    [solveForJ0Func, ~, ~] = lambertWHelpers();
    
    % Use Lambert W function method to estimate initial parameters
    fprintf('Using Lambert W function to estimate initial parameters...\n');
    
    % Estimate initial Rs value - using slope in high forward voltage region
    pos_idx = find(data_V > 0.2);
    if length(pos_idx) >= 2
        % Use high voltage region to estimate Rs
        [~, max_idx] = max(data_V);
        if max_idx > 1
            Rs_est = abs((data_V(max_idx) - data_V(max_idx-1)) / (data_JD(max_idx) - data_JD(max_idx-1)));
        else
            Rs_est = 1e3; % Default value
        end
    else
        Rs_est = 1e3; % Default value
    end
    
    % Estimate Rsh - using slope in negative voltage region
    neg_idx = find(data_V < -0.2);
    if length(neg_idx) >= 2
        % Use linear fitting to estimate Rsh
        p = polyfit(data_V(neg_idx), data_JD(neg_idx), 1);
        Rsh_est = abs(1/p(1));
    else
        Rsh_est = 1e7; % Default value
    end
    
    % Set ideality factor n
    n = config.physics.n;
    
    % Calculate thermal voltage
    V_th = config.physics.kb * config.physics.T / config.physics.q;
    
    % Use Lambert W method to solve for J0
    % Calculate J0 for each point in the measurement data, then take the average
    J0_values = [];
    
    % Select data points in forward region suitable for J0 estimation
    pos_idx = find(data_V > 0.1 & data_V < 0.25); % Select moderate forward bias region
    
    if ~isempty(pos_idx)
        for i = pos_idx'
            V = data_V(i);
            I = data_JD(i);
            
            % Use Lambert W function formula to calculate J0
            % I = (n*V_th/Rs) * LambertW( (Rs*Rsh*J0)/(n*V_th*(Rs+Rsh)) * exp( (Rsh*(Rs*J0+V))/(n*V_th*(Rs+Rsh)) ) ) - (Rsh*J0-V)/(Rs+Rsh)
            
            % Simplified calculation, ignoring I_ph (photogenerated current)
            % Use iterative method to estimate J0
            J0_est = 1e-9; % Initial guess
            
            % Define objective function: given I and V, solve for J0
            function_handle = @(J0) solveForJ0Func(J0, V, I, Rs_est, Rsh_est, n, V_th);
            
            % Use optimization algorithm to solve for J0
            options = optimoptions('lsqnonlin', 'Display', 'off', 'Algorithm', 'levenberg-marquardt');
            J0_opt = lsqnonlin(function_handle, J0_est, 1e-12, 1e-6, options);
            
            J0_values = [J0_values; J0_opt];
        end
        
        % Calculate average J0
        J0_est = median(J0_values);
    else
        J0_est = 1e-9; % Default value
    end
    
    % Estimate non-ohmic coefficient k - based on negative voltage region data
    k_est = 5e-7; % Initial value, will be optimized later
    
    % Output estimated parameters
    fprintf('Parameters estimated by Lambert W method:\n');
    fprintf('J0 = %.6e A\n', J0_est);
    fprintf('Rs = %.6e Ohm\n', Rs_est);
    fprintf('Rsh = %.6e Ohm\n', Rsh_est);
    
    % Set parameters
    params.x0 = [J0_est, Rs_est, Rsh_est, k_est];
    
    % Parameter range settings
    params.ub = [1e-6, 1e4, 1e10, 1e-5];    % Upper bounds
    params.lb = [1e-12, 1e1, 1e5, 1e-10];    % Lower bounds
    
    % Ensure initial values are within range
    params.x0 = min(max(params.x0, params.lb), params.ub);
    
    % Scale factors to make parameters of similar magnitude
    params.scaleFactors = [1e-9, 1e3, 1e7, 1e-8];
end

% Helper function: objective function for solving J0
function residual = solveForJ0(J0, V, I, Rs, Rsh, n, V_th)
    % Calculate current
    I_calc = calculateCurrentWithLambertW(V, J0, Rs, Rsh, n, V_th);
    
    % Calculate residual
    residual = I_calc - I;
end

% Calculate current using Lambert W function
function I = calculateCurrentWithLambertW(V, J0, Rs, Rsh, n, V_th)
    % Use Lambert W function to calculate current
    % Simplified formula: I = (n*V_th/Rs) * LambertW( (Rs*Rsh*J0)/(n*V_th*(Rs+Rsh)) * exp( (Rsh*(Rs*J0+V))/(n*V_th*(Rs+Rsh)) ) ) - (Rsh*J0-V)/(Rs+Rsh)
    
    % Calculate parameter for Lambert W function
    x = (Rs*Rsh*J0)/(n*V_th*(Rs+Rsh)) * exp((Rsh*(Rs*J0+V))/(n*V_th*(Rs+Rsh)));
    
    % Use MATLAB's lambertw function if available, otherwise use approximation
    if exist('lambertw', 'file')
        w = lambertw(x);
    else
        w = approximateLambertW(x);
    end
    
    % Calculate current
    I = (n*V_th/Rs) * w - (Rsh*J0-V)/(Rs+Rsh);
end

% Approximation implementation of Lambert W function
function w = approximateLambertW(x)
    % Approximate calculation for x >= 0
    if x < 0
        w = 0; % Negative values are not defined, return 0 as default
    elseif x == 0
        w = 0;
    elseif x < 1
        % Use series expansion approximation
        w = x * (1 - x + 1.5*x^2 - 2.667*x^3 + 5.208*x^4);
    else
        % Use iterative method
        % Initial guess
        if x < 3
            w = 0.5;
        else
            w = log(x) - log(log(x));
        end
        
        % Newton iteration
        for i = 1:10
            exp_w = exp(w);
            w_next = w - (w*exp_w - x)/(exp_w + w*exp_w);
            
            % Check convergence
            if abs(w_next - w) < 1e-10
                w = w_next;
                break;
            end
            
            w = w_next;
        end
    end
end