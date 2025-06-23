function JD = diodeModel(V, params, config)
    % Diode model function: Calculate current density for given voltages
    %
    % Input parameters:
    % V - Voltage array
    % params - Physical parameter array [J0, Rs, Rsh, k]
    %   params(1): J0 - Saturation current density
    %   params(2): Rs - Series resistance
    %   params(3): Rsh - Shunt resistance
    %   params(4): k - Non-ohmic coefficient
    % config - Configuration structure
    %
    % Output parameters:
    % JD - Calculated current density array
    
    % Initialize JD array
    JD = zeros(size(V));
    
    % Check if Rs parameter is positive
    if params(2) <= 0
        error('Physical parameter error: Rs must be positive (current value: %.6e)', params(2));
    end
    
    % Set optimization options
    fsolve_options = optimoptions('fsolve', ...
        'Display', 'off', ...
        'MaxIterations', 2000, ...           
        'FunctionTolerance', 1e-14, ...      
        'OptimalityTolerance', 1e-14, ...    
        'StepTolerance', 1e-14, ...          
        'Algorithm', 'trust-region-dogleg'); 
    
    % Set particle swarm optimization options
    pso_options = optimoptions('particleswarm', ...
        'Display', 'off', ...
        'MaxIterations', 200, ...            
        'SwarmSize', 30, ...                 
        'FunctionTolerance', 1e-10, ...      
        'UseParallel', true, ...             
        'MinNeighborsFraction', 0.2, ...     
        'SelfAdjustmentWeight', 1.5);        
    
    % Set local optimization options
    local_options = optimoptions('fmincon', ...
        'Display', 'off', ...
        'Algorithm', 'interior-point', ...
        'MaxIterations', 500, ...            
        'OptimalityTolerance', 1e-12, ...    
        'StepTolerance', 1e-12, ...          
        'FunctionTolerance', 1e-12, ...      
        'ConstraintTolerance', 1e-10);       
    
    % 删除这里重复的定义
    % 以下代码应该被删除，因为它们会覆盖上面更优的设置
    % pso_options = optimoptions('particleswarm', ...
    %     'Display', 'off', ...
    %     'MaxIterations', 100, ...
    %     'SwarmSize', 20, ...
    %     'FunctionTolerance', 1e-8, ...
    %     'UseParallel', true);
    % 
    % local_options = optimoptions('fmincon', ...
    %     'Display', 'off', ...
    %     'Algorithm', 'interior-point', ...
    %     'MaxIterations', 200, ...
    %     'OptimalityTolerance', 1e-10, ...
    %     'StepTolerance', 1e-10);
    
    % Solve for each voltage point
    for i = 1:length(V)
        % Use different initial guess strategies for negative and positive voltage regions
        if V(i) < 0
            % Negative voltage region, non-ohmic and ohmic terms may be more important
            if i > 1 && V(i-1) < 0
                initial_guess = JD(i-1);
            else
                % Initial guess for negative voltage region, mainly considering non-ohmic and ohmic components
                initial_guess = (V(i) / params(3)) + params(4) * (abs(V(i))^config.physics.m) * sign(V(i));
            end
        else
            % Positive voltage region
            if i > 1
                initial_guess = JD(i-1);
            else
                % Use simple diode term estimate for positive voltage region
                initial_guess = params(1) * (exp(config.physics.A * V(i) / config.physics.n) - 1);
            end
        end

        % Define diode equation
        diode_func = @(J) params(1) * (exp(config.physics.A * (V(i) - J * params(2)) / config.physics.n) - 1) + ...
                    (V(i) - J * params(2)) / params(3) + ...
                    params(4) * (abs(V(i) - J * params(2))^config.physics.m) * sign(V(i) - J * params(2)) - J;
        
        % Define objective function (squared residual)
        obj_func = @(J) diode_func(J)^2;
        
        % Set search bounds
        if V(i) < -0.2
            % Strong negative voltage region, wider search range
            lb = initial_guess * 0.1;  % Lower bound
            ub = initial_guess * 10.0; % Upper bound
            if lb >= 0 && initial_guess < 0
                lb = initial_guess * 10.0;
                ub = initial_guess * 0.1;
            end
        elseif V(i) < 0
            % Weak negative voltage region
            lb = initial_guess * 0.2;  % Lower bound
            ub = initial_guess * 5.0;  % Upper bound
            if lb >= 0 && initial_guess < 0
                lb = initial_guess * 5.0;
                ub = initial_guess * 0.2;
            end
        else
            % Positive voltage region, narrower search range
            lb = max(0, initial_guess * 0.5);  % Lower bound, ensure non-negative
            ub = initial_guess * 2.0;          % Upper bound
        end
        
        % Try using fsolve for quick solution (for simple cases)
        try
            [J_fsolve, ~, exitflag, ~] = fsolve(diode_func, initial_guess, fsolve_options);
            residual_fsolve = abs(diode_func(J_fsolve));
            
            % If fsolve converges with small residual, use the result directly
            if exitflag > 0 && residual_fsolve < 1e-8
                JD(i) = J_fsolve;
                continue;  % Skip subsequent global optimization
            end
        catch
            % If fsolve fails, continue with global optimization
        end
        
        % Use particle swarm optimization for global search
        [J_pso, ~] = particleswarm(obj_func, 1, lb, ub, pso_options);
        
        % Use local optimization to further refine the result
        try
            [J_final, residual_final] = fmincon(obj_func, J_pso, [], [], [], [], lb, ub, [], local_options);
            
            % Check final result residual
            if abs(diode_func(J_final)) < 1e-6
                JD(i) = J_final;
            else
                % If residual is still large, try multiple starting points with fsolve
                best_J = J_final;
                best_residual = abs(diode_func(J_final));
                
                % Try 5 different starting points
                for attempt = 1:5
                    if attempt > 1
                        % Random perturbation of starting point
                        if V(i) < 0
                            rand_factor = 0.5 + rand() * 1.0;
                        else
                            rand_factor = 0.9 + rand() * 0.2;
                        end
                        curr_guess = J_pso * rand_factor;
                    else
                        curr_guess = J_pso;
                    end
                    
                    [J_curr, ~, exitflag] = fsolve(diode_func, curr_guess, fsolve_options);
                    residual = abs(diode_func(J_curr));
                    
                    if residual < best_residual && exitflag > 0
                        best_J = J_curr;
                        best_residual = residual;
                    end
                end
                
                JD(i) = best_J;
            end
        catch
            % If local optimization fails, use particle swarm result
            JD(i) = J_pso;
        end
    end
end