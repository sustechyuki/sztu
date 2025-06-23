function [optimized_params, fit_results] = performFitting(data_V, data_JD, params, config)
    % Advanced Diode Model Parameter Fitting
    %
    % This function performs multi-stage optimization to fit diode model parameters
    % to experimental I-V data. It uses a staged approach to optimize different
    % parameters in different voltage regions, followed by global optimization.
    %
    % Features:
    %   - Multi-stage optimization process targeting specific parameters in each region
    %   - Negative voltage region fitting for Rsh and non-ohmic coefficient k
    %   - Low positive voltage region fitting for saturation current J0
    %   - High positive voltage region fitting for series resistance Rs
    %   - Global optimization using both Levenberg-Marquardt and Trust-Region-Reflective algorithms
    %   - Automatic parameter validation to ensure physical reasonability
    %   - Detailed error analysis for different voltage regions
    %
    % Inputs:
    %   data_V - Measured voltage data
    %   data_JD - Measured current density data
    %   params - Structure containing initial parameters, bounds and scale factors
    %   config - Configuration structure with physical constants
    %
    % Outputs:
    %   optimized_params - Optimized parameter vector [J0, Rs, Rsh, k]
    %   fit_results - Structure containing fitting results and statistics
    
    try
        % Scale initial parameters
        x0_scaled = params.x0 ./ params.scaleFactors;
        
        % Print initial parameters
        fprintf('\nInitial parameters (estimated using Lambert W function):\n');
        fprintf('J0 = %.6e A\n', params.x0(1));
        fprintf('Rs = %.6e Ohm\n', params.x0(2));
        fprintf('Rsh = %.6e Ohm\n', params.x0(3));
        fprintf('k = %.6e\n', params.x0(4));
        
        % Ensure physical reasonability of parameters
        if params.x0(2) <= 0  % Rs must be positive
            fprintf('Warning: Initial Rs is negative or zero, automatically adjusted to positive value\n');
            params.x0(2) = max(params.lb(2), abs(params.x0(2)));
            x0_scaled = params.x0 ./ params.scaleFactors;
        end
        
        % % Set optimization options - try different algorithms
        % fprintf('\nStarting fitting with levenberg-marquardt algorithm...\n');
        % options_lm = optimoptions('lsqnonlin', ...
        %     'Display', 'iter-detailed', ...
        %     'Algorithm', 'levenberg-marquardt', ...
        %     'FunctionTolerance', 1e-10, ...
        %     'OptimalityTolerance', 1e-10, ...
        %     'StepTolerance', 1e-10, ...
        %     'MaxFunctionEvaluations', 8000, ...
        %     'FiniteDifferenceType', 'central', ... % Use central difference for better accuracy
        %     'FiniteDifferenceStepSize', 1e-6, ... % Set appropriate difference step size
        %     'DiffMaxChange', 1e-1, ...
        %     'DiffMinChange', 1e-8, ...
        %     'MaxIterations', 4000);
        
                % Set optimization options - try different algorithms
                fprintf('\nStarting fitting with levenberg-marquardt algorithm...\n');
                options_lm = optimoptions('lsqnonlin', ...
                    'Display', 'iter-detailed', ...
                    'Algorithm', 'levenberg-marquardt', ...
                    'ScaleProblem', 'jacobian',   ...     
                    'FunctionTolerance', 1e-8,   ...   
                    'OptimalityTolerance', 1e-8,...
                    'StepTolerance', 1e-6,    ...     
                    'InitDamping', 0.1,          ...    
                    'MaxFunctionEvaluations', 4000, ...
                    'FiniteDifferenceType', 'central',...
                    'FiniteDifferenceStepSize', 1e-7, ...
                    'DiffMaxChange', 0.2,        ...    
                    'DiffMinChange', 1e-8,       ...     
                    'MaxIterations', 5000);           
                

        % Stage 1: Fit Rsh and k using initial parameters - mainly optimize negative voltage region
        fprintf('\nStage 1: Optimizing Rsh and non-ohmic coefficient k...\n');
        neg_idx = find(data_V < -0.1);
        if ~isempty(neg_idx)
            % Select negative voltage region data
            neg_V = data_V(neg_idx);
            neg_JD = data_JD(neg_idx);
            
            % Define optimization parameter indices - only optimize Rsh and k
            x0_limited = x0_scaled;
            param_mask = [false, false, true, true]; % Only optimize Rsh and k
            
            % Create negative region error function
            neg_errFun = @(x_opt) errorFunctionPartial(x_opt, x0_limited, param_mask, neg_V, neg_JD, params, config);
            
            % Prepare initial parameters containing only Rsh and k
            x0_neg_opt = x0_scaled(param_mask);
            lb_neg = params.lb(param_mask) ./ params.scaleFactors(param_mask);
            ub_neg = params.ub(param_mask) ./ params.scaleFactors(param_mask);
            
            % 执行优化前打印信息
            fprintf('正在执行负电压区域优化，优化参数: Rsh和k...\n');
            
            % Execute optimization
            [x_neg_opt, ~, ~, ~, ~] = lsqnonlin(neg_errFun, x0_neg_opt, lb_neg, ub_neg, options_lm);
            
            fprintf('负电压区域优化完成\n');
            
            % Update optimization results back to complete parameter set
            x0_scaled(param_mask) = x_neg_opt;
            
            % Calculate negative region fitting results
            x_actual_neg = x0_scaled .* params.scaleFactors;
            fit_JD_neg = diodeModel(neg_V, x_actual_neg, config);
            
            % Calculate negative region fitting error
            neg_rel_errors = abs((fit_JD_neg - neg_JD) ./ (abs(neg_JD) + eps)) * 100;
            fprintf('Negative voltage region fitting: Average relative error = %.2f%%\n', mean(neg_rel_errors));
            
            % Output optimization after Rsh and k
            fprintf('Optimized Rsh = %.6e Ohm\n', x_actual_neg(3));
            fprintf('Optimized k = %.6e\n', x_actual_neg(4));
        end
        
        % Stage 2: Subdivide positive voltage region optimization
        fprintf('\nStage 2: Subdivide positive voltage region optimization...\n');
        
        % Distinguish low positive voltage and high positive voltage regions separately
        % low_pos_idx = find(data_V > 0 & data_V <= 0.15);  % Low positive voltage region
        % high_pos_idx = find(data_V > 0.15);               % High positive voltage region
        low_pos_idx = find(data_V > 0 & data_V <= 0.2);  % Low positive voltage region
        high_pos_idx = find(data_V > 0.2);               % High positive voltage region
        
        
        % First optimize low positive voltage region (mainly J0)
        if ~isempty(low_pos_idx)
            fprintf('Optimize low positive voltage region (0 - 0.15V)...\n');
            low_pos_V = data_V(low_pos_idx);
            low_pos_JD = data_JD(low_pos_idx);
            
            % Define parameter mask - mainly optimize J0
            param_mask = [true, false, false, false]; % Only optimize J0
            
            % Create area error function
            low_pos_errFun = @(x_opt) errorFunctionPartial(x_opt, x0_scaled, param_mask, low_pos_V, low_pos_JD, params, config);
            
            % Prepare parameters
            x0_low_pos_opt = x0_scaled(param_mask);
            lb_low_pos = params.lb(param_mask) ./ params.scaleFactors(param_mask);
            ub_low_pos = params.ub(param_mask) ./ params.scaleFactors(param_mask);
            
            % Execute optimization
            [x_low_pos_opt, ~, ~, ~, ~] = lsqnonlin(low_pos_errFun, x0_low_pos_opt, lb_low_pos, ub_low_pos, options_lm);
            
            % Update parameters
            x0_scaled(param_mask) = x_low_pos_opt;
            
            % Calculate low positive voltage region fitting effects
            x_actual_low_pos = x0_scaled .* params.scaleFactors;
            fit_JD_low_pos = diodeModel(low_pos_V, x_actual_low_pos, config);
            low_pos_rel_errors = abs((fit_JD_low_pos - low_pos_JD) ./ (abs(low_pos_JD) + eps)) * 100;
            fprintf('Low positive voltage region fitting: Average relative error = %.2f%%\n', mean(low_pos_rel_errors));
            fprintf('Optimized J0 = %.6e A\n', x_actual_low_pos(1));
        end
        
        % Then optimize high positive voltage region (mainly Rs)
        if ~isempty(high_pos_idx)
            fprintf('Optimize high positive voltage region (>0.15V)...\n');
            high_pos_V = data_V(high_pos_idx);
            high_pos_JD = data_JD(high_pos_idx);
            
            % Define parameter mask - mainly optimize Rs
            param_mask = [false, true, false, false]; % Only optimize Rs
            
            % Create area error function
            high_pos_errFun = @(x_opt) errorFunctionPartial(x_opt, x0_scaled, param_mask, high_pos_V, high_pos_JD, params, config);
            
            % Prepare parameters
            x0_high_pos_opt = x0_scaled(param_mask);
            lb_high_pos = params.lb(param_mask) ./ params.scaleFactors(param_mask);
            ub_high_pos = params.ub(param_mask) ./ params.scaleFactors(param_mask);
            
            % Execute optimization
            [x_high_pos_opt, ~, ~, ~, ~] = lsqnonlin(high_pos_errFun, x0_high_pos_opt, lb_high_pos, ub_high_pos, options_lm);
            
            % Update parameters
            x0_scaled(param_mask) = x_high_pos_opt;
            
            % Validate Rs as positive
            if x0_scaled(2) * params.scaleFactors(2) <= 0
                fprintf('Warning: Rs is negative or zero, adjusting to positive value\n');
                % Use lower boundary value as alternative
                x0_scaled(2) = params.lb(2) / params.scaleFactors(2);
            end
            
            % Calculate high positive voltage region fitting effects
            x_actual_high_pos = x0_scaled .* params.scaleFactors;
            fit_JD_high_pos = diodeModel(high_pos_V, x_actual_high_pos, config);
            high_pos_rel_errors = abs((fit_JD_high_pos - high_pos_JD) ./ (abs(high_pos_JD) + eps)) * 100;
            fprintf('High positive voltage region fitting: Average relative error = %.2f%%\n', mean(high_pos_rel_errors));
            fprintf('Optimized Rs = %.6e Ohm\n', x_actual_high_pos(2));
        end
        
        % Comprehensive optimization of positive voltage region (J0 and Rs)
        pos_idx = find(data_V > 0);
        if ~isempty(pos_idx)
            fprintf('\nComprehensive optimization of positive voltage region...\n');
            % Select positive voltage region data
            pos_V = data_V(pos_idx);
            pos_JD = data_JD(pos_idx);
            
            % Define optimization parameter indices - optimize J0 and Rs
            param_mask = [true, true, false, false]; % Optimize J0 and Rs
            
            % Create positive area error function
            pos_errFun = @(x_opt) errorFunctionPartial(x_opt, x0_scaled, param_mask, pos_V, pos_JD, params, config);
            
            % Prepare parameters
            x0_pos_opt = x0_scaled(param_mask);
            lb_pos = params.lb(param_mask) ./ params.scaleFactors(param_mask);
            ub_pos = params.ub(param_mask) ./ params.scaleFactors(param_mask);
            
            % 执行优化前打印信息
            fprintf('正在执行正偏压区域综合优化，优化参数: J0和Rs...\n');
            
            % Execute optimization
            [x_pos_opt, ~, ~, ~, ~] = lsqnonlin(pos_errFun, x0_pos_opt, lb_pos, ub_pos, options_lm);
            
            fprintf('正偏压区域综合优化完成\n');
            
            % Update optimization results back to complete parameter set
            x0_scaled(param_mask) = x_pos_opt;
            
            % Validate Rs as positive
            if x0_scaled(2) * params.scaleFactors(2) <= 0
                fprintf('Warning: Rs is negative or zero, adjusting to positive value\n');
                % Use lower boundary value as alternative
                x0_scaled(2) = params.lb(2) / params.scaleFactors(2);
            end
            
            % Calculate positive region fitting results
            x_actual_pos = x0_scaled .* params.scaleFactors;
            fit_JD_pos = diodeModel(pos_V, x_actual_pos, config);
            
            % Calculate positive region fitting error
            pos_rel_errors = abs((fit_JD_pos - pos_JD) ./ (abs(pos_JD) + eps)) * 100;
            fprintf('Positive voltage region comprehensive fitting: Average relative error = %.2f%%\n', mean(pos_rel_errors));
            
            % Output optimization after J0 and Rs
            fprintf('Optimized J0 = %.6e A\n', x_actual_pos(1));
            fprintf('Optimized Rs = %.6e Ohm\n', x_actual_pos(2));
        end
        
        % Complete fitting
        fprintf('\nThird stage: Complete fitting...\n');
        % Create error function handle
        fprintf('\n创建全局误差函数句柄...\n');
        errFun = @(x) errorFunction(x, data_V, data_JD, params, config);
        fprintf('全局误差函数句柄创建完成\n');
        
        % ====================== GA遗传算法优化部分 ======================
        fprintf('\n开始执行遗传算法(GA)全局优化...\n');
        % 设置GA优化选项
        ga_options = optimoptions('ga', ...
            'Display', 'iter', ...
            'PopulationSize', 30, ...
            'MaxGenerations', 50, ...
            'EliteCount', 1, ...
            'CrossoverFraction', 0.8, ...
            'MutationFcn', @mutationadaptfeasible, ...
            'FunctionTolerance', 1e-4, ...
            'MaxStallGenerations', 10, ...
            'UseParallel', true);
        
        % 创建GA优化的目标函数（平方和误差）
        ga_objFun = @(x) sum(errFun(x).^2);
        
        % 在GA之前添加多起点尝试
        num_starting_points = 3;
        best_error = Inf;
        best_params = [];
        
        for i = 1:num_starting_points
            if i == 1
                current_x0 = x0_scaled; % 使用原始起点
            else
                % 随机扰动起点，但保持在合理范围内
                perturbation = 0.3 * (2*rand(size(x0_scaled)) - 1); % ±30%扰动
                current_x0 = x0_scaled .* (1 + perturbation);
                % 确保在边界内
                current_x0 = min(max(current_x0, params.lb ./ params.scaleFactors), params.ub ./ params.scaleFactors);
            end
            
            fprintf('尝试起点 #%d...\n', i);
            
            % 使用当前起点执行GA
            [x_ga_i, ~, ~, ~, ~, ~] = ga(ga_objFun, length(current_x0), [], [], [], [], ...
                params.lb ./ params.scaleFactors, params.ub ./ params.scaleFactors, [], ga_options);
            
            % 计算误差
            temp_params = x_ga_i .* params.scaleFactors;
            temp_JD = diodeModel(data_V, temp_params, config);
            temp_errors = abs((temp_JD - data_JD) ./ (abs(data_JD) + eps)) * 100;
            
            % 比较结果
            if mean(temp_errors) < best_error
                best_error = mean(temp_errors);
                best_params = x_ga_i;
                fprintf('找到更好的起点，平均误差: %.4f%%\n', best_error);
            end
        end
        
        % 使用最佳起点继续GA优化
        x_ga = best_params;
        
        % 检查GA优化结果
        optimized_params_ga = x_ga .* params.scaleFactors;
        fit_results_ga.JD = diodeModel(data_V, optimized_params_ga, config);
        relative_errors_ga = abs((fit_results_ga.JD - data_JD) ./ (abs(data_JD) + eps)) * 100;
        
        fprintf('遗传算法优化完成\n');
        fprintf('GA优化平均相对误差: %.4f%%\n', mean(relative_errors_ga));
        fprintf('参数优化结果:\n');
        fprintf('J0 = %.6e A\n', optimized_params_ga(1));
        fprintf('Rs = %.6e Ohm\n', optimized_params_ga(2));
        fprintf('Rsh = %.6e Ohm\n', optimized_params_ga(3));
        fprintf('k = %.6e\n', optimized_params_ga(4));
        
        % 使用GA结果作为Levenberg-Marquardt的初始值
        fprintf('将使用GA优化结果作为Levenberg-Marquardt优化的初始值\n');
        x0_scaled = x_ga;
        
        % 执行优化前打印信息
        fprintf('正在执行Levenberg-Marquardt全局优化，优化所有参数: J0, Rs, Rsh, k...\n');
        
        % Execute optimization - Levenberg-Marquardt algorithm
        [x_scaled_optimized, resnorm_lm, residual_lm, exitflag_lm, output_lm] = ...
            lsqnonlin(errFun, x0_scaled, [], [], options_lm);
            
        fprintf('Levenberg-Marquardt全局优化完成\n');
        
        % Validate Rs as positive
        if x_scaled_optimized(2) * params.scaleFactors(2) <= 0
            fprintf('Warning: LM algorithm produced negative value or zero of Rs, adjusting to positive value\n');
            % Use lower boundary value as alternative
            x_scaled_optimized(2) = params.lb(2) / params.scaleFactors(2);
        end
        
        % Store L-M algorithm results
        optimized_params_lm = x_scaled_optimized .* params.scaleFactors;
        fit_results_lm.JD = diodeModel(data_V, optimized_params_lm, config);
        fit_results_lm.resnorm = resnorm_lm;
        relative_errors_lm = abs((fit_results_lm.JD - data_JD) ./ (abs(data_JD) + eps)) * 100;
        
        % % Try to follow trust-region-reflective algorithm
        % fprintf('\nFourth stage: using trust-region-reflective algorithm for fitting...\n');
        % options_tr = optimoptions('lsqnonlin', ...
        %     'Display', 'iter-detailed', ...
        %     'Algorithm', 'trust-region-reflective', ...
        %     'FunctionTolerance', 1e-12, ...
        %     'OptimalityTolerance', 1e-12, ...
        %     'StepTolerance', 1e-12, ...
        %     'MaxFunctionEvaluations', 8000, ...
        %     'FiniteDifferenceType', 'central', ... % Use center difference
        %     'FiniteDifferenceStepSize', 1e-6, ... % Set difference step size
        %     'MaxIterations', 4000);
                % Try to follow trust-region-reflective algorithm
                fprintf('\nFourth stage: using trust-region-reflective algorithm for fitting...\n');
                options_tr = optimoptions('lsqnonlin', ...
                    'Display', 'iter-detailed', ...
                    'Algorithm', 'trust-region-reflective', ...
                    'FunctionTolerance', 1e-10,     ...  % 提高函数容差
                    'OptimalityTolerance', 1e-10,   ...  % 提高最优性容差
                    'StepTolerance', 1e-10,         ...   % 减小步长容差
                    'MaxFunctionEvaluations', 4000,  ...% 增加最大函数评估次数
                    'FiniteDifferenceType', 'central',...
                    'FiniteDifferenceStepSize', 1e-6, ...% 减小差分步长提高精度
                    'SubproblemAlgorithm', 'factorization', ...% 使用因式分解求解子问题
                    'MaxIterations', 5000);           % 增加最大迭代次数
                

        % Make sure boundary of Rs is positive
        if params.lb(2) <= 0
            fprintf('Warning: Rs boundary is negative or zero, adjusting to positive value\n');
            params.lb(2) = 10; % Set to minimum 10 Ohm
        end
        
        % Execute optimization - trust-region-reflective algorithm, using L-M result as initial value
        fprintf('\n正在执行Trust-Region-Reflective全局优化，优化所有参数: J0, Rs, Rsh, k...\n');
        fprintf('使用Levenberg-Marquardt优化结果作为初始值...\n');
        [x_scaled_optimized_tr, resnorm_tr, residual_tr, exitflag_tr, output_tr] = ...
            lsqnonlin(errFun, x_scaled_optimized, params.lb ./ params.scaleFactors, params.ub ./ params.scaleFactors, options_tr);
            
        fprintf('Trust-Region-Reflective全局优化完成\n');
        fprintf('优化退出标志: %d, 迭代次数: %d, 函数评估次数: %d\n', exitflag_tr, output_tr.iterations, output_tr.funcCount);
        
        % Validate Rs as positive (for trust-region-reflective algorithm, this should automatically be guaranteed as we set the lower boundary)
        if x_scaled_optimized_tr(2) * params.scaleFactors(2) <= 0
            fprintf('Warning: TR algorithm produced negative value or zero of Rs, adjusting to positive value\n');
            % Use lower boundary value as alternative
            x_scaled_optimized_tr(2) = params.lb(2) / params.scaleFactors(2);
        end
        
        % Compute trust-region algorithm results
        optimized_params_tr = x_scaled_optimized_tr .* params.scaleFactors;
        fit_results_tr.JD = diodeModel(data_V, optimized_params_tr, config);
        fit_results_tr.resnorm = resnorm_tr;
        relative_errors_tr = abs((fit_results_tr.JD - data_JD) ./ (abs(data_JD) + eps)) * 100;
        
        % 在执行BFGS优化前，比较LM和TR结果，选择更好的那个作为BFGS输入
        fprintf('\n比较Levenberg-Marquardt和Trust-Region-Reflective算法结果...\n');
        fprintf('Levenberg-Marquardt平均相对误差: %.4f%%\n', mean(relative_errors_lm));
        fprintf('Trust-Region-Reflective平均相对误差: %.4f%%\n', mean(relative_errors_tr));
        
        if mean(relative_errors_lm) < mean(relative_errors_tr)
            fprintf('Levenberg-Marquardt算法效果更好，将用作BFGS优化的初始值\n');
            best_pre_bfgs_x = x_scaled_optimized;
            best_pre_bfgs_params = optimized_params_lm;
        else
            fprintf('Trust-Region-Reflective算法效果更好，将用作BFGS优化的初始值\n');
            best_pre_bfgs_x = x_scaled_optimized_tr;
            best_pre_bfgs_params = optimized_params_tr; 
        end
        
        % ====================== BFGS局部优化部分 ======================
        fprintf('\n开始执行BFGS局部优化...\n');
        % 创建BFGS目标函数（平方和）
        bfgs_objFun = @(x) sum(errFun(x).^2);
        
        % 设置BFGS优化选项
        options_bfgs = optimoptions('fminunc', ...
            'Display', 'iter', ...
            'Algorithm', 'quasi-newton', ...
            'HessUpdate', 'bfgs', ...
            'OptimalityTolerance', 1e-12, ...    % 提高最优性容差
            'StepTolerance', 1e-10,   ...        % 减小步长容差
            'FunctionTolerance', 1e-12,   ...    % 提高函数容差
            'MaxFunctionEvaluations', 3000,  ... % 增加最大函数评估次数
            'MaxIterations', 3000,       ...     % 增加最大迭代次数
            'FiniteDifferenceType', 'central',... % 使用中心差分
            'FiniteDifferenceStepSize', 1e-7,... % 减小差分步长提高精度
            'CheckGradients', true,     ...      % 检查梯度计算
            'SpecifyObjectiveGradient', false); % 不使用提供的梯度(使用数值梯度)
        
        % 执行BFGS优化，使用LM和TR中更好的那个作为初始值
        fprintf('正在执行BFGS局部优化，优化所有参数...\n');
        [x_scaled_optimized_bfgs, ~, exitflag_bfgs, output_bfgs] = fminunc(bfgs_objFun, best_pre_bfgs_x, options_bfgs);
        
        fprintf('BFGS局部优化完成\n');
        fprintf('优化退出标志: %d, 迭代次数: %d, 函数评估次数: %d\n', exitflag_bfgs, output_bfgs.iterations, output_bfgs.funcCount);
        
        % 检查BFGS优化结果，确保Rs为正
        if x_scaled_optimized_bfgs(2) * params.scaleFactors(2) <= 0
            fprintf('Warning: BFGS algorithm produced negative or zero Rs, adjusting to positive value\n');
            x_scaled_optimized_bfgs(2) = params.lb(2) / params.scaleFactors(2);
        end
        
        % 计算BFGS优化结果
        optimized_params_bfgs = x_scaled_optimized_bfgs .* params.scaleFactors;
        fit_results_bfgs.JD = diodeModel(data_V, optimized_params_bfgs, config);
        fit_results_bfgs.resnorm = sum((fit_results_bfgs.JD - data_JD).^2);
        relative_errors_bfgs = abs((fit_results_bfgs.JD - data_JD) ./ (abs(data_JD) + eps)) * 100;
        
        fprintf('BFGS优化平均相对误差: %.4f%%\n', mean(relative_errors_bfgs));
        
        % Compare two algorithms results, choose the best one
        fprintf('\nCompare two algorithms results:\n');
        fprintf('Levenberg-Marquardt: Average relative error = %.2f%%\n', mean(relative_errors_lm));
        fprintf('Trust-Region-Reflective: Average relative error = %.2f%%\n', mean(relative_errors_tr));
        fprintf('BFGS: Average relative error = %.2f%%\n', mean(relative_errors_bfgs));
        
        % Check the fitting effects of each voltage region
        neg_idx = find(data_V < 0);
        pos_idx = find(data_V > 0);
        
        % LM algorithm performance in different regions
        neg_err_lm = relative_errors_lm(neg_idx);
        pos_err_lm = relative_errors_lm(pos_idx);
        
        % TR algorithm performance in different regions
        neg_err_tr = relative_errors_tr(neg_idx);
        pos_err_tr = relative_errors_tr(pos_idx);
        
        % BFGS algorithm performance in different regions
        neg_err_bfgs = relative_errors_bfgs(neg_idx);
        pos_err_bfgs = relative_errors_bfgs(pos_idx);
        
        fprintf('\nLevenberg-Marquardt: Negative region error = %.2f%%, Positive region error = %.2f%%\n', ...
            mean(neg_err_lm), mean(pos_err_lm));
        fprintf('Trust-Region-Reflective: Negative region error = %.2f%%, Positive region error = %.2f%%\n', ...
            mean(neg_err_tr), mean(pos_err_tr));
        fprintf('BFGS: Negative region error = %.2f%%, Positive region error = %.2f%%\n', ...
            mean(neg_err_bfgs), mean(pos_err_bfgs));
        
        % Adaptively select better algorithm results
        fprintf('\n正在比较各算法的性能...\n');
        fprintf('GA平均相对误差: %.4f%%\n', mean(relative_errors_ga));
        fprintf('Levenberg-Marquardt平均相对误差: %.4f%%\n', mean(relative_errors_lm));
        fprintf('Trust-Region-Reflective平均相对误差: %.4f%%\n', mean(relative_errors_tr));
        fprintf('BFGS平均相对误差: %.4f%%\n', mean(relative_errors_bfgs));
        
        % 直接使用BFGS优化结果
        fprintf('使用BFGS优化结果作为最终拟合参数\n');
        optimized_params = optimized_params_bfgs;
        fit_results.JD = fit_results_bfgs.JD;
        fit_results.resnorm = fit_results_bfgs.resnorm;
        fit_results.residual = (fit_results_bfgs.JD - data_JD);
        fit_results.exitflag = exitflag_bfgs;
        fit_results.output = output_bfgs;
        relative_errors = relative_errors_bfgs;
        
        % Final validation that all parameters are within physically reasonable ranges
        if optimized_params(2) <= 0  % Rs must be positive
            fprintf('Warning: Rs in final fitting result is negative or zero, adjusted to positive value\n');
            optimized_params(2) = max(params.lb(2), 10);  % Ensure greater than zero
            
            % Recalculate fitting results
            fit_results.JD = diodeModel(data_V, optimized_params, config);
            relative_errors = abs((fit_results.JD - data_JD) ./ (abs(data_JD) + eps)) * 100;
        end
        
        % Check if there is a significant gap between positive and negative regions
        fprintf('\n检查正负电压区域拟合效果差异...\n');
        neg_errors = relative_errors(neg_idx);
        pos_errors = relative_errors(pos_idx);
        fprintf('负电压区域平均相对误差: %.4f%%\n', mean(neg_errors));
        fprintf('正电压区域平均相对误差: %.4f%%\n', mean(pos_errors));
        
        % If positive region fitting effect is significantly worse than negative region, try to optimize positive region separately
        %%
        high_pos_idx = find(data_V > 0.20);
        if ~isempty(high_pos_idx)
            high_pos_errors = relative_errors(high_pos_idx);
            fprintf('高正偏压区域(>0.20V)平均相对误差: %.4f%%\n', mean(high_pos_errors));
            
            if mean(high_pos_errors) > 5  % 如果高正偏压区域误差大于15%
                fprintf('\n========== 开始执行: 高正偏压区域(>0.20V)专门优化 ==========\n');
                fprintf('高正偏压区域误差大于5%%，需要专门优化\n');
                
                % 选择高正偏压区域数据
                high_pos_V = data_V(high_pos_idx);
                high_pos_JD = data_JD(high_pos_idx);
                
                % 只优化Rs和J0
                param_mask = [true, true, false, false];
                
                % 创建专门的高正偏压误差函数
                fprintf('创建高正偏压区域专用误差函数...\n');
                high_pos_errFun = @(x_opt) errorFunctionStrongPositive(x_opt, optimized_params ./ params.scaleFactors, param_mask, data_V, data_JD, params, config);
                fprintf('高正偏压区域专用误差函数创建完成\n');
                
                % 准备初始参数
                x0_high_opt = optimized_params(param_mask) ./ params.scaleFactors(param_mask);
                lb_high = params.lb(param_mask) ./ params.scaleFactors(param_mask);
                ub_high = params.ub(param_mask) ./ params.scaleFactors(param_mask);
                
                % 执行优化
                options_high = optimoptions('lsqnonlin', ...
                    'Display', 'iter-detailed', ...
                    'Algorithm', 'trust-region-reflective', ...
                    'FunctionTolerance', 1e-12, ...
                    'OptimalityTolerance', 1e-12, ...
                    'StepTolerance', 1e-10, ...
                    'MaxFunctionEvaluations', 5000, ...
                    'MaxIterations', 3000);
                
                fprintf('开始执行高正偏压区域优化...\n');
                [x_high_opt, resnorm_high, ~, exitflag_high, output_high] = lsqnonlin(high_pos_errFun, x0_high_opt, lb_high, ub_high, options_high);
                fprintf('高正偏压区域优化完成\n');
                fprintf('优化退出标志: %d, 迭代次数: %d, 函数评估次数: %d\n', exitflag_high, output_high.iterations, output_high.funcCount);
                
                % 更新参数
                fprintf('更新参数...\n');
                x_scaled_high = optimized_params ./ params.scaleFactors;
                x_scaled_high(param_mask) = x_high_opt;
                optimized_params_high = x_scaled_high .* params.scaleFactors;
                fprintf('参数更新完成\n');
                
                % 计算优化后的拟合结果
                fprintf('计算高正偏压区域优化后的拟合结果...\n');
                fit_results_high.JD = diodeModel(data_V, optimized_params_high, config);
                relative_errors_high = abs((fit_results_high.JD - data_JD) ./ (abs(data_JD) + eps)) * 100;
                fprintf('拟合结果计算完成\n');
                
                % 检查优化效果
                high_errors_before = high_pos_errors;
                high_errors_after = relative_errors_high(high_pos_idx);
                
                fprintf('高正偏压优化前: 平均相对误差 = %.2f%%\n', mean(high_errors_before));
                fprintf('高正偏压优化后: 平均相对误差 = %.2f%%\n', mean(high_errors_after));
                
                % 如果优化确实改善了高正偏压区域且没有显著损害其他区域
                if mean(high_errors_after) < mean(high_errors_before) * 0.8
                    fprintf('使用高正偏压优化结果\n');
                    optimized_params = optimized_params_high;
                    fit_results.JD = fit_results_high.JD;
                    relative_errors = relative_errors_high;
                else
                    fprintf('保持原始优化结果\n');
                end
            end
        end
        %%
        fprintf('\n检查是否需要执行正偏压区域单独优化...\n');
        fprintf('正偏压区域误差与负偏压区域误差比值: %.4f\n', mean(pos_errors)/mean(neg_errors));
        
        if mean(pos_errors) > 2*mean(neg_errors) && mean(pos_errors) > 5
            fprintf('\n========== 开始执行: 正偏压区域单独优化 ==========\n');
            fprintf('正偏压区域误差大于负偏压区域误差的2倍且大于10%%，需要单独优化\n');
            
            % Define optimization parameter indices - only optimize J0 and Rs
            param_mask = [true, true, false, false]; % Only optimize J0 and Rs
            
            % Create error function with enhanced weights for positive region
            fprintf('创建正偏压区域增强误差函数...\n');
            pos_errFun = @(x_opt) errorFunctionEnhancedPositive(x_opt, optimized_params ./ params.scaleFactors, param_mask, data_V, data_JD, params, config);
            fprintf('正偏压区域增强误差函数创建完成\n');
            
            % Prepare initial parameters containing only J0 and Rs
            x0_pos_opt = optimized_params(param_mask) ./ params.scaleFactors(param_mask);
            lb_pos = params.lb(param_mask) ./ params.scaleFactors(param_mask);
            ub_pos = params.ub(param_mask) ./ params.scaleFactors(param_mask);
            
            % Execute optimization
            % options_pos = optimoptions('lsqnonlin', ...
            %     'Display', 'iter-detailed', ...
            %     'Algorithm', 'levenberg-marquardt', ...
            %     'FunctionTolerance', 1e-12, ...
            %     'OptimalityTolerance', 1e-12, ...
            %     'StepTolerance', 1e-12, ...
            %     'MaxFunctionEvaluations', 3000, ...
            %     'MaxIterations', 2000);
            options_pos = optimoptions('lsqnonlin', ...
                'Display', 'iter-detailed', ...
                'Algorithm', 'levenberg-marquardt', ...
                'FunctionTolerance', 1e-10, ...
                'OptimalityTolerance', 1e-10, ...
                'StepTolerance', 1e-10, ...
                'MaxFunctionEvaluations', 5000, ...
                'MaxIterations', 2000);

            fprintf('开始执行正偏压区域单独优化...\n');
            [x_pos_opt, resnorm_pos, ~, exitflag_pos, output_pos] = lsqnonlin(pos_errFun, x0_pos_opt, [], [], options_pos);
            
            fprintf('正偏压区域单独优化完成\n');
            fprintf('优化退出标志: %d, 迭代次数: %d, 函数评估次数: %d\n', exitflag_pos, output_pos.iterations, output_pos.funcCount);
            
            % Validate Rs as positive
            if x_pos_opt(2) * params.scaleFactors(2) <= 0
                fprintf('Warning: Positive region optimization produced negative or zero Rs, adjusting to positive value\n');
                % Use positive value as alternative
                x_pos_opt(2) = lb_pos(2);
            end
            
            % Update optimization results back to complete parameter set
            fprintf('更新参数...\n');
            x_scaled_enhanced = optimized_params ./ params.scaleFactors;
            x_scaled_enhanced(param_mask) = x_pos_opt;
            optimized_params_enhanced = x_scaled_enhanced .* params.scaleFactors;
            fprintf('参数更新完成\n');
            
            % Compute enhanced fitting results
            fprintf('计算正偏压区域优化后的拟合结果...\n');
            fit_results_enhanced.JD = diodeModel(data_V, optimized_params_enhanced, config);
            relative_errors_enhanced = abs((fit_results_enhanced.JD - data_JD) ./ (abs(data_JD) + eps)) * 100;
            fprintf('拟合结果计算完成\n');
            
            % Check enhanced fitting effects
            pos_errors_enhanced = relative_errors_enhanced(pos_idx);
            neg_errors_enhanced = relative_errors_enhanced(neg_idx);
            
            fprintf('Enhanced before: Positive region error = %.2f%%, Negative region error = %.2f%%\n', ...
                mean(pos_errors), mean(neg_errors));
            fprintf('Enhanced after: Positive region error = %.2f%%, Negative region error = %.2f%%\n', ...
                mean(pos_errors_enhanced), mean(neg_errors_enhanced));
            
            % If enhanced really improved positive region while not damaging negative region
            if mean(pos_errors_enhanced) < mean(pos_errors) && mean(neg_errors_enhanced) < 2*mean(neg_errors)
                fprintf('Using enhanced optimized results\n');
                optimized_params = optimized_params_enhanced;
                fit_results.JD = fit_results_enhanced.JD;
                relative_errors = relative_errors_enhanced;
            else
                fprintf('Keep original optimization result\n');
            end
        end
        
        % Compute negative voltage region error
        neg_idx = find(data_V < -0.1);
        neg_errors = relative_errors(neg_idx);
        
        % Output error statistics
        fprintf('\nRelative error statistics for each point:\n');
        fprintf('Maximum relative error: %.2f%%\n', max(relative_errors));
        fprintf('Average relative error: %.2f%%\n', mean(relative_errors));
        fprintf('Median relative error: %.2f%%\n', median(relative_errors));
        fprintf('Negative voltage region average relative error: %.2f%%\n', mean(neg_errors));
        
        % Special check for strong negative voltage region (-0.5 to -0.2)
        %
        % strong_neg_idx = find(data_V < -0.2 & data_V >= -0.5);
        %
        strong_neg_idx = find(data_V < -0.4 & data_V >= -0.5);
        if ~isempty(strong_neg_idx)
            strong_neg_errors = relative_errors(strong_neg_idx);
            fprintf('Strong negative voltage region (-0.5 to -0.2V) average relative error: %.2f%%\n', mean(strong_neg_errors));
        end
        
        % Special check for positive voltage region
        pos_idx = find(data_V > 0);
        if ~isempty(pos_idx)
            pos_errors = relative_errors(pos_idx);
            fprintf('Positive voltage region average relative error: %.2f%%\n', mean(pos_errors));
            
            % Further subdivide positive voltage region
            low_pos_idx = find(data_V > 0 & data_V <= 0.15);
            %mid_pos_idx = find(data_V > 0.15 & data_V <= 0.25);
            mid_pos_idx = find(data_V > 0.15 & data_V <= 0.20);
            high_pos_idx = find(data_V > 0.20);
            %high_pos_idx = find(data_V > 0.25);
            
            if ~isempty(low_pos_idx)
                low_pos_errors = relative_errors(low_pos_idx);
                fprintf('Low positive voltage region (0-0.15V) average relative error: %.2f%%\n', mean(low_pos_errors));
            end
            
            if ~isempty(mid_pos_idx)
                mid_pos_errors = relative_errors(mid_pos_idx);
                fprintf('Medium positive voltage region (0.15-0.25V) average relative error: %.2f%%\n', mean(mid_pos_errors));
            end
            
            if ~isempty(high_pos_idx)
                high_pos_errors = relative_errors(high_pos_idx);
                fprintf('High positive voltage region (>0.25V) average relative error: %.2f%%\n', mean(high_pos_errors));
            end
        end
        
        % Find error of the point with the maximum error
        [max_error, max_error_idx] = max(relative_errors);
        fprintf('\nPoint with maximum error:\n');
        fprintf('Voltage: %.3f V\n', data_V(max_error_idx));
        fprintf('Measured current: %.3e A\n', data_JD(max_error_idx));
        fprintf('Fitted current: %.3e A\n', fit_results.JD(max_error_idx));
        fprintf('Relative error: %.2f%%\n', max_error);
        
        % Output final fitting parameters
        fprintf('\nFitting parameters:\n');
        fprintf('J0 = %.6e A\n', optimized_params(1));
        fprintf('Rs = %.6e Ohm\n', optimized_params(2));
        fprintf('Rsh = %.6e Ohm\n', optimized_params(3));
        fprintf('k = %.6e\n', optimized_params(4));
        
    catch ME
        error('Error in fitting process: %s', ME.message);
    end
end

% Error function for negative voltage region
function err = errorFunctionNegative(x, data_V, data_JD, params, config)
    % Unscale parameters
    x_actual = x .* params.scaleFactors;
    
    % Calculate model predictions
    predicted = diodeModel(data_V, x_actual, config);
    
    % Calculate errors
    err = zeros(size(data_JD));
    
    for i = 1:length(data_JD)
        actual_abs = abs(data_JD(i));
        pred_abs = abs(predicted(i));
        
        threshold = 1e-12;
        
        if actual_abs < threshold || pred_abs < threshold
            err(i) = (predicted(i) - data_JD(i)) / max(1e-12, max(max(abs(data_JD))));
        else
            % Use logarithmic error
            log_actual = log10(actual_abs);
            log_pred = log10(pred_abs);
            err(i) = log_pred - log_actual;
            
            % Maintain sign consistency
            if sign(predicted(i)) ~= sign(data_JD(i))
                err(i) = err(i) * 4;
            end
        end
        
        % Add weight for more negative voltage regions
        
        if data_V(i) < -0.3
            err(i) = err(i) * 3;
        end
    end
end

% 优化 errorFunctionPartial 函数，提高拟合精度
function err = errorFunctionPartial(x_opt, x0, param_mask, data_V, data_JD, params, config)
    % Build complete parameter vector
    fprintf('执行errorFunctionPartial函数，构建完整参数向量...\n');
    x_full = x0;
    x_full(param_mask) = x_opt;
    
    % Unscale parameters
    fprintf('解除参数缩放...\n');
    x_actual = x_full .* params.scaleFactors;
    
    % Calculate model predictions
    fprintf('计算模型预测值...\n');
    predicted = diodeModel(data_V, x_actual, config);
    
    % Create voltage masks
    fprintf('创建电压区域掩码...\n');
    masks = createVoltageMasks(data_V, data_JD, predicted);
    
    % Calculate weighted error
    err = calculateWeightedError(data_V, data_JD, predicted, masks, 'all');
end
% Partial parameter optimization error function
% function err = errorFunctionPartial(x_opt, x0, param_mask, data_V, data_JD, params, config)
%     % Build complete parameter vector
%     x_full = x0;
%     x_full(param_mask) = x_opt;
    
%     % Unscale parameters
%     x_actual = x_full .* params.scaleFactors;
    
%     % Calculate model predictions
%     predicted = diodeModel(data_V, x_actual, config);
    
%     % Calculate errors
%     err = zeros(size(data_JD));
    
%     for i = 1:length(data_JD)
%         actual_abs = abs(data_JD(i));
%         pred_abs = abs(predicted(i));
        
%         threshold = 1e-12;
        
%         if actual_abs < threshold || pred_abs < threshold
%             err(i) = (predicted(i) - data_JD(i)) / max(1e-12, max(max(abs(data_JD))));
%         else
%             % Use logarithmic error
%             log_actual = log10(actual_abs);
%             log_pred = log10(pred_abs);
%             err(i) = log_pred - log_actual;
            
%             % Maintain sign consistency
%             if sign(predicted(i)) ~= sign(data_JD(i))
%                 err(i) = err(i) * 4;
%             end
%         end
        
%         % Apply different weights for different voltage regions
%         if data_V(i) < -0.3
%             % err(i) = err(i) * 5; % Strong negative voltage region
%             err(i) = err(i) * 8; % Strong negative voltage region
%         elseif data_V(i) < -0.1
%             err(i) = err(i) * 5; % Weak negative voltage region
%             % err(i) = err(i) * 3;
%         elseif data_V(i) < 0.1
%             err(i) = err(i) * 2; % Near zero point
%         end
%     end
% end

% 优化 errorFunctionEnhancedPositive 函数，提高高正偏压区域拟合效果
function err = errorFunctionEnhancedPositive(x_opt, x0, param_mask, data_V, data_JD, params, config)
    % Build complete parameter vector
    x_full = x0;
    x_full(param_mask) = x_opt;
    
    % Unscale parameters
    x_actual = x_full .* params.scaleFactors;
    
    % Calculate model predictions
    predicted = diodeModel(data_V, x_actual, config);
    
    % Create voltage masks
    masks = createVoltageMasks(data_V, data_JD, predicted);
    
    % Calculate weighted error with enhanced positive region weights
    err = calculateWeightedError(data_V, data_JD, predicted, masks, 'enhanced_positive');
end
% Enhanced positive voltage region fitting error function
% function err = errorFunctionEnhancedPositive(x_opt, x0, param_mask, data_V, data_JD, params, config)
%     % Build complete parameter vector
%     x_full = x0;
%     x_full(param_mask) = x_opt;
    
%     % Unscale parameters
%     x_actual = x_full .* params.scaleFactors;
    
%     % Calculate model predictions
%     predicted = diodeModel(data_V, x_actual, config);
    
%     % Calculate errors
%     err = zeros(size(data_JD));
    
%     for i = 1:length(data_JD)
%         actual_abs = abs(data_JD(i));
%         pred_abs = abs(predicted(i));
        
%         threshold = 1e-12;
        
%         if actual_abs < threshold || pred_abs < threshold
%             err(i) = (predicted(i) - data_JD(i)) / max(1e-12, max(max(abs(data_JD))));
%         else
%             % Use logarithmic error
%             log_actual = log10(actual_abs);
%             log_pred = log10(pred_abs);
%             err(i) = log_pred - log_actual;
            
%             % Maintain sign consistency
%             if sign(predicted(i)) ~= sign(data_JD(i))
%                 err(i) = err(i) * 4;
%             end
%         end
%         % Specially enhance positive voltage region weights, especially for poorly fitted high voltage regions
%         % if data_V(i) > 0.25
%         if data_V(i) > 0.20    
%             err(i) = err(i) * 30;  % Especially emphasize high positive voltage region
%             % err(i) = err(i) * 15;
%         elseif data_V(i) > 0.15
%             err(i) = err(i) * 5;  % Medium positive voltage region
%         elseif data_V(i) > 0
%             err(i) = err(i) * 3;  % Low positive voltage region
%         elseif data_V(i) < -0.3
%             err(i) = err(i) * 2; % Weaken the influence of strong negative voltage region
%             % err(i) = err(i) * 0.5;
%         end
%     end
% end

% 优化 errorFunction 函数，提高拟合精度
function err = errorFunction(x, data_V, data_JD, params, config)
    % Unscale parameters
    x_actual = x .* params.scaleFactors;
    
    % Check if Rs is negative, if so automatically adjust to positive value
    if x_actual(2) <= 0
        % Output warning instead of throwing error
        warning('Rs is negative (%.2e), automatically adjusting to positive value', x_actual(2));
        
        % Adjust Rs to a small positive value (use lower bound or 10 ohms, whichever is larger)
        x_actual(2) = max(params.lb(2), 10);
        
        % Update scaled parameter
        x(2) = x_actual(2) / params.scaleFactors(2);
    end
    
    % Calculate model predictions
    predicted = diodeModel(data_V, x_actual, config);
    
    % 使用向量化操作提高计算效率
    actual_abs = abs(data_JD);
    pred_abs = abs(predicted);
    
    threshold = 1e-12;
    
    % 初始化误差向量
    err = zeros(size(data_JD));
    
    % 处理小于阈值的情况
    small_vals = (actual_abs < threshold) | (pred_abs < threshold);
    err(small_vals) = (predicted(small_vals) - data_JD(small_vals)) / max(1e-12, max(max(abs(data_JD))));
    
    % 处理正常值的情况
    normal_vals = ~small_vals;
    if any(normal_vals)
        % 使用对数误差
        log_actual = log10(actual_abs(normal_vals));
        log_pred = log10(pred_abs(normal_vals));
        err(normal_vals) = log_pred - log_actual;
        
        % 符号不一致的惩罚
        sign_mismatch = sign(predicted(normal_vals)) ~= sign(data_JD(normal_vals));
        err(normal_vals(sign_mismatch)) = err(normal_vals(sign_mismatch)) * 5;  % 增加符号不匹配惩罚
    end
    
    % 为不同电压区域应用权重 - 平衡各区域拟合效果
    high_pos = data_V > 0.20;
    mid_pos = data_V > 0.15 & data_V <= 0.20;
    low_pos = data_V > 0 & data_V <= 0.15;
    strong_neg = data_V < -0.3;
    
    err(high_pos) = err(high_pos) * 2.5;  % 增加高正偏压区域权重
    err(mid_pos) = err(mid_pos) * 2.0;    % 增加中正偏压区域权重
    err(low_pos) = err(low_pos) * 1.5;    % 增加低正偏压区域权重
    err(strong_neg) = err(strong_neg) * 2.0; % 增加强负电压区域权重
    
    % 如果Rs接近下界，添加惩罚以防止优化器将其推向更小的值
    if x_actual(2) < params.lb(2) * 1.2
        err = err * (1 + 0.15 * (params.lb(2) * 1.2 - x_actual(2)) / params.lb(2));  % 增加惩罚系数
    end
end
% Standard error function - modified to automatically handle negative Rs values
% function err = errorFunction(x, data_V, data_JD, params, config)
%     % Unscale parameters
%     x_actual = x .* params.scaleFactors;
    
%     % Check if Rs is negative, if so automatically adjust to positive value
%     if x_actual(2) <= 0
%         % Output warning instead of throwing error
%         warning('Rs is negative (%.2e), automatically adjusting to positive value', x_actual(2));
        
%         % Adjust Rs to a small positive value (use lower bound or 10 ohms, whichever is larger)
%         x_actual(2) = max(params.lb(2), 10);
        
%         % Update scaled parameter
%         x(2) = x_actual(2) / params.scaleFactors(2);
%     end
    
%     % Calculate model predictions
%     predicted = diodeModel(data_V, x_actual, config);
    
%     % Calculate errors
%     err = zeros(size(data_JD));
    
%     for i = 1:length(data_JD)
%         actual_abs = abs(data_JD(i));
%         pred_abs = abs(predicted(i));
        
%         threshold = 1e-12;
        
%         if actual_abs < threshold || pred_abs < threshold
%             err(i) = (predicted(i) - data_JD(i)) / max(1e-12, max(max(abs(data_JD))));
%         else
%             % Use logarithmic error
%             log_actual = log10(actual_abs);
%             log_pred = log10(pred_abs);
%             err(i) = log_pred - log_actual;
            
%             % Maintain sign consistency
%             if sign(predicted(i)) ~= sign(data_JD(i))
%                 err(i) = err(i) * 4;
%             end
%         end
%     end
    
%     % If Rs is close to lower bound, add penalty to prevent optimizer from pushing it to smaller values
%     if x_actual(2) < params.lb(2) * 1.1
%         err = err * (1 + 0.1 * (params.lb(2) * 1.1 - x_actual(2)) / params.lb(2));
%     end
% end

% 优化 errorFunctionStrongPositive 函数，显著提高高正偏压区域拟合效果
function err = errorFunctionStrongPositive(x_opt, x0, param_mask, data_V, data_JD, params, config)
    % 构建完整参数向量
    x_full = x0;
    x_full(param_mask) = x_opt;
    
    % 解除参数缩放
    x_actual = x_full .* params.scaleFactors;
    
    % 计算模型预测
    predicted = diodeModel(data_V, x_actual, config);
    
    % 使用向量化操作提高计算效率
    threshold = 1e-12;
    actual_abs = abs(data_JD);
    pred_abs = abs(predicted);
    
    % 初始化误差向量
    err = zeros(size(data_JD));
    
    % 处理低于阈值的值
    below_threshold_mask = (actual_abs < threshold) | (pred_abs < threshold);
    err(below_threshold_mask) = (predicted(below_threshold_mask) - data_JD(below_threshold_mask)) / max(1e-12, max(abs(data_JD)));
    
    % 计算对数空间误差
    above_threshold_mask = ~below_threshold_mask;
    if any(above_threshold_mask)
        log_actual = log10(actual_abs(above_threshold_mask));
        log_pred = log10(pred_abs(above_threshold_mask));
        err(above_threshold_mask) = log_pred - log_actual;
    end
    
    % 符号不匹配的惩罚
    sign_mismatch_mask = sign(predicted) ~= sign(data_JD) & actual_abs > threshold & pred_abs > threshold;
    err(sign_mismatch_mask) = err(sign_mismatch_mask) * 6;  % 增加符号不匹配惩罚
    
    % 为不同电压区域应用权重 - 特别强调高正偏压区域
    low_pos_mask = data_V >= 0 & data_V < 0.15;
    mid_pos_mask = data_V >= 0.15 & data_V < 0.20;
    high_pos_mask = data_V >= 0.20;
    
    % 增加权重以提高拟合效果
    err(low_pos_mask) = err(low_pos_mask) * 8;       % 增加低正偏压区域权重
    err(mid_pos_mask) = err(mid_pos_mask) * 25;      % 显著增加中正偏压区域权重
    err(high_pos_mask) = err(high_pos_mask) * 50;    % 显著增加高正偏压权重
    
    % 降低负电压区域的权重
    neg_mask = data_V < 0;
    err(neg_mask) = err(neg_mask) * 0.05;  % 进一步降低负电压区域权重
    
    % 添加额外的高正偏压区域细分，更精细地控制拟合
    very_high_pos_mask = data_V >= 0.25;
    if any(very_high_pos_mask)
        err(very_high_pos_mask) = err(very_high_pos_mask) * 1.5;  % 对最高正偏压区域额外加权
    end
end

function masks = createVoltageMasks(data_V, data_JD, predicted)
    % 创建不同电压区域的掩码
    masks.high_pos = data_V > 0.20;
    masks.mid_pos = data_V > 0.15 & data_V <= 0.20;
    masks.low_pos = data_V > 0 & data_V <= 0.15;
    masks.neg = data_V < 0;
    masks.strong_neg = data_V < -0.3;
    % 添加其他需要的掩码
end

function err = calculateWeightedError(data_V, data_JD, predicted, masks, mode)
    % 计算加权误差
    actual_abs = abs(data_JD);
    pred_abs = abs(predicted);
    threshold = 1e-12;
    
    % 初始化误差向量
    err = zeros(size(data_JD));
    
    % 处理小于阈值的情况
    small_vals = (actual_abs < threshold) | (pred_abs < threshold);
    err(small_vals) = (predicted(small_vals) - data_JD(small_vals)) / max(1e-12, max(max(abs(data_JD))));
    
    % 处理正常值的情况
    normal_vals = ~small_vals;
    if any(normal_vals)
        % 使用对数误差
        log_actual = log10(actual_abs(normal_vals));
        log_pred = log10(pred_abs(normal_vals));
        err(normal_vals) = log_pred - log_actual;
        
        % 符号不一致的惩罚
        sign_mismatch = sign(predicted(normal_vals)) ~= sign(data_JD(normal_vals));
        err(normal_vals(sign_mismatch)) = err(normal_vals(sign_mismatch)) * 5;
    end
    
    % 根据模式应用不同权重
    if strcmp(mode, 'all')
        % 均衡权重
        err(masks.high_pos) = err(masks.high_pos) * 2.5;
        err(masks.mid_pos) = err(masks.mid_pos) * 2.0;
        err(masks.low_pos) = err(masks.low_pos) * 1.5;
        err(masks.strong_neg) = err(masks.strong_neg) * 2.0;
    elseif strcmp(mode, 'enhanced_positive')
        % 增强正偏压权重
        err(masks.high_pos) = err(masks.high_pos) * 50;
        err(masks.mid_pos) = err(masks.mid_pos) * 25;
        err(masks.low_pos) = err(masks.low_pos) * 8;
        err(masks.neg) = err(masks.neg) * 0.05;
    end
end
