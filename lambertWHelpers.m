function [solveForJ0Func, calcCurrentFunc, approxWFunc] = lambertWHelpers()
    % LAMBERTWHELPERS 返回与Lambert W函数相关的辅助函数句柄
    %   [SOLVEFUNJ0FUNC, CALCCURRENTFUNC, APPROXWFUNC] = LAMBERTWHELPERS() 返回
    %   三个函数句柄，用于二极管模型中的Lambert W函数计算。
    %
    %   输出:
    %     SOLVEFUNJ0FUNC - 用于求解J0的目标函数句柄
    %     CALCCURRENTFUNC - 使用Lambert W函数计算电流的函数句柄
    %     APPROXWFUNC - Lambert W函数的近似实现函数句柄

    % 返回函数句柄
    solveForJ0Func = @solveForJ0;
    calcCurrentFunc = @calculateCurrentWithLambertW;
    approxWFunc = @approximateLambertW;
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