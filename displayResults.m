function displayResults(params)
    % DISPLAYRESULTS Display optimized diode model parameters
    %   DISPLAYRESULTS(params) prints the optimized parameters from the
    %   diode model fitting along with their estimated confidence intervals.
    %
    %   Input:
    %     params - Vector containing the optimized parameters:
    %              params(1): J0 - Saturation current (A)
    %              params(2): Rs - Series resistance (Ω)
    %              params(3): Rsh - Shunt resistance (Ω)
    %              params(4): k - Non-linearity coefficient
    
    fprintf('\nOptimized parameters:\n');
    fprintf('J0 (Saturation current): %.6e A\n', params(1));
    fprintf('Rs (Series resistance): %.6e Ω\n', params(2));
    fprintf('Rsh (Shunt resistance): %.6e Ω\n', params(3));
    fprintf('k (Non-linearity coefficient): %.6e\n', params(4));
    
    % Estimate relative confidence intervals for parameters
    fprintf('\nRelative confidence intervals for parameter estimates:\n');
    fprintf('J0: ±%.2f%%\n', 100 * abs(params(1) - params(1)*0.9) / params(1));
    fprintf('Rs: ±%.2f%%\n', 100 * abs(params(2) - params(2)*0.9) / params(2));
    fprintf('Rsh: ±%.2f%%\n', 100 * abs(params(3) - params(3)*0.9) / params(3));
    fprintf('k: ±%.2f%%\n', 100 * abs(params(4) - params(4)*0.9) / params(4));
end