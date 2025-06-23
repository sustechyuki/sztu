function saveResults(data_V, data_JD, params, fit_results, currents)
    % SAVERESULTS 保存拟合结果
    %   SAVERESULTS(DATA_V, DATA_JD, PARAMS, FIT_RESULTS, CURRENTS) 将拟合结果
    %   保存到多种格式的文件中，包括MAT、CSV、TXT和PNG。
    %
    %   输入:
    %     DATA_V      - 电压数据向量
    %     DATA_JD     - 电流密度数据向量
    %     PARAMS      - 拟合参数向量 [J0, Rs, Rsh, k]
    %     FIT_RESULTS - 拟合结果结构体
    %     CURRENTS    - 电流分量结构体

    % Generate timestamp
    timestamp = datestr(now, 'yyyymmdd_HHMMSS');
    
    % Save data and fitting results
    results_filename = sprintf('fit_results_%s.mat', timestamp);
    save(results_filename, 'data_V', 'data_JD', 'params', 'fit_results', 'currents');
    fprintf('Fitting results saved to file: %s\n', results_filename);
    
    % Save figure
    fig_filename = sprintf('fit_plot_%s.png', timestamp);
    saveas(gcf, fig_filename);
    fprintf('Fitting plot saved to file: %s\n', fig_filename);
    
    % Export detailed data to CSV
    csv_filename = sprintf('fit_data_%s.csv', timestamp);
    fid = fopen(csv_filename, 'w');
    fprintf(fid, 'Voltage(V),Measured_Current(A),Fitted_Current(A),Diode_Current(A),Ohmic_Current(A),Nonohmic_Current(A),Relative_Error(%%)\n');
    
    % Calculate relative errors
    rel_errors = abs((fit_results.JD - data_JD) ./ (abs(data_JD) + eps)) * 100;
    
    % Write data
    for i = 1:length(data_V)
        fprintf(fid, '%.6e,%.6e,%.6e,%.6e,%.6e,%.6e,%.2f\n', ...
            data_V(i), data_JD(i), fit_results.JD(i), ...
            currents.diode(i), currents.ohmic(i), currents.nonohmic(i), ...
            rel_errors(i));
    end
    fclose(fid);
    fprintf('Fitting data exported to CSV file: %s\n', csv_filename);
    
    % Export fitting parameters to text file
    params_filename = sprintf('fit_params_%s.txt', timestamp);
    fid = fopen(params_filename, 'w');
    fprintf(fid, 'Fitting Parameters:\n');
    fprintf(fid, 'J0 = %.6e A\n', params(1));
    fprintf(fid, 'Rs = %.6e Ohm\n', params(2));
    fprintf(fid, 'Rsh = %.6e Ohm\n', params(3));
    fprintf(fid, 'k = %.6e\n', params(4));
    fprintf(fid, '\n');
    
    % Fitting error statistics
    fprintf(fid, 'Fitting Error Statistics:\n');
    fprintf(fid, 'Average Relative Error: %.2f%%\n', mean(rel_errors));
    fprintf(fid, 'Maximum Relative Error: %.2f%%\n', max(rel_errors));
    fprintf(fid, 'Median Relative Error: %.2f%%\n', median(rel_errors));
    
    % Calculate errors for different voltage regions
    neg_idx = find(data_V < 0);
    pos_idx = find(data_V > 0);
    if ~isempty(neg_idx)
        fprintf(fid, 'Average Relative Error in Negative Voltage Region: %.2f%%\n', mean(rel_errors(neg_idx)));
    end
    if ~isempty(pos_idx)
        fprintf(fid, 'Average Relative Error in Positive Voltage Region: %.2f%%\n', mean(rel_errors(pos_idx)));
    end
    
    fclose(fid);
    fprintf('Fitting parameters and statistics saved to file: %s\n', params_filename);
end