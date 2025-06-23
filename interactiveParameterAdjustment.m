function [adjusted_params, fit_results] = interactiveParameterAdjustment(data_V, data_JD, initial_params, config)
    % INTERACTIVEPARAMETERADJUSTMENT 交互式参数调整
    %   [ADJUSTED_PARAMS, FIT_RESULTS] = INTERACTIVEPARAMETERADJUSTMENT(DATA_V, DATA_JD, INITIAL_PARAMS, CONFIG)
    %   提供交互式界面，允许用户手动调整拟合参数并实时查看拟合效果。
    %
    %   输入:
    %     DATA_V         - 电压数据向量
    %     DATA_JD        - 电流密度数据向量
    %     INITIAL_PARAMS - 初始参数向量 [J0, Rs, Rsh, k]
    %     CONFIG         - 配置结构体
    %
    %   输出:
    %     ADJUSTED_PARAMS - 调整后的参数向量
    %     FIT_RESULTS     - 拟合结果结构体

    % Copy initial parameters
    adjusted_params = initial_params;
    
    % Initialize parameter history for undo function
    params_history = initial_params;
    history_index = 1;
    
    % Store the original initial parameters (for exit without saving)
    original_params = initial_params;
    
    % Ensure physical reasonability of initial parameters
    if adjusted_params(2) <= 0  % Rs must be positive
        fprintf('Warning: Initial Rs is negative or zero, automatically adjusted to positive value\n');
        adjusted_params(2) = 10; % Use a reasonable default value
        original_params = adjusted_params; % Update original params as well
    end
    
    % Calculate initial fit and errors
    fit_results.JD = diodeModel(data_V, adjusted_params, config);
    errors = abs((fit_results.JD - data_JD) ./ (abs(data_JD) + eps)) * 100;
    
    % Calculate current components for CSV export
    currents = calculateCurrents(data_V, adjusted_params, config);
    
    % Create real-time updated chart
    fig_handle = figure('Name', 'Interactive Parameter Adjustment', 'Position', [100, 100, 1200, 800]);
    
    % Define subplot structure
    subplot(2,1,1);
    h_data = semilogy(data_V, abs(data_JD), 'bo', 'DisplayName', 'Measured Data');
    hold on;
    h_fit = semilogy(data_V, abs(fit_results.JD), 'ro', 'DisplayName', 'Fitted Result');
    xlabel('Voltage (V)');
    ylabel('Current Density (A)');
    title('Current-Voltage Characteristics (Log Scale)');
    legend('Location', 'best');
    grid on;
    
    subplot(2,1,2);
    h_error = plot(data_V, errors, 'b.-');
    xlabel('Voltage (V)');
    ylabel('Relative Error (%)');
    title(sprintf('Fitting Error (Average: %.2f%%)', mean(errors)));
    grid on;
    
    % Display current parameter values
    annotation('textbox', [0.01, 0.01, 0.98, 0.08], ...
        'String', sprintf('J0: %.2e A   Rs: %.2e Ohm   Rsh: %.2e Ohm   k: %.2e', ...
        adjusted_params(1), adjusted_params(2), adjusted_params(3), adjusted_params(4)), ...
        'EdgeColor', 'none', 'FontSize', 10, 'HorizontalAlignment', 'center');
    
    % Initialize adaptive step factors
    adapt_factors = [0.1, 0.1, 0.1, 0.1]; % For J0, Rs, Rsh, k
    
    % Flag to track if user wants to save adjusted parameters
    save_adjustments = true;
    
    % Continue adjusting until user is satisfied
    while true
        % Display adjustment options
        fprintf('\nCurrent parameters: J0=%.2e, Rs=%.2e, Rsh=%.2e, k=%.2e\n', ...
            adjusted_params(1), adjusted_params(2), adjusted_params(3), adjusted_params(4));
        fprintf('Average relative error: %.2f%%\n', mean(errors));
        
        % Display optimization recommendations
        displayOptimizationRecommendations(adjusted_params, errors, data_V, data_JD);
        
        fprintf('\nParameter adjustment options:\n');
        fprintf('1: Increase J0  2: Decrease J0\n');
        fprintf('3: Increase Rs  4: Decrease Rs\n');
        fprintf('5: Increase Rsh 6: Decrease Rsh\n');
        fprintf('7: Increase k   8: Decrease k\n');
        fprintf('9: Undo last adjustment\n');
        fprintf('10: Save current components to CSV file\n');
        fprintf('0: Exit without saving\n');
        
        % Get user input and ensure it's a numeric type
        choice_str = input('Please select operation (0-10): ', 's');
        choice = str2double(choice_str);
        
        % Check if it's a valid numeric input
        if isnan(choice)
            fprintf('Please enter a valid number (0-10)\n');
            continue;
        end
        
        if choice == 0
            % Exit without saving changes
            fprintf('Exiting without saving parameter adjustments.\n');
            save_adjustments = false;
            break;
        elseif choice == 9
            % Undo function
            if history_index > 1
                % Return to previous parameter set
                history_index = history_index - 1;
                adjusted_params = params_history(history_index, :);
                
                % Recalculate fit and errors
                fit_results.JD = diodeModel(data_V, adjusted_params, config);
                errors = abs((fit_results.JD - data_JD) ./ (abs(data_JD) + eps)) * 100;
                currents = calculateCurrents(data_V, adjusted_params, config);
                
                % Update chart
                set(h_fit, 'YData', abs(fit_results.JD));
                set(h_error, 'YData', errors);
                title(subplot(2,1,2), sprintf('Fitting Error (Average: %.2f%%)', mean(errors)));
                
                % Update parameter display
                delete(findall(gcf, 'Type', 'annotation'));
                annotation('textbox', [0.01, 0.01, 0.98, 0.08], ...
                    'String', sprintf('J0: %.2e A   Rs: %.2e Ohm   Rsh: %.2e Ohm   k: %.2e', ...
                    adjusted_params(1), adjusted_params(2), adjusted_params(3), adjusted_params(4)), ...
                    'EdgeColor', 'none', 'FontSize', 10, 'HorizontalAlignment', 'center');
                
                fprintf('Returned to previous state\n');
                drawnow;
            else
                fprintf('Cannot undo, already at initial state\n');
            end
            continue;
        elseif choice == 10
            % Save current components to CSV file
            saveCurrentComponentsToCSV(data_V, data_JD, currents);
            continue;
        elseif choice >= 1 && choice <= 8
            % Save current parameters to history
            history_index = history_index + 1;
            params_history(history_index, :) = adjusted_params;
            
            % Determine parameter index to adjust
            param_idx = ceil(choice / 2);
            
            % Determine adjustment direction
            if mod(choice, 2) == 1
                direction = 1;
            else
                direction = -1;
            end
            
            % Use adaptive step size for adjustment
            delta = adjusted_params(param_idx) * adapt_factors(param_idx) * direction;
            
            % Save old parameters and error for comparison
            old_params = adjusted_params;
            old_error = mean(errors);
            
            % Update parameter
            adjusted_params(param_idx) = adjusted_params(param_idx) + delta;
            
            % Ensure parameters are within reasonable ranges
            if param_idx == 1 % J0
                adjusted_params(param_idx) = max(1e-12, adjusted_params(param_idx));
            elseif param_idx == 2 % Rs - especially emphasize must be positive
                adjusted_params(param_idx) = max(1, adjusted_params(param_idx));
                if adjusted_params(param_idx) <= 0
                    fprintf('Warning: Rs cannot be negative or zero. Adjusted to positive value.\n');
                    adjusted_params(param_idx) = 1; % Ensure positive value
                end
            elseif param_idx == 3 % Rsh
                adjusted_params(param_idx) = max(1e4, adjusted_params(param_idx));
            elseif param_idx == 4 % k
                adjusted_params(param_idx) = max(1e-10, adjusted_params(param_idx));
            end
            
            % Recalculate fit and errors
            fit_results.JD = diodeModel(data_V, adjusted_params, config);
            errors = abs((fit_results.JD - data_JD) ./ (abs(data_JD) + eps)) * 100;
            currents = calculateCurrents(data_V, adjusted_params, config);
            
            % Adaptively adjust step factors
            new_error = mean(errors);
            if new_error < old_error
                % If error decreases, slightly increase step size
                adapt_factors(param_idx) = adapt_factors(param_idx) * 1.2;
                fprintf('Parameter adjustment effective, adaptive step increased\n');
            else
                % If error increases, decrease step size
                adapt_factors(param_idx) = adapt_factors(param_idx) * 0.5;
                fprintf('Parameter adjustment not effective, adaptive step decreased\n');
            end
            
            % Limit adaptive step factors to reasonable range
            adapt_factors(param_idx) = min(0.5, max(0.01, adapt_factors(param_idx)));
            
            % Update chart
            set(h_fit, 'YData', abs(fit_results.JD));
            set(h_error, 'YData', errors);
            title(subplot(2,1,2), sprintf('Fitting Error (Average: %.2f%%)', mean(errors)));
            
            % Update parameter display
            delete(findall(gcf, 'Type', 'annotation'));
            annotation('textbox', [0.01, 0.01, 0.98, 0.08], ...
                'String', sprintf('J0: %.2e A   Rs: %.2e Ohm   Rsh: %.2e Ohm   k: %.2e', ...
                adjusted_params(1), adjusted_params(2), adjusted_params(3), adjusted_params(4)), ...
                'EdgeColor', 'none', 'FontSize', 10, 'HorizontalAlignment', 'center');
            
            drawnow;
        else
            fprintf('Invalid selection, please enter a number between 0-10\n');
        end
    end
    
    % Close the figure
    if ishandle(fig_handle)
        close(fig_handle);
    end
    
    % If user chose to exit without saving, reset parameters to original values
    if ~save_adjustments
        adjusted_params = original_params;
        % Set empty fit results to indicate no adjustments were saved
        fit_results = struct();
        return;
    end
    
    % Final check to ensure Rs is positive
    if adjusted_params(2) <= 0
        fprintf('Warning: Final Rs is negative or zero, adjusted to positive value\n');
        adjusted_params(2) = 1; % Set to a reasonable small positive value
    end
    
    % Calculate final fitting results
    fit_results.JD = diodeModel(data_V, adjusted_params, config);
    fit_results.resnorm = sum(((fit_results.JD - data_JD) ./ (abs(data_JD) + eps)).^2);
    fit_results.currents = calculateCurrents(data_V, adjusted_params, config);
end

function saveCurrentComponentsToCSV(V, measured_JD, currents)
    % Save current components to CSV file
    [filename, pathname] = uiputfile('*.csv', 'Save Current Components to CSV File');
    
    if isequal(filename, 0) || isequal(pathname, 0)
        fprintf('Save cancelled\n');
        return;
    end
    
    fullpath = fullfile(pathname, filename);
    
    % Create table with measured data included
    T = table(V, measured_JD, currents.total, currents.diode, currents.ohmic, currents.nonohmic, ...
        'VariableNames', {'Voltage', 'Measured_Current', 'Total_Current', 'Diode_Component', 'Ohmic_Component', 'NonOhmic_Component'});
    
    % Write to CSV file
    writetable(T, fullpath);
    
    fprintf('Current components successfully saved to: %s\n', fullpath);
end

function displayOptimizationRecommendations(params, errors, data_V, data_JD)
    % Display optimization recommendations based on current fitting state
    fprintf('\nRecommendation based on current fitting:\n');
    
    % Region-specific error analysis
    neg_idx = find(data_V < 0);
    low_pos_idx = find(data_V >= 0 & data_V <= 0.15);
    high_pos_idx = find(data_V > 0.15);
    
    if ~isempty(neg_idx)
        neg_errors = errors(neg_idx);
        fprintf('- Negative voltage region error: %.2f%%\n', mean(neg_errors));
    end
    
    if ~isempty(low_pos_idx)
        low_pos_errors = errors(low_pos_idx);
        fprintf('- Low positive voltage (0-0.15V) error: %.2f%%\n', mean(low_pos_errors));
    end
    
    if ~isempty(high_pos_idx)
        high_pos_errors = errors(high_pos_idx);
        fprintf('- High positive voltage (>0.15V) error: %.2f%%\n', mean(high_pos_errors));
    end
    
    % Parameter-specific recommendations
    if ~isempty(low_pos_idx) && mean(low_pos_errors) > 10
        fprintf('  * Consider adjusting J0 to improve low positive voltage fit\n');
    end
    
    if ~isempty(high_pos_idx) && mean(high_pos_errors) > 10
        fprintf('  * Consider adjusting Rs to improve high positive voltage fit\n');
    end
    
    if ~isempty(neg_idx) && mean(neg_errors) > 10
        fprintf('  * Consider adjusting Rsh and k to improve negative voltage fit\n');
    end
    
    % Parameter value recommendations
    if params(1) < 1e-10 || params(1) > 1e-6
        if params(1) < 1e-10
            range_str = 'below';
        else
            range_str = 'above';
        end
        fprintf('  * J0 (%.2e) is %s typical range, may need adjustment\n', params(1), range_str);
    end
    
    if params(2) < 10 || params(2) > 1e4
        if params(2) < 10
            range_str = 'below';
        else
            range_str = 'above';
        end
        fprintf('  * Rs (%.2e) is %s typical range, may need adjustment\n', params(2), range_str);
    end
    
    if params(3) < 1e5 || params(3) > 1e8
        if params(3) < 1e5
            range_str = 'below';
        else
            range_str = 'above';
        end
        fprintf('  * Rsh (%.2e) is %s typical range, may need adjustment\n', params(3), range_str);
    end
    
    % Overall accuracy assessment
    if mean(errors) < 5
        fprintf('- Excellent fit achieved (error < 5%%)\n');
    elseif mean(errors) < 10
        fprintf('- Good fit achieved (error < 10%%)\n');
    elseif mean(errors) < 20
        fprintf('- Moderate fit achieved (error < 20%%)\n');
    else
        fprintf('- Poor fit, significant improvement needed\n');
    end
end