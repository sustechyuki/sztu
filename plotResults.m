function plotResults(V, JD_measured, fit_results, currents)
    % PLOTRESULTS Generate comprehensive visualization of diode fitting results
    %   PLOTRESULTS(V, JD_measured, fit_results, currents) creates a multi-panel
    %   figure displaying the diode I-V characteristics with measured data and
    %   fitted model components.
    %
    %   Inputs:
    %     V           - Voltage data vector (V)
    %     JD_measured - Measured current density data (A)
    %     fit_results - Structure containing fitting parameters
    %     currents    - Structure with current components:
    %                   .total    - Total fitted current
    %                   .diode    - Diode current component
    %                   .ohmic    - Ohmic current component
    %                   .nonohmic - Non-ohmic current component
    
    % Create figure window for detailed analysis
    figure('Position', [100, 100, 1200, 800]);
    % set(gcf, 'DefaultAxesFontName', 'Times New Roman');

    % First subplot: I-V characteristics (log scale) with current components
    subplot(2,2,[1,2]);

    % Plot measured data and fitting results
    semilogy(V, abs(JD_measured), 'bo', 'DisplayName', 'measure data', 'MarkerSize', 6);
    hold on;
    semilogy(V, abs(currents.total), 'ro', 'DisplayName', 'total fit', 'MarkerSize', 6);
    semilogy(V, abs(currents.diode), 'b--', 'DisplayName', 'diodecurrent', 'LineWidth', 1.5);
    semilogy(V, abs(currents.ohmic), 'g--', 'DisplayName', 'Ohmic loss', 'LineWidth', 1.5);
    semilogy(V, abs(currents.nonohmic), 'm--', 'DisplayName', 'non-Ohmic loss', 'LineWidth', 1.5);
    
    % % Alternative plotting style with consistent naming convention
    % semilogy(V, abs(JD_measured), 'bo', 'DisplayName', 'Measured Data', 'MarkerSize', 6, 'MarkerFaceColor', 'w');
    % hold on;
    % semilogy(V, abs(currents.total), 'r^', 'DisplayName', 'Total Fit', 'MarkerSize', 5, 'MarkerFaceColor', [1 0.5 0]);
    % semilogy(V, abs(currents.diode), 'b-', 'DisplayName', 'Diode Current', 'LineWidth', 1.8);
    % semilogy(V, abs(currents.ohmic), 'g:', 'DisplayName', 'Ohmic Loss', 'LineWidth', 2.2);
    % semilogy(V, abs(currents.nonohmic), 'm-.', 'DisplayName', 'Recombination', 'LineWidth', 1.8);    
    
    % Set axis labels and properties
    xlabel('Voltage (V)', 'FontSize', 12);
    ylabel('Current Density (A)', 'FontSize', 12);
    title('Diode I-V Characteristics and Current Components (semi Log)', 'FontSize', 14);
    grid on;
    legend('Location', 'best');
    
    % Set appropriate axis limits
    xlim([min(V) max(V)]);
    ylim([min(abs(JD_measured))*0.1 max(abs(JD_measured))*10]);
    
    % Second subplot: Relative error (linear scale)
    subplot(2,2,3);
    relative_error = abs((currents.total - JD_measured) ./ (abs(JD_measured) + eps)) * 100;
    plot(V, relative_error, 'b.-', 'LineWidth', 1);
    xlabel('Voltage (V)');
    ylabel('Relative Error (%)');
    title('Fitting Relative Error');
    grid on;
    
    % Third subplot: I-V characteristics (linear scale)
    subplot(2,2,4);
    plot(V, JD_measured, 'bo', 'DisplayName', 'Measured Data');
    hold on;
    plot(V, currents.total, 'ro', 'DisplayName', 'Fitting Result');
    xlabel('Voltage (V)');
    ylabel('Current Density (A)');
    title('I-V Characteristics (Linear Scale)');
    grid on;
    legend('Location', 'best');
    
    % Add overall title
    sgtitle('Diode Characteristics Fitting Analysis', 'FontSize', 14);
    
    % Add fitting parameter information on the figure
    annotation('textbox', [0.02, 0.02, 0.4, 0.05], ...
        'String', sprintf('Max Relative Error: %.2f%%  Mean Relative Error: %.2f%%', ...
        max(relative_error), mean(relative_error)), ...
        'EdgeColor', 'none');
    
    % Print current component percentages
    fprintf('Diode current percentage: %.2f%%\n', mean(abs(currents.diode ./ currents.total)) * 100);
    fprintf('Ohmic current percentage: %.2f%%\n', mean(abs(currents.ohmic ./ currents.total)) * 100);
    fprintf('Non-ohmic current percentage: %.2f%%\n', mean(abs(currents.nonohmic ./ currents.total)) * 100);
end