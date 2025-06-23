function currents = calculateCurrents(V, x, config, varargin)
    % CALCULATECURRENTS Calculate current components for diode model
    %   CALCULATECURRENTS(V, x, config) calculates the total current and its
    %   components (diode, ohmic, non-ohmic) based on the diode model parameters.
    %
    %   Inputs:
    %     V      - Voltage data vector (V)
    %     x      - Parameter vector:
    %              x(1): J0 - Saturation current (A)
    %              x(2): Rs - Series resistance (Ω)
    %              x(3): Rsh - Shunt resistance (Ω)
    %              x(4): k - Non-linearity coefficient
    %     config - Configuration structure with physics parameters
    %
    %   Optional Input:
    %     varargin{1} - Boolean flag to export results to Excel
    %
    %   Output:
    %     currents - Structure containing current components:
    %                .total    - Total fitted current
    %                .diode    - Diode current component
    %                .ohmic    - Ohmic current component
    %                .nonohmic - Non-ohmic current component
    %                .sum      - Sum of components (for verification)
    
    % Calculate total current from model
    currents.total = diodeModel(V, x, config);
    
    % Calculate voltage drop
    V_drop = V - currents.total .* x(2);
    
    % Calculate individual current components 
    currents.diode = x(1) * (exp(config.physics.A * V_drop / config.physics.n) - 1);
    currents.ohmic = V_drop / x(3);
    currents.nonohmic = x(4) * (abs(V_drop).^config.physics.m) .* sign(V_drop);
    
    % Calculate sum of components (for verification)
    currents.sum = currents.diode + currents.ohmic + currents.nonohmic;
    
    % Check if Excel output is requested
    if nargin > 3 && ~isempty(varargin{1}) && varargin{1}
        % Create output data table
        outputTable = table(V, V_drop, currents.total, currents.diode, currents.ohmic, currents.nonohmic, currents.sum, ...
            'VariableNames', {'Voltage_V', 'Voltage_Drop', 'Total_Current', 'Diode_Current', 'Ohmic_Current', 'NonOhmic_Current', 'Sum_Components'});
        
        % Generate filename with timestamp to avoid overwriting
        timestamp = datestr(now, 'yyyymmdd_HHMMSS');
        filename = fullfile(pwd, ['Fitting_Currents_', timestamp, '.xlsx']);
        
        % Write to Excel file
        writetable(outputTable, filename);
        disp(['Current components saved to: ', filename]);
    end
end