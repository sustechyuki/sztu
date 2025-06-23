function validateInputData(V, JD)
    % VALIDATEINPUTDATA Validate voltage and current density data
    %   VALIDATEINPUTDATA(V, JD) validates that the input voltage and
    %   current density data meet basic requirements for diode model fitting.
    %
    %   Validation checks:
    %   - Data is not empty
    %   - Vectors have matching lengths
    %   - No NaN values present
    %   - No Inf values present
    %   - Data is numeric type
    
    % Check if input data is empty
    if isempty(V) || isempty(JD)
        error('Input data cannot be empty');
    end
    
    % Check if data lengths match
    if length(V) ~= length(JD)
        error('Voltage and current data must have the same length');
    end
    
    % Check for NaN values
    if any(isnan(V)) || any(isnan(JD))
        error('Input data contains NaN values');
    end
    
    % Check for Inf values
    if any(isinf(V)) || any(isinf(JD))
        error('Input data contains infinite values');
    end
    
    % Check data type
    if ~isnumeric(V) || ~isnumeric(JD)
        error('Input data must be numeric type');
    end
end