function err = errorFunction(x, data_V, data_JD, params, config, options)
    % ERRORFUNCTION Universal error function for diode model fitting
    %
    % This function calculates weighted error between measured data and model predictions.
    % It supports different error calculation modes and region-specific optimizations.
    %
    % Inputs:
    %   x - Scaled parameter vector [J0, Rs, Rsh, k]
    %   data_V - Measured voltage data
    %   data_JD - Measured current density data
    %   params - Parameter structure containing scaling factors
    %   config - Configuration structure
    %   options - Optional structure with fields:
    %       .mode - Error calculation mode ('standard', 'partial', 'enhanced_positive', etc.)
    %       .region_weights - Custom region weights (struct)
    %       .param_mask - Boolean mask for parameters to optimize
    %       .x_full - Full parameter vector (for partial optimization)
    %
    % Output:
    %   err - Weighted error vector
    
    % Default options
    if nargin < 6
        options = struct();
    end
    
    % Set default mode
    if ~isfield(options, 'mode')
        options.mode = 'standard';
    end
    
    % Handle partial parameter optimization
    if isfield(options, 'param_mask') && isfield(options, 'x_full')
        x_full = options.x_full;
        x_full(options.param_mask) = x;
        x = x_full;
    end
    
    % Rescale parameters to actual values
    x_actual = x .* params.scaleFactors;
    
    % Validate physical parameters
    x_actual = validateParameters(x_actual);
    
    % Calculate model predictions
    predicted = diodeModel(data_V, x_actual, config);
    
    % Create voltage region masks
    masks = createVoltageMasks(data_V, data_JD, predicted);
    
    % Apply mode-specific processing
    switch options.mode
        case 'partial'
            % For optimization of specific regions
            if isfield(options, 'region_mask')
                selected_data = options.region_mask;
                err = zeros(size(data_JD));
                err(selected_data) = calculateWeightedError(data_JD(selected_data), predicted(selected_data), filterMasks(masks, selected_data));
                err(~selected_data) = 0;
            else
                err = calculateWeightedError(data_JD, predicted, masks);
            end
            
        case 'enhanced_positive'
            % Enhanced optimization for positive voltage region
            pos_weights = struct('strong_neg', 0.05, 'weak_neg', 0.05, 'low_pos', 8, 'mid_pos', 25, 'high_pos', 50, 'very_high_pos', 75);
            err = calculateWeightedError(data_JD, predicted, masks, pos_weights);
            
        case 'strong_positive'
            % Special optimization for high positive voltage region
            high_pos_weights = struct('strong_neg', 0.01, 'weak_neg', 0.01, 'low_pos', 0.1, 'mid_pos', 0.5, 'high_pos', 10, 'very_high_pos', 20);
            err = calculateWeightedError(data_JD, predicted, masks, high_pos_weights);
            
        case 'negative'
            % Focused on negative voltage region
            neg_weights = struct('strong_neg', 10, 'weak_neg', 5, 'mid_neg', 8, 'low_pos', 0.1, 'mid_pos', 0.1, 'high_pos', 0.1);
            err = calculateWeightedError(data_JD, predicted, masks, neg_weights);
            
        case 'balanced'
            % Balanced weights for all regions
            balanced_weights = struct('strong_neg', 5, 'weak_neg', 5, 'low_pos', 5, 'mid_pos', 5, 'high_pos', 5);
            err = calculateWeightedError(data_JD, predicted, masks, balanced_weights);
            
        case 'custom'
            % Custom weights provided by caller
            if isfield(options, 'region_weights')
                err = calculateWeightedError(data_JD, predicted, masks, options.region_weights);
            else
                err = calculateWeightedError(data_JD, predicted, masks);
            end
            
        otherwise
            % Standard error calculation
            err = calculateWeightedError(data_JD, predicted, masks);
    end
end

function filtered_masks = filterMasks(masks, region_mask)
    % Filter existing masks to apply only to selected region
    filtered_masks = masks;
    fields = fieldnames(masks);
    for i = 1:length(fields)
        field = fields{i};
        filtered_masks.(field) = masks.(field) & region_mask;
    end
end

function err = calculateWeightedError(data_JD, predicted, masks, custom_weights)
    % Calculate weighted error with customizable weights
    
    % Initialize error vector
    err = zeros(size(data_JD));
    threshold = 1e-12;
    
    % Calculate logarithmic error for values above threshold
    valid_idx = (abs(data_JD) > threshold) & (abs(predicted) > threshold);
    if any(valid_idx)
        % Logarithmic error component
        log_actual = log10(abs(data_JD(valid_idx)));
        log_pred = log10(abs(predicted(valid_idx)));
        err(valid_idx) = log_pred - log_actual;
        
        % Add relative error component for large values
        large_idx = valid_idx & (abs(data_JD) > 1e-9);
        if any(large_idx)
            rel_err = (predicted(large_idx) - data_JD(large_idx)) ./ abs(data_JD(large_idx));
            err(large_idx) = err(large_idx) + 0.2 * rel_err; % Add 20% relative error component
        end
    end
    
    % Handle values below threshold
    below_idx = ~valid_idx;
    if any(below_idx)
        err(below_idx) = (predicted(below_idx) - data_JD(below_idx)) / threshold;
    end
    
    % Apply sign mismatch penalty
    sign_mismatch = sign(predicted) ~= sign(data_JD) & abs(data_JD) > threshold;
    err(sign_mismatch) = err(sign_mismatch) * 4;
    
    % Apply region-specific weights
    % Default weights
    weights = struct(...
        'strong_neg', 5.0, ...
        'mid_neg', 4.0, ...
        'weak_neg', 3.0, ...
        'near_zero', 2.0, ...
        'low_pos', 3.0, ...
        'mid_pos', 4.0, ...
        'high_pos', 5.0, ...
        'very_high_pos', 6.0);
    
    % Override with custom weights if provided
    if nargin > 3 && isstruct(custom_weights)
        fields = fieldnames(custom_weights);
        for i = 1:length(fields)
            field = fields{i};
            if isfield(weights, field)
                weights.(field) = custom_weights.(field);
            end
        end
    end
    
    % Apply weights to appropriate regions
    fields = fieldnames(masks);
    for i = 1:length(fields)
        field = fields{i};
        if isfield(weights, field) && any(masks.(field))
            err(masks.(field)) = err(masks.(field)) * weights.(field);
        end
    end
end

function masks = createVoltageMasks(data_V, data_JD, predicted)
    % Create masks for different voltage regions
    
    % Calculate thresholds and absolute values
    threshold = 1e-12;
    actual_abs = abs(data_JD);
    pred_abs = abs(predicted);
    
    % Create basic masks
    masks = struct();
    masks.below_threshold = (actual_abs < threshold) | (pred_abs < threshold);
    masks.above_threshold = ~masks.below_threshold;
    masks.sign_mismatch = sign(predicted) ~= sign(data_JD) & actual_abs > threshold & pred_abs > threshold;
    
    % Negative voltage region masks
    masks.strong_neg = data_V < -0.3;
    masks.mid_neg = data_V >= -0.3 & data_V < -0.15;
    masks.weak_neg = data_V >= -0.15 & data_V < 0;
    masks.near_zero = abs(data_V) < 0.05;
    
    % Positive voltage region masks
    masks.low_pos = data_V >= 0 & data_V < 0.15;
    masks.mid_pos = data_V >= 0.15 & data_V < 0.20;
    masks.high_pos = data_V >= 0.20 & data_V < 0.25;
    masks.very_high_pos = data_V >= 0.25;
    
    % Additional region segments can be defined here as needed
end