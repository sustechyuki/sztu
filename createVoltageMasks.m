function masks = createVoltageMasks(data_V, data_JD, predicted)
    % CREATEVOLTAGEMASKS 创建不同电压区域的掩码
    %   MASKS = CREATEVOLTAGEMASKS(DATA_V, DATA_JD, PREDICTED) 根据电压和电流数据
    %   创建不同区域的掩码，用于误差计算和权重应用。
    %
    %   输入:
    %     DATA_V - 电压数据向量
    %     DATA_JD - 测量的电流密度数据向量
    %     PREDICTED - 预测的电流密度数据向量
    %
    %   输出:
    %     MASKS - 包含不同区域掩码的结构体
    
    % 预计算绝对值和阈值
    threshold = 1e-12;
    actual_abs = abs(data_JD);
    pred_abs = abs(predicted);
    
    % 创建基本掩码
    masks = struct();
    masks.below_threshold = (actual_abs < threshold) | (pred_abs < threshold);
    masks.above_threshold = ~masks.below_threshold;
    masks.sign_mismatch = sign(predicted) ~= sign(data_JD) & actual_abs > threshold & pred_abs > threshold;
    
    % 负电压区域掩码
    masks.strong_neg = data_V < -0.2;
    masks.weak_neg = data_V >= -0.2 & data_V < 0;
    masks.near_zero = abs(data_V) < 0.05;
    
    % 正电压区域掩码
    masks.low_pos = data_V >= 0 & data_V < 0.15;
    masks.mid_pos = data_V >= 0.15 & data_V < 0.20;
    masks.high_pos = data_V >= 0.20;
    masks.very_high_pos = data_V >= 0.25;
    
    % 分段处理负电压区域
    neg_voltages = unique(data_V(data_V < 0));
    n_neg_segments = 4; % 将负电压范围分为4段
    
    if ~isempty(neg_voltages)
        neg_segment_size = ceil(length(neg_voltages) / n_neg_segments);
        
        % 识别最负电压段（第一段）
        if length(neg_voltages) >= neg_segment_size
            most_neg_voltages = neg_voltages(1:min(neg_segment_size, length(neg_voltages)));
            
            % 创建最负电压点的掩码
            most_neg_mask = false(size(data_V));
            for v = most_neg_voltages'
                most_neg_mask = most_neg_mask | (data_V == v);
            end
            
            masks.most_neg = most_neg_mask;
        end
        
        % 创建其他负电压段的掩码
        for i = 2:n_neg_segments
            start_idx = (i-1) * neg_segment_size + 1;
            end_idx = min(i * neg_segment_size, length(neg_voltages));
            
            if start_idx <= length(neg_voltages)
                segment_voltages = neg_voltages(start_idx:end_idx);
                segment_mask = false(size(data_V));
                
                for v = segment_voltages'
                    segment_mask = segment_mask | (data_V == v);
                end
                
                field_name = sprintf('neg_segment_%d', i);
                masks.(field_name) = segment_mask;
            end
        end
    end
end