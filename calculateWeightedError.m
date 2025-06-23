function err = calculateWeightedError(data_V, data_JD, predicted, masks, region_type)
    % CALCULATEWEIGHTEDERROR 计算加权误差
    %   ERR = CALCULATEWEIGHTEDERROR(DATA_V, DATA_JD, PREDICTED, MASKS, REGION_TYPE)
    %   根据不同区域的权重计算测量数据和预测数据之间的误差。
    %
    %   输入:
    %     DATA_V - 电压数据向量
    %     DATA_JD - 测量的电流密度数据向量
    %     PREDICTED - 预测的电流密度数据向量
    %     MASKS - 由createVoltageMasks函数创建的掩码结构体
    %     REGION_TYPE - 区域类型，可选值：'all', 'negative', 'positive', 'enhanced_positive'
    %
    %   输出:
    %     ERR - 加权误差向量
    
    if nargin < 5
        region_type = 'all';
    end
    
    % 预计算
    max_data_abs = max(abs(data_JD));
    threshold = 1e-12;
    
    % 初始化误差向量
    err = zeros(size(data_JD));
    
    % 处理低于阈值的值
    if isfield(masks, 'below_threshold')
        err(masks.below_threshold) = (predicted(masks.below_threshold) - data_JD(masks.below_threshold)) / max(threshold, max_data_abs);
    end
    
    % 计算对数空间误差
    if isfield(masks, 'above_threshold') && any(masks.above_threshold)
        actual_abs = abs(data_JD(masks.above_threshold));
        pred_abs = abs(predicted(masks.above_threshold));
        log_actual = log10(actual_abs);
        log_pred = log10(pred_abs);
        err(masks.above_threshold) = log_pred - log_actual;
    end
    
    % 应用符号不匹配权重
    if isfield(masks, 'sign_mismatch')
        err(masks.sign_mismatch) = err(masks.sign_mismatch) * 4;
    end
    
    % 根据区域类型应用不同的权重策略
    switch region_type
        case 'negative'
            % 负电压区域优化
            if isfield(masks, 'strong_neg')
                err(masks.strong_neg) = err(masks.strong_neg) * 8;
            end
            if isfield(masks, 'weak_neg')
                err(masks.weak_neg) = err(masks.weak_neg) * 5;
            end
            if isfield(masks, 'most_neg')
                err(masks.most_neg) = err(masks.most_neg) * 3;
            end
            
        case 'positive'
            % 正电压区域优化
            pos_mask = data_V >= 0;
            err(pos_mask) = err(pos_mask) * 3;
            
            % 细分正电压区域权重
            low_pos_mask = data_V >= 0 & data_V < 0.15;
            mid_pos_mask = data_V >= 0.15 & data_V < 0.20;
            high_pos_mask = data_V >= 0.20;
            
            err(low_pos_mask) = err(low_pos_mask) * 2;
            err(mid_pos_mask) = err(mid_pos_mask) * 3;
            err(high_pos_mask) = err(high_pos_mask) * 4;
            
        case 'enhanced_positive'
            % 增强的正电压区域优化
            low_pos_mask = data_V >= 0 & data_V < 0.15;
            mid_pos_mask = data_V >= 0.15 & data_V < 0.20;
            high_pos_mask = data_V >= 0.20;
            
            err(low_pos_mask) = err(low_pos_mask) * 8;
            err(mid_pos_mask) = err(mid_pos_mask) * 25;
            err(high_pos_mask) = err(high_pos_mask) * 50;
            
            % 降低负电压区域权重
            neg_mask = data_V < 0;
            err(neg_mask) = err(neg_mask) * 0.05;
            
            % 最高正偏压区域额外加权
            very_high_pos_mask = data_V >= 0.25;
            if any(very_high_pos_mask)
                err(very_high_pos_mask) = err(very_high_pos_mask) * 1.5;
            end
            
        otherwise
            % 默认权重
            if isfield(masks, 'strong_neg')
                err(masks.strong_neg) = err(masks.strong_neg) * 5.0;
            end
            if isfield(masks, 'weak_neg')
                err(masks.weak_neg) = err(masks.weak_neg) * 3.0;
            end
            if isfield(masks, 'most_neg')
                err(masks.most_neg) = err(masks.most_neg) * 2.0;
            end
            if isfield(masks, 'near_zero')
                err(masks.near_zero) = err(masks.near_zero) * 2.0;
            end
    end
end
