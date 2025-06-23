function params = validateParameters(params, param_names)
    % VALIDATEPARAMETERS 验证参数的物理合理性
    %   PARAMS = VALIDATEPARAMETERS(PARAMS, PARAM_NAMES) 检查参数是否在物理上合理，
    %   并在必要时进行调整。
    %
    %   输入:
    %     PARAMS - 参数向量 [J0, Rs, Rsh, k]
    %     PARAM_NAMES - 参数名称的元胞数组（可选）
    %
    %   输出:
    %     PARAMS - 调整后的参数向量
    
    if nargin < 2
        param_names = {'J0', 'Rs', 'Rsh', 'k'};
    end
    
    % 检查Rs是否为正值
    if params(2) <= 0
        fprintf('警告: %s 为负值或零，自动调整为正值\n', param_names{2});
        params(2) = 10; % 使用合理的默认值
    end
    
    % 检查Rsh是否为正值
    if params(3) <= 0
        fprintf('警告: %s 为负值或零，自动调整为正值\n', param_names{3});
        params(3) = 1e6; % 使用合理的默认值
    end
    
    % 检查J0是否为正值
    if params(1) <= 0
        fprintf('警告: %s 为负值或零，自动调整为正值\n', param_names{1});
        params(1) = 1e-9; % 使用合理的默认值
    end
    
    % 检查k是否为正值
    %if params(4) <= 0
        %fprintf('警告: %s 为负值或零，自动调整为正值\n', param_names{4});
        %params(4) = 1e-7; % 使用合理的默认值
    %end
end