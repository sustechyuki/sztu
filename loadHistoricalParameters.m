function params = loadHistoricalParameters(data_V, data_JD, config)
    % LOADHISTORICALPARAMETERS 从历史参数文件中加载初始参数
    %   PARAMS = LOADHISTORICALPARAMETERS(DATA_V, DATA_JD, CONFIG) 列出可用的历史参数文件
    %   并允许用户选择一个作为初始参数。如果没有可用的文件或用户取消选择，则使用
    %   Lambert W 函数估计初始参数。
    %
    %   输入:
    %     DATA_V  - 电压数据向量
    %     DATA_JD - 电流密度数据向量
    %     CONFIG  - 配置结构体
    %
    %   输出:
    %     PARAMS  - 参数结构体，包含初始参数和边界

    % List available parameter files
    mat_files = dir('adjusted_params_*.mat');
    txt_files = dir('adjusted_params_*.txt');
    
    if isempty(mat_files)
        fprintf('No historical parameter files found, will use Lambert W function to estimate initial parameters\n');
        params = initializeParameters(data_V, data_JD, config);
    else
        % Display available files
        fprintf('Available parameter files:\n');
        for i = 1:length(mat_files)
            fprintf('%d: %s\n', i, mat_files(i).name);
        end
        
        % Let user select a file
        file_idx = input('Please select a file number (enter 0 to cancel): ');
        if file_idx > 0 && file_idx <= length(mat_files)
            % Load the selected file
            load_file = mat_files(file_idx).name;
            loaded_data = load(load_file);
            
            % Extract parameters
            if isfield(loaded_data, 'params')
                fprintf('Loading parameters from file %s\n', load_file);
                
                % Create parameter structure
                params = struct();
                params.x0 = loaded_data.params;
                
                % Parameter range settings
                params.ub = [1e-6, 1e5, 1e10, 1e-4];    % Upper bounds
                params.lb = [1e-12, 1e1, 1e5, 1e-10];    % Lower bounds
                
                % Ensure initial values are within range
                params.x0 = min(max(params.x0, params.lb), params.ub);
                
                % Scale factors
                params.scaleFactors = [1e-9, 1e3, 1e7, 1e-8];
                
                % Display loaded parameters
                fprintf('Loaded parameters:\n');
                fprintf('J0 = %.6e A\n', params.x0(1));
                fprintf('Rs = %.6e Ohm\n', params.x0(2));
                fprintf('Rsh = %.6e Ohm\n', params.x0(3));
                fprintf('k = %.6e\n', params.x0(4));
            else
                fprintf('File format error, will use Lambert W function to estimate initial parameters\n');
                params = initializeParameters(data_V, data_JD, config);
            end
        else
            fprintf('Cancelled loading parameter file, will use Lambert W function to estimate initial parameters\n');
            params = initializeParameters(data_V, data_JD, config);
        end
    end
end