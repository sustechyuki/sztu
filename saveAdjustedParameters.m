function saveAdjustedParameters(params)
    % SAVEADJUSTEDPARAMETERS 保存调整后的参数
    %   SAVEADJUSTEDPARAMETERS(PARAMS) 将调整后的参数保存到MAT文件和TXT文件中，
    %   文件名包含时间戳以避免覆盖。
    %
    %   输入:
    %     PARAMS - 参数向量 [J0, Rs, Rsh, k]

    % Generate timestamp
    timestamp = datestr(now, 'yyyymmdd_HHMMSS');
    filename = sprintf('adjusted_params_%s.mat', timestamp);
    
    % Save to file
    save(filename, 'params');
    fprintf('Parameters saved to file: %s\n', filename);
    
    % Also export as text file for easy viewing
    txt_filename = sprintf('adjusted_params_%s.txt', timestamp);
    fid = fopen(txt_filename, 'w');
    fprintf(fid, 'J0 = %.6e A\n', params(1));
    fprintf(fid, 'Rs = %.6e Ohm\n', params(2));
    fprintf(fid, 'Rsh = %.6e Ohm\n', params(3));
    fprintf(fid, 'k = %.6e\n', params(4));
    fclose(fid);
    fprintf('Parameters exported as text file: %s\n', txt_filename);
end