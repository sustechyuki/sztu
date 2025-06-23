function [data_V, data_JD] = loadData()
    % LOADDATA Load I-V measurement data from Excel or CSV file
    %   [data_V, data_JD] = LOADDATA() prompts user to select a data file and
    %   returns voltage and current density data for diode characterization.
    %
    %   Outputs:
    %       data_V  - Voltage data vector (V)
    %       data_JD - Current density data vector (A/cm²)
    %
    %   File Requirements:
    %       - Supported formats: Excel (.xlsx, .xls) or CSV
    %       - Data should be arranged in columns
    %       - User will be prompted to specify the cell range
    %
    %   Note: Currently uses a predefined voltage range (-0.5V to 0.3V)
    %   for data interpolation.

    % Open file dialog to select file
    [filename, pathname] = uigetfile({'*.xlsx;*.xls', 'Excel Files (*.xlsx, *.xls)'; ...
                                     '*.csv', 'CSV Files (*.csv)'; ...
                                     '*.*', 'All Files (*.*)'}, ...
                                     'Select file containing IV data');
    % Build full file path
    fullpath = fullfile(pathname, filename);
    fprintf('Reading file: %s\n', fullpath);
        
    % Prompt user for cell range
    range_str = input('Enter cell range to read (e.g.: A2:B82): ', 's');
    
    % Read specified data range
    data = readmatrix(fullpath, 'Range', range_str);
    fprintf('Data successfully read, size: %d rows x %d columns\n', size(data, 1), size(data, 2));
    
    % Voltage data (predefined range)
    data_V = -0.5:0.01:0.3;
    data_JD = data;
    % data_JD = []
    
    % Convert to column vectors
    data_V = data_V(:);
    data_JD = data_JD(:);
end