% LE 6th January 2021

%% Calculates statistics for gene level data
function data = calculate_genes(T, rows, assigned_reads)
    
    % Gather data
    if ~isequal(rows, 0)
        raw_data = table2cell(T([rows{:}]-1, :)); % All data option
    else
        raw_data = sortrows(table2cell(T), 1); % Select genes/rows option
    end
    
    % Initialise headers
    data = [[T.Properties.VariableDescriptions, {'Length (kb)'}, {'RPK'}, {'Assigned reads'}, {'Assigned reads (millions)'}, {'RPKM'}];...
            [raw_data, cell(size(raw_data, 1),5)]];
    data{1,7} = 'Counts';
    
    % Perform calculations
    for i = 2:size(data, 1)
        data{i, 8} = data{i, 6}/1000;
        data{i, 9} = data{i, 7}/data{i,8};
        data{i, 10} = assigned_reads;
        data{i, 11} = assigned_reads/1000000;
        data{i, 12} = data{i, 9}/data{i, 11};
    end
end
