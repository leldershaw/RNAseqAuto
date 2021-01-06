% LE 6th January 2021

%% Imports and calculates statistics for gene and exon level RNAseq data
function import_and_calculate(analysis_type, files, files_path, rows_file, rows_path, output_type)
    warning('off','MATLAB:table:ModifiedAndSavedVarnames'); %Suppress table creation warning
    
    % Import rows file
    if ~isequal(rows_file, 0)
        rows = readcell(fullfile(rows_path, rows_file));
    end
    
    pair_data = cell(size(files,2)/2, 3);
    avg_exon = cell(size(files,2)/2, 2);
    
    % Separate summary and count files
    summ_files = files(contains(files, "Summary"));
    count_files = files(~contains(files, "Summary"));
    
    for i = 1:size(summ_files, 2)
        % Import summary file
        T = readtable(fullfile(files_path, summ_files{i}), 'FileType', 'text', 'ReadVariableNames', true);
       
        pair_data{i,1} = T.Properties.VariableDescriptions(2); % data ID
        avg_exon{i,1} = T.Properties.VariableDescriptions(2); % data ID
        pair_data{i,2} = T{1,2}; % assigned/mapped reads
        
        % Find count file paired with summary file
        s = split(summ_files{i}, ["[","]", "__"]);
        if length(find(contains(count_files, s{2}))) ~= 1
            error("MATLAB:find_exons_and_calculate:invalidPair",...
                                "ERROR: '%s' does not have exactly 1 count file.", s{2});
        else
            count = count_files{contains(count_files, s{2})};
        end
        
        % Import count file
        T = readtable(fullfile(files_path, count), 'FileType', 'text', 'ReadVariableNames', true);
        
        % Perform Gene and Exon analysis
        if isequal(analysis_type, 'Gene_select') % GENE - select genes/rows option
            try
                pair_data{i,3} = calculate_genes(T, rows, pair_data{i,2});
            catch ME % catch errors in row file
                if (contains(ME.identifier,'MATLAB:badsubscript')) % eg. row 1 is the header so invalid
                    error('MATLAB:import_and_calculate:invalidRow',...
                        'ERROR: A row listed in the Rows Text File is not valid.');
                else % rethrow error if something else
                    rethrow(ME);
                end
            end
        elseif isequal(analysis_type, 'Gene_all') % GENE - all data option
            pair_data{i,3} = calculate_genes(T, 0, pair_data{i,2});
        elseif isequal(analysis_type, 'Exon') % EXON
            [pair_data{i,3}, avg_exon{i,2}] = calculate_exons(T, pair_data{i,2});
        end 
    end
    
    warning('on','MATLAB:table:ModifiedAndSavedVarnames');

    %% Create output folder
    folder_name = strcat(datestr(now,'yyyy-mm-dd HH.MM.SS PM'), '  RNASeqAuto');
    mkdir(files_path, folder_name);
    out_path = fullfile(files_path, folder_name);
    
    %% Save execution log
    c_files = strings(size(files, 2),1);
    for n = 1:size(files, 2)
        c_files(n) = files{n};
    end
    
    out_file = [];
    if output_type(1)
        out_file = [out_file; "Single Sheet"];
    end
    if output_type(2)
        out_file = [out_file; "Sorted by GeneIDs"];
    end
    if isequal(out_file, [])
        out_file = "None";
    end
    
    data_type = 'All data';
    if contains(analysis_type, 'Gene')
        level_type = 'Gene';
        if contains(analysis_type, 'select')
            data_type = 'Select genes/rows';
        end
    else
        level_type = 'Exon';
    end
    
    executionFile = [["Analysis Type:" level_type];
                     ["Data Analysed:" data_type];
                     [["Output File Type:"; strings(size(out_file, 1)-1)] out_file];
                     ["Output Folder:" out_path];
                     ["Input Folder:" files_path];
                     [["Counts + Matched Summary Files:"; strings(size(c_files, 1)-1,1)] sort(c_files)]];
           
    if isequal(analysis_type, 'Gene_select')
        executionFile = [executionFile; ["Rows Text File" fullfile(rows_path,rows_file)]];
    end

    writematrix(executionFile, fullfile(out_path, 'execution_log.csv'));
    
    %% Save average exon calculations to excel file
    if isequal(analysis_type, 'Exon')
        dataFile3 = [{'Counts File ID'}, {'Geneid'}, {'Length'}, {'Counts'}, {'Average RPKM'}];
        for m = 1:size(avg_exon,1)
            dataFile3 = [dataFile3; [repmat(avg_exon{m,1},size(avg_exon{m,2}, 1)-1,1), avg_exon{m,2}(2:end,:)]];
        end
        writecell(dataFile3, fullfile(out_path, 'output_averageExon.xlsx'));
    end

    %% Save to excel file - single sheet
    if output_type(1)
        dataFile = ['Counts File ID', pair_data{1,3}(1,:)];
        for m = 1:size(pair_data,1)
            dataFile = [dataFile; [repmat(pair_data{m,1},size(pair_data{m,3}, 1)-1,1), pair_data{m,3}(2:end,:)]];
        end

        writecell(dataFile, fullfile(out_path, 'output_singleSheet.xlsx'));
    end
    
    %% Save to excel file - sorted by GeneIDs into sheets
    if output_type(2)
        for m = 1:size(pair_data,1)
            for n = 2:size(pair_data{m,3},1)
                dataFile2 = [pair_data{m,1}, pair_data{m,3}(n,2:end)];
                writecell(['Counts File ID', pair_data{1,3}(1,2:end)], ...
                    fullfile(out_path, 'output_sortedByGeneIDs.xlsx'), ...
                    'Sheet', pair_data{m,3}{n,1}, 'Range', 'A1');
                writecell(dataFile2, fullfile(out_path, 'output_sortedByGeneIDs.xlsx'), ...
                    'Sheet', pair_data{m,3}{n,1}, 'WriteMode', 'append');
            end
        end
    end
    
end

    