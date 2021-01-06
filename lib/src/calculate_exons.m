% LE 6th January 2020

%% Finds exons and calculates statistics for exon and subsequent gene level data
function [superexons, avg_exon] = calculate_exons(T, assigned_reads)
    
    % Gather sorted data + headers
    data = [T.Properties.VariableDescriptions; sortrows(table2cell(T), [1,2,5,3])];
    
    % Initialise headers
    exlist = cell(size(data, 1), 15);
    exlist(1,:) = [T.Properties.VariableDescriptions(1:6), {'Counts'}, {'Length (kb)'}, ...
        {'RPK'}, {'Assigned reads'}, {'Assigned reads (millions)'}, {'RPKM'}, ...
        {'Gene Length'}, {'Gene Counts'}, {'Average RPKM'}];
    rpkm_sum = 0;
    exon_num = 0;
    
    j = 2;
    for i = 2:size(data,1)

        if i==2 % initial step
            exlist{j, 1} = data{i, 1};
            exlist{j, 2} = data{i, 2};
            exlist{j, 3} = data{i, 3};
            exlist{j, 4} = data{i, 4};
            exlist{j, 5} = data{i, 5};
            exlist{j, 6} = exlist{j, 4} - exlist{j, 3} + 1;
            exlist{j, 7} = data{i, 7};
            gene_length = exlist{j, 6};
            gene_counts = exlist{j, 7};
        % If exon does not overlap with previous exon or has different GeneID, then new gene    
        elseif data{i, 3} > exlist{j, 4} || ~strcmp(data{i, 1}, exlist{j, 1})                           
            exlist{j, 6} = exlist{j, 4} - exlist{j, 3} + 1;
            
            % Exon Level Calculations
            exlist{j, 8} = exlist{j, 6}/1000;
            exlist{j, 9} = exlist{j, 7}/exlist{j, 8};
            exlist{j, 10} = assigned_reads;
            exlist{j, 11} = assigned_reads/1000000;
            exlist{j, 12} = exlist{j, 9}/exlist{j, 11};

            if ~strcmp(exlist{j-1, 1}, exlist{j, 1}) % if different GeneID
                if j ~= 2 % not first entry
                    % Average Exon Calculations
                    exlist{j-1, 13} = gene_length;
                    exlist{j-1, 14} = gene_counts;
                    exlist{j-1, 15} = rpkm_sum/exon_num;
                    rpkm_sum = 0;
                    exon_num = 0;
                end
                gene_length = exlist{j, 6};
                gene_counts = exlist{j, 7};            
            else
                gene_length = gene_length + exlist{j, 6};
                gene_counts = gene_counts + exlist{j, 7};
            end            
            rpkm_sum = rpkm_sum + exlist{j,12};
            
            j = j+1;
            exlist{j, 1} = data{i, 1};
            exlist{j, 2} = data{i, 2};
            exlist{j, 3} = data{i, 3};
            exlist{j, 4} = data{i, 4};
            exlist{j, 5} = data{i, 5};
            exlist{j, 6} = exlist{j, 4} - exlist{j, 3} + 1;
            exlist{j, 7} = data{i, 7};
            exon_num = exon_num + 1;            
        else % same gene
            exlist{j, 3} = min(exlist{j, 3}, data{i, 3});
            exlist{j, 4} = max(exlist{j, 4}, data{i, 4});
            exlist{j, 7} = exlist{j, 7} + data{i, 7};
        end
    end
    exlist{j, 6} = exlist{j, 4} - exlist{j, 3} + 1;
    
    % Final data entry + calculations
    % Exon Level Calculations
    exlist{j, 8} = exlist{j, 6}/1000;
    exlist{j, 9} = exlist{j, 7}/exlist{j, 8};
    exlist{j, 10} = assigned_reads;
    exlist{j, 11} = assigned_reads/1000000;
    exlist{j, 12} = exlist{j, 9}/exlist{j, 11};
    rpkm_sum = rpkm_sum + exlist{j,12};
        
    if ~strcmp(exlist{j-1, 1}, exlist{j, 1})
        if j ~= 2
            exlist{j-1, 13} = gene_length;
            exlist{j-1, 14} = gene_counts;
            exlist{j-1, 15} = rpkm_sum/exon_num;
            rpkm_sum = 0;
            exon_num = 0;
        end
        gene_length = exlist{j, 6};
        gene_counts = exlist{j, 7};
        rpkm_sum = exlist{j,12};
        exon_num = exon_num + 1;        
    else        
        gene_length = gene_length + exlist{j, 6};
        gene_counts = gene_counts + exlist{j, 7};
        exon_num = exon_num + 1;
    end

    % Average Exon Calculations
    exlist{j, 13} = gene_length; 
    exlist{j, 14} = gene_counts;
    exlist{j, 15} = rpkm_sum/exon_num;
    
    exlist = exlist(~all(cellfun(@isempty, exlist), 2), :); % remove empty rows
    superexons = exlist(:, 1:12); % superexons and RPKM calculations (single sheet output)
    avg_exon = exlist(~cellfun(@isempty,exlist(:,13)), [1,13:15]); % average for each gene
