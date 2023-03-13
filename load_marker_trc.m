function out = load_marker_trc(filename)
file = filename;
[file_data,s_data]= readtext(file, '\t', '', '', 'empty2NaN');
[m,n] = size(file_data);

%% get order of markers
markerOrder = file_data(4,:);
data_label{1} = 'Frame';
data_label{2} = 'Time';
for i = 3 : 3: numel(markerOrder)
    data_label{i} = [markerOrder{i} '_X'];
    data_label{i+1} = [markerOrder{i} '_Y'];
    data_label{i+2} = [markerOrder{i} '_Z'];
end


% %% find cell containing 'endheader' 
% exact_match_mask = strcmp(file_data, 'endheader');
% endheader_location = find(exact_match_mask);
% 
% %% find the column headings (row directly below endheader)
% data_label_index = endheader_location + 1;
% data_label = file_data(data_label_index,:);

%% create an array with all of the numerical data
num_data = cell(m - 6, n);
for row = 1 : m - 6
    for col = 1 : n
        num_data{row, col} = file_data{row+6, col};
    end
end


%% check for duplicate lables  
% go through the data labels and find any that are duplicates (this occurs
% in the ground reaction force data where each forceplate has the same
% column headings) and add a number to distinguish the duplicates.

for col = 1:length(data_label)
    tf = strcmp(data_label(col),data_label);
    c = find(tf>0);
    if length(c) > 1
        for j = 1:length(c)
            data_label(c(j)) = cellstr([data_label{c(j)} num2str(j)]);
        end
    end
end
%% create output structure 
% field names from  data labels
% corresponding data from  columns of the data array

data_label = strrep(data_label, '*', 'not_labelled_');

for col = 1:length(data_label)
    f_name = data_label{col};
    % find any spaces and replace with underscore
    e = strfind(f_name, ' ');
    if ~isempty(e)
        f_name(e) = '_';
    end
    e = strfind(f_name, '.');
    if ~isempty(e)
        f_name(e) = '_';
    end
    if ~isempty(str2num(f_name(1)))
        f_name = ['N' f_name];
    end
    out.(f_name) = num_data(:,col);
end