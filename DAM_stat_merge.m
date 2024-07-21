% Merge data generated by DAM_stat.m across different trials

tic;
clear;
close all;
warning off;
addpath('functions');

%% parameters
save_path = '/Users/apple/Desktop/MyRun14-15'; % path to save the merged data
save_name = 'merged_stat'; % name of the merged data file

% load data
data_sources = { % paths to the stat files generated by DAM_stat.m
    % '/Users/apple/Desktop/MyRun14/actograph_11260430-12121800';
    % '/Users/apple/Desktop/MyRun15_1/actograph_12170430-01021800';
    '/Users/apple/Desktop/141516/stat_01120430-01140400';
    '/Users/apple/Desktop/141516/stat_11280430-11300400';
    '/Users/apple/Desktop/141516/stat_12190430-12210400';
};

prompt = 'Merge mode (inner/outer): ';
merge_mode = questdlg(prompt, 'Merge mode', 'inner', 'outer', 'inner'); % 'inner' or 'outer' join mode for merging data, 'inner' is an intersection of all group names, 'outer' is a union of all group names
% merge_mode = 'inner';

picture_format = 'pdf'; % 'png' and 'pdf' are recommanded 
picture_resolution = 150; % resolution of the picture
picture_fontsize = 16;

Y_max_activity = 100; % Y-axis limit of the activity bar plot, set as 'auto' for automatic scaling

% bar plot parameters
bar_width = 1; % width of the bars in the bar plot
bar_facecolor = [0.5 0.5 0.5]; % color of the bars
bar_edgecolor = [0 0 0]; % edge color of the bars
bar_facecolor_middle = 'auto'; % 'auto' or color array/string

%% init
CheckDir(save_path);
time_stamp = string(datetime("now",'Format','yyyyMMddHHmmss'));
save_folder = fullfile(save_path, save_name + "_" + time_stamp);
mkdir(save_folder);

source_N = length(data_sources);
group_name = cell(source_N, 1);
data_duration = zeros(source_N, 1);
fprintf('Validating data sources...\n');
for i = 1:source_N
    fprintf('\t');
    CheckDir(data_sources{i});
    file_list = dir(fullfile(data_sources{i}, 'info_*.mat'));
    LL_info = load(fullfile(data_sources{i}, file_list(1).name));
    group_name{i} = LL_info.group_names;
    data_duration(i) = hours(datetime(LL_info.time_end) - datetime(LL_info.time_start)) + 0.5;
    if i == 1 && strcmp(bar_facecolor_middle, 'auto')
        light_on = LL_info.light_on;
        if ~isempty(light_on)
            bar_facecolor_middle = [1 1 1]; % color of the bars from T6 to T18
        else
            bar_facecolor_middle = [0.7 0.7 0.7];
        end
    end
end

% get intersection of group names
group_names = group_name{1};
if strcmp(merge_mode, 'inner')
    for i = 2:source_N
        group_names = intersect(group_names, group_name{i});
    end
elseif strcmp(merge_mode, 'outer')
    for i = 2:source_N
        group_names = union(group_names, group_name{i});
    end
else
    error('Invalid merge mode');
end

group_N = length(group_names);
if group_N == 0
    error('No common group names found');
end

% compare data duration
if length(unique(data_duration)) > 1
    error('Data lengths are not consistent across different trials');
end

% write data source to txt
fid = fopen(fullfile(save_folder, 'data_sources.txt'), 'w');
fprintf(fid, 'Data sources:\n');
for i = 1:source_N
    fprintf(fid, '\t%s\n', data_sources{i});
end
fprintf(fid, '\n');
fprintf(fid, 'merge_mode: %s\n', merge_mode);
fprintf(fid, '\n');
fprintf(fid, 'Group names:\n');
for g = 1:group_N
    fprintf(fid, '\t%s\n', group_names{g});
end
fprintf(fid, '\n');
fclose(fid);

%% merge activity_plot
fprintf('Valid group names: \n');
for g = 1:group_N
    fprintf('\t%s\n', group_names{g});
end
fprintf('Merging stat data...\n');
for g = 1:group_N
    act_table = table();
    for i = 1:source_N
        % load data

        % find the file starts with 'activity_plot_' in the folder
        file_list = dir(fullfile(data_sources{i}, 'group', 'activity_plot_*.xlsx'));
        plt_file = fullfile(data_sources{i}, 'group', file_list(1).name);

        T = readTableByKeyword(plt_file, group_names{g}+"_flies");
        if isempty(T)
            continue;
        end

        % cat table horizontally
        act_table = cat(2, act_table, T);
    end
    act_g = table2array(act_table);
    act_g = act_g(1:end-1,:); % remove the row of mean
    [fig_act, act_plt, err, N] = plotActivityBars(act_g, "facecolor", bar_facecolor, "facecolor_middle", bar_facecolor_middle, "bar_width", bar_width, "edgecolor", bar_edgecolor, "Ymax", Y_max_activity, "FontSize", picture_fontsize);
    title("Average activity of " + group_names{g} + " (N = " + N + ")");
    print(fig_act, fullfile(save_folder, group_names{g} + "_activity"), "-d"+picture_format, "-r"+picture_resolution);
    close(fig_act);
    table_act_plt = array2table([act_plt', err'], 'VariableNames', {'Mean', 'SEM'});
    writetable(table_act_plt, fullfile(save_folder, "activity_plot.xlsx"), 'Sheet', "g" + g + "_" + group_names{g}, 'WriteRowNames', true);
    writetable(act_table, fullfile(save_folder, "activity_plot.xlsx"), 'Sheet', "g" + g + "_" + group_names{g} + "_flies", 'WriteRowNames', true);
end

%% merge period data and death time
period_table = table();
death_time_table = table();
for g = 1:group_N
    period_table = table();
    for i = 1:source_N
        % merge period data
        file_list = dir(fullfile(data_sources{i}, 'group', 'period_*.xlsx'));
        plt_file = fullfile(data_sources{i}, 'group', file_list(1).name);

        T = readTableByKeyword(plt_file, group_names{g});
        if isempty(T)
            continue;
        end

        period_table = cat(1, period_table, T);

        % merge death time
        file_list = dir(fullfile(data_sources{i}, 'group', 'death_time_*.xlsx'));
        plt_file = fullfile(data_sources{i}, 'group', file_list(1).name);
        
        T = readTableByKeyword(plt_file, group_names{g});
        if isnumeric(T.DeathTime)
            T.DeathTime = num2cell(T.DeathTime);
        end
        
        if isempty(T)
            continue;
        end
        death_time_table = cat(1, death_time_table, T);
    end
    writetable(period_table, fullfile(save_folder, "period.xlsx"), 'Sheet', "g" + g + "_" + group_names{g}, 'WriteRowNames', true);
    writetable(death_time_table, fullfile(save_folder, "death_time.xlsx"), 'Sheet', "g" + g + "_" + group_names{g}, 'WriteRowNames', true);
end


%%
toc

function T = readTableByKeyword(file, keyword)
    % read table from the file that contains the keyword
    % file: path to the file
    % keyword: keyword to search for
    % return: table
    file_sheetnames = sheetnames(file);
    sheet_name = "";
    count = 0;
    for s = 1:length(file_sheetnames)
        if endsWith(file_sheetnames{s}, keyword)
            count = count + 1;
            if count == 1
                sheet_name = file_sheetnames{s};
            end
        end
    end
    if count > 1
        fprintf('Multiple sheets containing %s found in %s\n', keyword, file);
        for s = 1:length(file_sheetnames)
            fprintf('\t%s\n', file_sheetnames{s});
        end
        fprintf('The first one is read: %s\n', sheet_name);
    end
    if sheet_name == ""
        fprintf('Sheet name containing %s not found in %s\n', keyword, file);
        T = table();
    else
        T = readtable(file, 'Sheet', sheet_name, 'ReadRowNames', true);
    end
end
