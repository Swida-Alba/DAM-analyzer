% Performs statistics and visualization of DAM data.
% Input data should be in the format of merged 30-min activities collected by DAM system, saved in .txt format.
% All files in the same data folder will be processed. The files should be collected in the same experiment.
% 
% Created by KRL, on Feb 27, 2024, E-mail: krleng@pku.edu.cn
% Last modified by KRL, on Mar 3, 2024
%
% Requirements:
%   MATLAB R2022b or later
%   MATLAB Parallel Computing Toolbox
%   MATLAB Statistics and Machine Learning Toolbox
%   functions used:
%       DirDiagnostics.m
%       CheckDir.m
%       CircadianPeriodogram.m
%       parsave.m
%       plotPeriodogram.m
%       RoundByMagnitude.m
%       GetDeathTimeInDAM.m
%       plotActivityBars.m
%       ReadFileDAM.m

tic;
clear;
close all;
warning off;
addpath('functions');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% input parameters

% data path 
datafolder = "/Users/apple/Desktop/MyRun14/"; % the folder containing the .txt files of the data.

% file parameters
% group = {
%     {1:16, 17:32};
%     {1:20, 21:32};
%     {[],   23:32};
%     {1:16, 17:32};
%     {1:20, 21:32};
%     {1:22,    []};
%     {1:22, 23:32};
%     {1:22, 23:32};
% };
group = {1:10,11:32};
group_names = {'WT', 'per0'};

time_start = '2023-11-26 4:30:00'; % 'yyyy-MM-dd HH:mm:ss' % start time of the data to be analyzed. This point is included in the analysis.
time_end   = '2023-11-28 4:00:00'; % 'yyyy-MM-dd HH:mm:ss' % end time of the data to be analyzed. This point is also included in the analysis.
light_on = '10:00:00'; % for morning/evening anticipation analysis
% light_on = ''; % for DD, set light_on as empty, and the code will not analyze morning/evening anticipation.
light_on_duration = 12; % light ON duration in hours

period_method = "greedy Chi-Square"; % "Enright" ("Chi-Square"), "Lomb-Scargle", "greedy Chi-Square", "Cosinor"
confidence_level = 0.9999994; % confidence level for statistics. % 0.95 for 2 sigma, 0.997 for 3 sigma, 0.99994 for 4 sigma, 0.9999994 for 5 sigma

isToExcludeDeath = true; % whether to exclude the flies that died during the time range.

picture_format = 'pdf'; % 'png' and 'pdf' are recommanded 
picture_resolution = 150; % resolution of the picture
picture_fontsize = 16;

Y_max_activity = 100; % Y-axis limit of the activity bar plot, set as 'auto' for automatic scaling

% bar plot parameters
bar_width = 1; % width of the bars in the bar plot
bar_facecolor = [0.5 0.5 0.5]; % color of the bars
bar_edgecolor = [0 0 0]; % edge color of the bars

if ~isempty(light_on)
    bar_facecolor_middle = [1 1 1]; % color of the bars from T6 to T18
else
    bar_facecolor_middle = [0.7 0.7 0.7];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% initialization

% check data directory
datafolder = DirDiagnostics(datafolder); % check if the folder exists. if not, you can reset it or select the folder manually in the dialog box.

% load data folder
filelist = dir(fullfile(datafolder, '*.txt'));
file_N = length(filelist);
if file_N == 0
    error('No .txt file found in the folder.');
end
file_1 = ReadFileDAM(fullfile(datafolder, filelist(1).name));
t_dataEnd = file_1.Datetime(end);

% get group information
group_N = length(group_names);
if ~isa(group{1}, 'cell')
    group = {group};
end

if length(group) == 1
    for g = 1:file_N
        group{g} = group{1};
    end
else
    if length(group) ~= file_N
        error('Multiple group indices are provided, but the number of groups does not match the number of files. Please check the `group` variable. file_N = %d, group_N = %d', file_N, length(group));
    end
end

% get time information and data range
t_start = datetime(time_start, 'Format', 'yyyy-MM-dd HH:mm:ss');
t_end = datetime(time_end, 'Format', 'yyyy-MM-dd HH:mm:ss');

day_N = days(t_end - t_start + hours(0.5));
if day_N < 8
    warning on;
    warning('The time range is less than 8 days. The period analysis may not be accurate.');
    warning off;
end

if t_start < file_1.Datetime(1) || t_start > file_1.Datetime(end) || t_end < file_1.Datetime(1) || t_end > file_1.Datetime(end)
    error("Selected time range exceeds the data range. The data range is from "+string(file_1.Datetime(1))+" to "+string(t_dataEnd)+".");
end

date_info = string(datetime(string(t_start), 'Format', 'MMdd')) + "-" + string(datetime(string(t_end), 'Format', 'MMdd'));
date_info_long = string(datetime(string(t_start), 'Format', 'MMddHHmm')) + "-" + string(datetime(string(t_end), 'Format', 'MMddHHmm'));
idx_range = find(file_1.Datetime >= t_start & file_1.Datetime <= t_end); % index of the time range

if isempty(idx_range)
    error('No data found in the selected time range. Please check `time_start` and `time_end` variables');
end

% get days in the time range for analysis of Morning/Evening anticipation
light_on_time = NaT;
if ~isempty(light_on)
    day_start = t_start;
    day_start.Format = 'yyyy-MM-dd';
    day_end = t_end;
    day_end.Format = 'yyyy-MM-dd';
    if t_start > datetime(string(day_start)+" "+light_on) - hours(5.5) % for the morning anticipation of the first day
        day_start = day_start + days(1);
    end
    if t_end < datetime(string(day_end)+" "+light_on) + hours(light_on_duration) % for the evening anticipation of the last day
        day_end = day_end - days(1);
    end

    days_series = day_start : (day_end+days(1)); % to include day_end, add one more day in this operation.
    % get light on time for each day if in LD cycles
    light_on_time = NaT(length(days_series), 1);
    for d = 1:length(days_series)
        light_on_d = datetime(string(days_series(d)) + " " + light_on, 'Format', 'yyyy-MM-dd HH:mm:ss');
        light_on_time(d) = light_on_d;
    end
    light_on_time.Format = 'yyyy-MM-dd HH:mm:ss';
end

% set up folders
folder_stat = fullfile(datafolder, 'stat_'+date_info_long);
folder_period = fullfile(folder_stat, 'period');
folder_bars = fullfile(folder_stat, 'activity_bars');
folder_data = fullfile(folder_stat, 'data');
folder_group = fullfile(folder_stat, 'group');
mkdir(folder_stat);
mkdir(folder_period);
mkdir(folder_bars);
mkdir(folder_data);
mkdir(folder_group);

% save parameters
save(fullfile(folder_stat, "info_"+date_info_long+".mat"), 'time_start', 'time_end', 'light_on', 'light_on_duration', 'light_on_time', 'group_names');

%% process data by file
fprintf('Processing files...\n');
parfor i = 1:file_N
    [~, filename, ext] = fileparts(filelist(i).name);
    file_i = ReadFileDAM(fullfile(datafolder, [filename, ext]));

    % extract activity data
    activity = file_i.data(idx_range,end-31:end); % the last 32 columns are monitor channels 1-32.
    activity_all = file_i.data(:,end-31:end); % all the data

    % get death time of each fly
    death_time = GetDeathTimeInDAM(file_i, t_start, 24); % get the death time of each fly
    if isToExcludeDeath
        activity(:, death_time <= t_end) = NaN; % exclude the flies that died during the time range.
        activity_all(:, death_time <= t_end) = NaN;
    end
    dead_idx = find(all(isnan(activity), 1)); % get the index of the dead flies
    survival_hours_since_time_start = hours(death_time - t_start);
    survival_hours = hours(death_time - file_i.Datetime(1));
    death_time = string(death_time);
    parsave(fullfile(folder_data, filename+"_activity_"+date_info+".mat"), activity, light_on, time_start, time_end, death_time, survival_hours_since_time_start, survival_hours, idx_range);

    % period analysis
    [stat, powerArray, confidenceLine, confidence, P_candi] = CircadianPeriodogram(activity', "Method", period_method, "ConfidenceLevel", confidence_level);
    parsave(fullfile(folder_data, filename+"_period_"+date_info+".mat"), stat, powerArray, confidenceLine, confidence, P_candi);

    % plot activity bars and periodogram for each fly in time range
    for fly_i = 1:32
        % skip the dead flies
        if ismember(fly_i, dead_idx)
            continue;
        end
        % get the assigned group of the fly
        group_idx = find(cellfun(@(x) ismember(fly_i, x), group{i}));
        if isempty(group_idx) % skip the flies not assigned to any group
            group_idx = 0;
            gName = "UnassignedGroup";
        else 
            gName = group_names{group_idx};
        end

        % plot activity bars
        title_name = "Average activity of " + gName + " fly " + fly_i;
        fig_name = "g" + group_idx + "_" + filename + "_fly" + fly_i + "_" + date_info;
        [fig_activity, act_plot, err] = plotActivityBars(activity(:, fly_i), "facecolor", bar_facecolor, "facecolor_middle", bar_facecolor_middle, "bar_width", bar_width, "edgecolor", bar_edgecolor, "Ymax", Y_max_activity, "FontSize", picture_fontsize);
        title(title_name);
        print(fig_activity, fullfile(folder_bars, fig_name), "-d"+picture_format, "-r"+picture_resolution);
        close(fig_activity);

        % period analysis and plot periodogram
        title_name = period_method + " Periodogram of " + gName + " fly " + fly_i;
        fig_name = "g" + group_idx + "_" + filename + "_fly" + fly_i + "_" + date_info;
        fig_periodogram = plotPeriodogram(P_candi, powerArray(fly_i,:), confidenceLine, confidence, 'ylabel', 'Statistic', 'titleName', title_name);
        print(fig_periodogram, fullfile(folder_period, fig_name), "-d"+picture_format, "-r"+picture_resolution);
        close(fig_periodogram);
    end

    % morning/evening anticipation analysis
    if ~isempty(light_on)
        MIs = NaN(length(days_series), 32); % morning anticipation index
        EIs = NaN(length(days_series), 32); % evening anticipation index
        for d = 1:length(light_on_time)
            light_on_d = light_on_time(d);
            light_off_d = light_on_d + hours(light_on_duration);

            morning_3h_activity = activity_all(file_i.Datetime > light_on_d - hours(3) & file_i.Datetime <= light_on_d, :);
            morning_6h_activity = activity_all(file_i.Datetime > light_on_d - hours(6) & file_i.Datetime <= light_on_d, :);
            evening_3h_activity = activity_all(file_i.Datetime > light_off_d - hours(3) & file_i.Datetime <= light_off_d, :);
            evening_6h_activity = activity_all(file_i.Datetime > light_off_d - hours(6) & file_i.Datetime <= light_off_d, :);

            MI_d = sum(morning_3h_activity, 1) ./ sum(morning_6h_activity, 1); % MI of each fly at day d
            EI_d = sum(evening_3h_activity, 1) ./ sum(evening_6h_activity, 1); % EI of each fly at day d

            MIs(d, :) = MI_d;
            EIs(d, :) = EI_d;
        end
        parsave(fullfile(folder_data, filename+"_anticipation_"+date_info+".mat"), MIs, EIs);
    end
    fprintf('File %d/%d processed.\n', i, file_N);

end

%% merge data by group
fprintf('Analyzing groups......\n');
for g = 1:group_N
    % group_flies = group{g};
    group_name = group_names{g};

    %%% merge activity data
    activity_g = [];
    data_source = {};
    for i = 1:file_N
        group_flies = group{i}{g};
        [~, filename, ~] = fileparts(filelist(i).name);
        LL_act = load(fullfile(folder_data, filename+"_activity_"+date_info+".mat"));
        activity_g = cat(2, activity_g, LL_act.activity(:, group_flies)); % concatenate activity data of the flies in the group, each column for a fly
        for fly_i = group_flies
            data_source = cat(2, data_source, filename+"_fly"+fly_i); % record the file source and fly NO. of the data
        end
    end

    % save activity data of the current group to .xlsx file
    row_names = string(file_1.Datetime(idx_range));
    table_activity = array2table(activity_g, 'VariableNames', data_source, 'RowNames', row_names, 'DimensionNames', {'Time', 'Fly'});
    writetable(table_activity, fullfile(folder_group, "activity_" + date_info + ".xlsx"), 'WriteRowNames', true, 'Sheet', "g" + g + "_" + group_name);

    % plot average activity of the group 
    [fig_act, act_plt, err, N, act_flies] = plotActivityBars(activity_g, "facecolor", bar_facecolor, "facecolor_middle", bar_facecolor_middle, "bar_width", bar_width, "edgecolor", bar_edgecolor, "Ymax", Y_max_activity, "FontSize", picture_fontsize);
    title("Average activity of " + group_name + " (N = " + N + ")");
    print(fig_act, fullfile(folder_group, "g" + g + "_" + group_name + "_activity_" + date_info), "-d"+picture_format, "-r"+picture_resolution);
    close(fig_act);
    table_act_plt = array2table([act_plt', err'], 'VariableNames', {'Mean', 'SEM'});
    writetable(table_act_plt, fullfile(folder_group, "activity_plot_" + date_info + ".xlsx"), 'Sheet', "g" + g + "_" + group_name);
    row_names = [string(0.5:0.5:24) + " hr", "Mean"]';
    act_flies = transpose(act_flies);
    table_act_flies = array2table([act_flies; mean(act_flies,1)], 'VariableNames', data_source, 'RowNames', row_names, 'DimensionNames', {'Time', 'Fly'});
    writetable(table_act_flies, fullfile(folder_group, "activity_plot_" + date_info + ".xlsx"), 'Sheet', "g" + g + "_" + group_name + "_flies", 'WriteRowNames', true);

    %%% merge period data
    period_g = [];
    peakPower_g = [];
    p_val_g = [];
    significance_g = [];
    for i = 1:file_N
        group_flies = group{i}{g};
        [~, filename, ~] = fileparts(filelist(i).name);
        LL_per = load(fullfile(folder_data, filename+"_period_"+date_info+".mat"));
        period_g = cat(1, period_g, LL_per.stat.Period(group_flies));
        peakPower_g = cat(1, peakPower_g, LL_per.stat.PeakPower(group_flies));
        p_val_g = cat(1, p_val_g, LL_per.stat.p_value(group_flies));
        significance_g = cat(1, significance_g, LL_per.stat.SignificanceFlag(group_flies));
    end

    % save period data of the current group to .xlsx file
    period_data = table(period_g, peakPower_g, p_val_g, significance_g, 'VariableNames', {'Period', 'PeakPower', 'p_value', 'SignificanceFlag'}, 'RowNames', data_source, 'DimensionNames', {'Fly', 'Variables'});
    writetable(period_data, fullfile(folder_group, "period_" + date_info + ".xlsx"), 'WriteRowNames', true, 'Sheet', "g" + g + "_" + group_name);

    %%% merge morning/evening anticipation data
    if ~isempty(light_on)
        light_off_time = light_on_time + hours(light_on_duration);
        light_off_time.Format = 'yyyy-MM-dd HH:mm:ss';
        MI_g = [];
        EI_g = [];
        for i = 1:file_N
            group_flies = group{i}{g};
            [~, filename, ~] = fileparts(filelist(i).name);
            LL_anticipation = load(fullfile(folder_data, filename+"_anticipation_"+date_info+".mat"));
            MI_g = cat(2, MI_g, mean(LL_anticipation.MIs(:, group_flies), 1, 'omitnan'));
            EI_g = cat(2, EI_g, mean(LL_anticipation.EIs(:, group_flies), 1, 'omitnan'));
        end

        % save morning/evening anticipation data of the current group to .xlsx file
        table_anticipation = array2table([MI_g', EI_g'], 'VariableNames', {'MI', 'EI'}, 'RowNames', data_source, 'DimensionNames', {'Fly', 'Variables'});
        writetable(table_anticipation, fullfile(folder_group, "anticipation_" + date_info + ".xlsx"), 'WriteRowNames', true, 'Sheet', "g" + g + "_" + group_name);
    end

    %%% merge death time data
    death_time_g = [];
    survival_hours_since_time_start_g = [];
    survival_hours_g = [];
    for i = 1:file_N
        group_flies = group{i}{g};
        [~, filename, ~] = fileparts(filelist(i).name);
        LL_death = load(fullfile(folder_data, filename+"_activity_"+date_info+".mat"));
        death_time_g = cat(1, death_time_g, LL_death.death_time(group_flies));
        survival_hours_since_time_start_g = cat(1, survival_hours_since_time_start_g, LL_death.survival_hours_since_time_start(group_flies));
        survival_hours_g = cat(1, survival_hours_g, LL_death.survival_hours(group_flies));
    end

    % save death time data of the current group to .xlsx file
    table_death = table(death_time_g, survival_hours_since_time_start_g, survival_hours_g, 'VariableNames', {'DeathTime', 'SurvivalHoursSinceStart', 'SurvivalHours'}, 'RowNames', data_source, 'DimensionNames', {'Fly', 'Variables'});
    writetable(table_death, fullfile(folder_group, "death_time_" + date_info + ".xlsx"), 'WriteRowNames', true, 'Sheet', "g" + g + "_" + group_name);
    
end

%% timer
toc
