% Draw double-plot actograph of DAM data
% Modified from DAM_stat.m and HE_visualstat_3.m
% Input data should be in the format of merged 30-min activities collected by DAM system, saved in .txt format.
% All files in the same data folder will be processed. The files should be collected in the same experiment.
%
% Created by KRL, on Feb 29, 2024. E-mail: krleng@pku.edu.cn
% Last modified by KRL, on Mar 3, 2024
%
% Requirements:
%   MATLAB R2022b or later
%   MATLAB Parallel Computing Toolbox
%   functions used:
%       actograph.m
%       DirDiagnostics.m
%       CheckDir.m
%       ReadFileDAM.m
%       ReshapeForActograph.m
%       GetDeathTimeInDAM.m
%       getPeakTime.m
%       plotPeakDrift.m
%       AddPatchMap.m
%       parsave.m

tic
clear
close all
warning on
addpath('functions');

%% file parameters

% data path 
datafolder = '/Users/apple/Desktop/0618/MyRun7'; % the folder containing the .txt files of the data.

% group information
% group = {
%     {[1 9 17 25 26], [2 3 4 10 18 19 20], [11 12 27], [15 16 30 31], [13 14 28 29], [5:8, 21:24]};
%     {[1 9 17 25 26], [2 3 4 10 18 19 20], [11 12 27], [15 16 30 31], [13 14 28 29], [5:8, 21:24]};
%     {[1 9 17 25 26], [2 3 4 10 18 19 20], [11 12 27], [15 16 30 31], [13 14 28 29], [5:8, 21:24]};
%     {[1 9 17 25 26], [2 3 4 10 18 19 20], [11 12 27], [15 16 30 31], [13 14 28 29], [5:8, 21:24]};
%     {[1 9 17 25 26], [2 3 4 10 18 19 20], [11 12 27], [15 16 30 31], [13 14 28 29], [5:8, 21:24]};    
% };
% group_names = {'WT','CHO', 'norpA', 'cry02', 'norpA-cry02', 'norpA-CHO'};

% group = {[1 8 9 17 25], [2 3 12 13 18 19], [20 21 24 32]};
% group_names = {'6586', '3917', '2425'};

group = {
{[23 24],[13 14 15 16]};
{[23 24],[13 14 15 16]};
{[23 24],[13 14 15 16]};
{[23 24],[13 14 15 16]};
{[23 24],[13 14 15 16]};
};

group_names = {'WT','13158'};



time_start = '2023-07-27 04:30:00'; % 'yyyy-MM-dd HH:mm:ss' 
% start time of the data to be analyzed. This point is included in the analysis. 
% Usually, the start time is 6 h before the light ON of LD1

light_on_duration = 12; % light on duration in hours

picture_format = 'png'; % 'png' and 'pdf' are recommanded for pixel and vector graphics, respectively.
picture_resolution = 150; % resolution of the picture

Y_max_actograph = 70; % maximum value of the y-axis in the actograph plot

isToExcludeDeath = true; % whether to exclude the flies that died during the time range.

peakline_colors = {'#000000', '#D04848', '#6895D2', '#F3B95F', '#1DA050', '#6040D0', '#FF5733', '#C70039', '#900C3F', '#581845', '#FFC300', '#DAF7A6'};

%%% paramters for actograph plot
protocol = ['D' 'LDLDLDLD' 'DDDDDDDDDDDDDDDDDDDDDDDD'];
shift_day = 4; % shift the light/dark cycle at the end of the day;
shift_time = 0; % shift time in hours; positive for delay, negative for advance.

colorNum = 1;
palette = {[0.6 0.6 0.6]}; % color palette

D_time = 6; % in hour; duration of the first 'D' in protocol to plot.

%%% parameters for peak drift plot
peaktime_file = '';
% peaktime_file = '/Users/apple/Desktop/MyRunl03/actograph_03110230-03201100/peaktime_file_0311-0320.xlsx'; % the file containing the peak time data and err for peak drift plot. If no existed data provided, set it to ''.
% This peaktime_file should be in the format of xlsx file with two sheets: 'peaktime' and 'err'.
% Each column in the 'peaktime' sheet should contain the peak time data of one group, and the 'err' sheet should match the 'peaktime' sheet in size.

peaktime_method = 'MeanOfMax'; % 'MeanOfMax', 'MaxOfMean', 'Cosinor'
left_pad = 1; % the number of D hours to pad the left of the L's for peak drift plot. no greater than D_time.
fontsize = 16;

evening_filter_file = fullfile(datafolder, 'peakfilter_example.xlsx');
% evening_filter_file = '';

printIndividual = false;
printGroups     = true;


%% protocol and color patch
%%%%% PROTOCOL AND COLOR PATCH %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~isa(protocol, 'char')
    error('The `protocol` should be a character array.');
end

protocol(protocol == ' ') = [];
protocol_len = length(protocol);
colorPatch = cell(1,colorNum);

for i = 1:length(colorPatch)
    colorPatch{i}.Color = palette{i};
    colorPatch{i}.ProtocolMat = repmat([12; 0], 1, protocol_len);
    colorPatch{i}.ProtocolMat(1,1) = D_time; % 6-hr Dark plot before the first LD cycle.
end

for i = 2:protocol_len % set light/dark duration
    if protocol(i) == 'D'
        colorPatch{1}.ProtocolMat(1,i) = 24 - light_on_duration;
    elseif protocol(i) == 'L'
        colorPatch{1}.ProtocolMat(1,i) = light_on_duration;
    end
end

for i = 1:protocol_len
    if protocol(i) == 'D'
        colorPatch{1}.ProtocolMat(2,i) = 1;
    elseif protocol(i) == 'd'
        colorPatch{1}.ProtocolMat(2,i) = 0.9;
    elseif protocol(i) == 'L'
        colorPatch{1}.ProtocolMat(2,i) = 0;
    end
end

% delay or advance of light/dark patches

if length(shift_day) ~= length(shift_time)
    error('The length of shift_day and shift_time should be the same.');
end

for i = 1:length(shift_day)
    colorPatch{1}.ProtocolMat(1, shift_day(i)*2+1) = colorPatch{1}.ProtocolMat(1, shift_day(i)*2+1) + shift_time(i);
end

%%%%% END OF INPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% generate patch map
colorPatch = AddPatchMap(colorPatch);

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
    for i = 1:file_N
        group{i} = group{1};
    end
else
    if length(group) ~= file_N
        error('Multiple group indices are provided, but the number of groups does not match the number of files. Please check the `group` variable. file_N = %d, group_N = %d', file_N, length(group));
    end
end

% check if the group indices are repeated among groups
for i = 1:file_N
    group_set = group{i};
    all_index = [group_set{:}];
    unique_index = unique(all_index);
    if length(all_index) ~= length(unique_index)
        error('The group indices are repeated among groups. Please check the `group` variable.');
    end
end

% get time range
protocol_time = sum(colorPatch{1}.ProtocolMat(1,:));
t_start = datetime(time_start, 'Format', 'yyyy-MM-dd HH:mm:ss');
t_end = t_start + hours(protocol_time - 0.5);
if t_end > t_dataEnd
    fprintf("Selected time range exceeds the data range. The end time is set to "+string(t_dataEnd)+".\n");
    t_end = t_dataEnd;
end

if t_start < file_1.Datetime(1) || t_start > file_1.Datetime(end) || t_end < file_1.Datetime(1) || t_end > file_1.Datetime(end)
    error("Selected time range exceeds the data range. The data range is from "+string(file_1.Datetime(1))+" to "+string(t_dataEnd)+".");
end

date_info = string(datetime(string(t_start), 'Format', 'MMdd')) + "-" + string(datetime(string(t_end), 'Format', 'MMdd'));
date_info_long = string(datetime(string(t_start), 'Format', 'MMddHHmm')) + "-" + string(datetime(string(t_end), 'Format', 'MMddHHmm'));
idx_range = find(file_1.Datetime >= t_start & file_1.Datetime <= t_end); % index of the time range

if isempty(idx_range)
    error('No data found in the selected time range. Please check `time_start` and `protocol` variables.');
end

% set up folders
folder_actograph = fullfile(datafolder, 'actograph_'+date_info_long);
folder_individual = fullfile(folder_actograph, 'individual');
mkdir(folder_actograph);
mkdir(folder_individual);

% read existed data for peak drift plot from excel
if ~isempty(peaktime_file) 
    if ~isfile(peaktime_file)
        error('The peaktime_file does not exist. Please check the file path.');
    end
    fprintf("Reading peak time data from "+peaktime_file+"...\n");
    peaktime_plt = readtable(peaktime_file, 'Sheet', 'peaktime_plt', 'ReadVariableNames', false, 'ReadRowNames', false);
    peaktime_err = readtable(peaktime_file, 'Sheet', 'peaktime_err', 'ReadVariableNames', false, 'ReadRowNames', false);
    peaktime_plt = table2array(peaktime_plt);
    peaktime_err = table2array(peaktime_err);
end

% save parameters
save(fullfile(folder_actograph, "info_"+date_info_long+".mat"), ...
        'group', 'group_names', 'time_start', 'light_on_duration', 'protocol', ...
        'shift_day', 'shift_time', 'colorPatch', 'isToExcludeDeath', ...
        'peaktime_file', 'D_time', 'left_pad');

%% plot actograph for individual flies in each file
fprintf("Analyzing DAM data files...\n");
activity_files = cell(1, file_N);
for i = 1:file_N
    [~, filename, ext] = fileparts(filelist(i).name);
    file_i = ReadFileDAM(fullfile(datafolder, [filename, ext]));
    
    %%% extract activity data
    activity = file_i.data(idx_range,end-31:end); % the last 32 columns are monitor channels 1-32.
    
    %%% plot actograph for each fly
    if printIndividual
        for fly_i = 1:32
            % get the assigned group of the fly
            group_idx = find(cellfun(@(x) ismember(fly_i, x), group{i}));
            if isempty(group_idx) % skip the flies not assigned to any group
                group_idx = 0;
                gName = "UnassignedGroup";
            else 
                gName = group_names{group_idx};
            end
    
            title_name = {
                "Activity Heatmap (g"+group_idx+" file"+i+" fly"+fly_i+")", ...
                "Activity Trendline", ...
                "Activity Actograph (g"+group_idx+" file"+i+" fly"+fly_i+")"
                };
            fig_name = "g" + group_idx + "_" + filename + "_fly" + fly_i;
            activity_plt = ReshapeForActograph(activity(:,fly_i));
            actAcg = actograph(activity_plt, 'Ymax', Y_max_actograph, 'titleName', title_name, 'ColorPatch', colorPatch);
            if strcmp(picture_format, 'pdf')
                print(actAcg.heatmap, fullfile(folder_individual, fig_name + "_heatmap_" + date_info), "-d"+picture_format, "-r"+picture_resolution, "-bestfit");
                print(actAcg.actogram, fullfile(folder_individual, fig_name + "_actograph_" + date_info), "-d"+picture_format, "-r"+picture_resolution, "-bestfit");
            else
                print(actAcg.heatmap, fullfile(folder_individual, fig_name + "_heatmap_" + date_info), "-d"+picture_format, "-r"+picture_resolution);
                print(actAcg.actogram, fullfile(folder_individual, fig_name + "_actograph_" + date_info), "-d"+picture_format, "-r"+picture_resolution);
            end
            close all;
        end
    end

    %%% Exclude the flies that died during the time range for group analysis
    death_time = GetDeathTimeInDAM(file_i, t_start, 24); % get the death time of each fly
    if isToExcludeDeath
        activity(:, death_time <= t_end) = NaN; % exclude the flies that died during the time range.
    end
    activity_files{i} = activity;

    fprintf("Individuals in file "+i+" processed.\n");
end

%% plot actograph for groups
fprintf("Processing groups......\n");
activity_groups = cell(1, group_N);
for g = 1:group_N
    group_name = group_names{g};

    activity_g = [];
    for i = 1:file_N
        group_flies = group{i}{g};
        activity = activity_files{i};
        activity_g = cat(2, activity_g, activity(:,group_flies));
    end
    dead_N = sum(all(isnan(activity_g), 1));
    N_val = size(activity_g, 2) - dead_N;
    activity_groups{g} = activity_g;
    activity_plt = ReshapeForActograph(mean(activity_g, 2, 'omitnan'));
    fig_name = "g" + g + "_" + group_name;
    if all(isnan(activity_g), 'all')
        fprintf("All flies in Group " + g + " " + group_name + " are dead.\n");
    end
    

    % plot actograph
    if printGroups
        title_name = {
            "Activity Heatmap of "+group_name+" (N = "+N_val+")", ...
            "Activity Trendline", ...
            "Activity Actograph of "+group_name+" (N = "+N_val+")"
            };
        
        actAcg = actograph(activity_plt, 'Ymax', Y_max_actograph, 'titleName', title_name, 'ColorPatch', colorPatch);
        if strcmp(picture_format, 'pdf')
            print(actAcg.heatmap, fullfile(folder_actograph, fig_name + "_heatmap_" + date_info), "-d"+picture_format, "-r"+picture_resolution, "-bestfit");
            print(actAcg.actogram, fullfile(folder_actograph, fig_name + "_actograph_" + date_info), "-d"+picture_format, "-r"+picture_resolution, "-bestfit");
        else
            print(actAcg.heatmap, fullfile(folder_actograph, fig_name + "_heatmap_" + date_info), "-d"+picture_format, "-r"+picture_resolution);
            print(actAcg.actogram, fullfile(folder_actograph, fig_name + "_actograph_" + date_info), "-d"+picture_format, "-r"+picture_resolution);
        end
        close all;
    end
    parsave(fullfile(folder_actograph, fig_name + "_data_" + date_info + ".mat"), activity_plt, N_val, activity_g);
    fprintf("Group "+g+" processed.\n");

end

%% plot peak time drifts and distributions

trim_time = D_time - left_pad;
if trim_time < 0
    error('In the plotting of peak drifts, left_pad should be less than D_time. Please check the two variables.');
end

% generate patch map for peak drifts
colorPatch_peak = colorPatch;
for i = 1:length(colorPatch_peak)
    colorPatch_peak{i}.ProtocolMat(1,1) = colorPatch_peak{i}.ProtocolMat(1,1) - trim_time;
end
colorPatch_peak = AddPatchMap(colorPatch_peak);

if exist(evening_filter_file, 'file')
    evening_filter = generatePeakFilter(evening_filter_file);
    fprintf("Evening filter loaded from --- "+evening_filter_file+".\n");
else

    evening_filter = [];
    fprintf("No evening filter applied.\n");
end

% plot peak drifts
peaktime_groups = [];
if isempty(peaktime_file)
    [peaktime_plt, peaktime_err, peaktime_groups] = getPeakTime(activity_groups, ...
    'TrimActivityHead', trim_time, ...
    'BlockNumPerDay', 48, ...
    'ReserveLastRow', false, ...
    'Method', peaktime_method, ...
    'EveningFilter', evening_filter, ...
    'ColorPatch', colorPatch_peak);
end

% save peak time data to .xlsx file for future use
peak_plt_table = array2table(peaktime_plt, 'VariableNames', group_names);
peak_err_table = array2table(peaktime_err, 'VariableNames', group_names);
filename_xlsx = "peaktime_file_"+date_info;
if exist(fullfile(folder_actograph, filename_xlsx+".xlsx"), 'file')
    warning("The file "+filename_xlsx+".xlsx already exists. A new file will be saved as "+filename_xlsx+"_new.xlsx.")
    filename_xlsx = filename_xlsx + "_new";
end
writetable(peak_plt_table, fullfile(folder_actograph, filename_xlsx+".xlsx"), 'Sheet', 'peaktime_plt', 'WriteMode', 'overwrite');
writetable(peak_err_table, fullfile(folder_actograph, filename_xlsx+".xlsx"), 'Sheet', 'peaktime_err', 'WriteMode', 'overwritesheet');
for g = 1:group_N
    peaktime_g = array2table(peaktime_groups{g});
    writetable(peaktime_g, fullfile(folder_actograph, filename_xlsx+".xlsx"), 'Sheet', "g"+g+"_"+group_names{g}, 'WriteMode', 'overwritesheet');
end

save(fullfile(folder_actograph, "peak_estimates_"+date_info+".mat"), 'peaktime_plt', 'peaktime_err', 'peaktime_groups', 'colorPatch_peak', 'peaktime_method', 'evening_filter', 'evening_filter_file');

fig_pt = plotPeakDrift(peaktime_plt, peaktime_err, colorPatch_peak, 'GroupNames', group_names, 'FontSize', fontsize, 'AxisVisible', 'on', 'Palettes', peakline_colors);
if strcmp(picture_format, 'pdf')
    print(fig_pt, fullfile(folder_actograph, "PeakDrift_"+date_info), "-d"+picture_format, "-r"+picture_resolution, "-bestfit");
else
    print(fig_pt, fullfile(folder_actograph, "PeakDrift_"+date_info), "-d"+picture_format, "-r"+picture_resolution);
end
% close(fig_pt);


%% timer
toc;

%% function 







