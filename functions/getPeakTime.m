function [varargout] = getPeakTime(activity_groups, options)
    arguments 
        activity_groups;
        options.BlockNumPerDay = 48; % number of 30-min blocks per day
        options.ReserveLastRow = false; % whether to reserve the last row of the actograph
        options.TrimActivityHead = 3; % in hour; Trimmed duration of the first 'D' in protocol to plot. For example, the original duration of the first 'D' is 6 h, and the TrimActivityHead is 2 h, then the left duration of the first 'D' is 4 h in the plot.
        options.Method = 'MeanOfMax'; % 'MeanOfMax', 'MaxOfMean', 'Cosinor'
        options.EveningFilter = [];
        options.ColorPatch = [];
    end

    group_N = length(activity_groups);

    % pre-process the activity data
    for g = 1:group_N
        activity_g = activity_groups{g};
        activity_g = activity_g(((options.TrimActivityHead)*2+1):end, :);
        activity_g = smoothdata(activity_g, 1, "loess", [3 3]);
        data_len = size(activity_g, 1);
        dayN = floor(data_len / options.BlockNumPerDay);
        data_len = dayN * options.BlockNumPerDay;
        activity_groups{g} = activity_g(1:data_len, :);
    end

    % filter activities by the LD pattern
    if ~isempty(options.EveningFilter)
        efilters = options.EveningFilter;
        len_f = length(efilters);
        if len_f == 1 % if only one evening filter is provided, use it for all groups
            efs = cell(group_N, 1);
            for g = 1:group_N
                efs{g} = options.EveningFilter{1};
            end
            efilters = efs;
        elseif len_f < group_N % if the number of evening filters is not equal to the number of groups, raise an error
            error('The number of evening filters should be 1 or equal to the number of groups.');
        elseif len_f > group_N
            warning('The number of evening filters is greater than the number of groups.')
        end

        for g = 1:group_N
            efilter = efilters{g};
            color_seq = reshape(transpose(efilter), [], 1);
            if length(color_seq) < data_len
                error('The length of the evening filter is shorter than the length of the activity data.');
            elseif length(color_seq) > data_len
                warning('The length of the evening filter is longer than the length of the activity data. Truncate the evening filter.');
            end
            color_seq = color_seq(1:data_len);
            activity_g = activity_groups{g};
            activity_g = activity_g .* color_seq; % filter the activities by the shifted color sequence
            activity_groups{g} = activity_g;
        end
    end

    
    peaktime_plt = [];  % peak activities of days averaged over flies, each group as a column
    peaktime_err = []; % standard error of the peak activities over flies, each group as a column
    peaktime_groups = cell(1, group_N);
    for g = 1:group_N
        activity_g = activity_groups{g};

        % if all activities are NaN, skip the group and fill the peaktime with NaN
        if all(isnan(activity_g), 'all') 
            peaktime_g = nan(size(activity_g));
            pt_avg = nan(dayN, 1);
            pt_err = nan(dayN, 1);
            peaktime_plt = cat(2, peaktime_plt, pt_avg);
            peaktime_err = cat(2, peaktime_err, pt_err);
            peaktime_groups{g} = peaktime_g;
            continue;
        end

        if strcmp(options.Method, 'MaxOfMean') % bootstrap
            boot_N = 2000;
            g_means = nan(size(activity_g, 1), boot_N);
            for b = 1:boot_N
                boot_idx = randi(size(activity_g, 2), 1, size(activity_g, 2));
                g_means(:,b) = mean(activity_g(:, boot_idx), 2, 'omitnan');
            end
            activity_g = g_means;
        end
        peaktime_g = [];
        for idx = 1:size(activity_g, 2)
            if all(isnan(activity_g(:,idx)))
                continue;
            end
            act_i = activity_g(:, idx);
            % reshape activity for actograph
            [~, act_fold] = ReshapeForActograph(act_i, 'BlockNumPerDay', options.BlockNumPerDay);
            act_fold(end,:) = []; % remove the last row of the actograph, which is the padding of NaNs

            % trim the last row if necessary
            if ~options.ReserveLastRow
                % if exists nan in the last row, remove the last row
                if any(isnan(act_fold(end,:)))
                    act_fold = act_fold(1:end-1,:);
                end
            end
            
            % get peak time of each day
            if strcmp(options.Method, 'MeanOfMax') || strcmp(options.Method, 'MaxOfMean')
                [~, peak_at] = max(act_fold,[],2);
                peak_at = peak_at / 2; % convert index of 1:48 of 30-min bins to 0.5:24 hour.
            elseif strcmp(options.Method, 'Cosinor')
                t = 0.5:0.5:24;
                day_N = size(act_fold, 1);
                peak_at = nan(day_N, 1);
                P = 24;
                for i = 1:day_N
                    [stat, ~, ~, phi] = cosinor(act_fold(i,:), t, P);
                    if stat.p < 1
                        peak_at(i) = (phi+pi) / (2 * pi) * P; % phi falls at the morning peak, so add pi to phi
                    end
                end
            end
            peaktime_g = cat(2, peaktime_g, peak_at);
        end
        if ~strcmp(options.Method, 'MaxOfMean')
            gfly_N = size(peaktime_g, 2) - sum(all(isnan(peaktime_g), 1));
            pt_avg = mean(peaktime_g, 2, 'omitnan');
            pt_err = std(peaktime_g, 0, 2, 'omitnan') / sqrt(gfly_N);
            peaktime_err = cat(2, peaktime_err, pt_err);
        else
            pt_avg = mean(peaktime_g, 2, 'omitnan');

            %%% standard error of the peak activities over flies
            pt_err = std(peaktime_g, 0, 2, 'omitnan');
            peaktime_err = cat(2, peaktime_err, pt_err);

            % %%% 95% confidence interval of the peak activities over flies
            % pt_interval = prctile(peaktime_g, [2.5, 97.5], 2);
            % pt_err = [pt_avg - pt_interval(:,1), pt_interval(:,2) - pt_avg];
            % peaktime_err = cat(3, peaktime_err, pt_err);
        end
        peaktime_plt = cat(2, peaktime_plt, pt_avg);
        peaktime_groups{g} = peaktime_g;
    end
    varargout = {peaktime_plt, peaktime_err, peaktime_groups};

end
