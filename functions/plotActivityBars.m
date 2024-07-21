function [fig, act_avg, act_err, n_col, act_flies] = plotActivityBars(activity, options)
    % plot average activities to 0-24 h
    arguments
        activity; % T x fly_N vector, T should be a multiple of 48
        options.facecolor = [0.5 0.5 0.5];
        options.facecolor_middle = [1 1 1];
        options.bar_width = 1;
        options.edgecolor = 'k';
        options.Visible = 'off';
        options.Ymax = [];
        options.FontSize = [];
    end

    [n_row, n_col] = size(activity);
    fly_N = sum(all(~isnan(activity), 1)); % number of flies (exclude NaNs)
    if mod(n_row, 48) ~= 0
        % trim the last part of the data
        activity = activity(1:floor(n_row/48)*48, :);
        warning('Input should be a column vector with a length of multiple of 48. The last part of the data is trimmed.');
    end

    if n_col > 1 % if there are multiple flies, calculate the average activity across flies and days
        act_flies = zeros(n_col, 48);
        for nn = 1:n_col
            act_nn = transpose(reshape(activity(:, nn), 48, [])); % reshape the activity of each fly to day_Num x 48
            act_nn_avg = mean(act_nn, 1, 'omitnan'); % average activity of each fly across days
            act_flies(nn, :) = act_nn_avg;
        end
        act_avg = mean(act_flies, 1, 'omitnan'); % average activity across flies
        act_err = std(act_flies, 0, 1, 'omitnan') / sqrt(fly_N); % standard error of the mean, calculated across flies
    else % if there is only one fly, calculate the average activity across days
        act_i = transpose(reshape(activity,48,[]));
        act_avg = mean(act_i, 1);
        act_err = std(act_i, 0, 1) / sqrt(size(act_i, 1)); % calculated across days
    end

    fig = figure('Visible', options.Visible);
    hold on;

    xt = (0.5:0.5:24) - 0.25;
    bar(xt, act_avg, options.bar_width, 'FaceColor', options.facecolor, 'EdgeColor', options.edgecolor);

    % change the color of the bars from T6 to T18
    bar(xt(13:36), act_avg(13:36), options.bar_width, 'FaceColor', options.facecolor_middle, 'EdgeColor', options.edgecolor);

    % add SEM as error bar
    errorbar(xt, act_avg, 0, act_err, 'k', 'linestyle', 'none'); % only positive part of the error bar
    hold off;

    xlim([0 24]);
    xticks(6:6:24);
    xticklabels({'6', '12', '18', '24'});
    xlabel('Time (h)');
    ylabel('Activity');

    if isa(options.Ymax, 'double')
        ylim([0 options.Ymax]);
    end

    % remove top and right border
    set(gca, 'box', 'off');

    % set font size
    if ~isempty(options.FontSize)
        set(gca, 'FontSize', options.FontSize);
    end
    
end