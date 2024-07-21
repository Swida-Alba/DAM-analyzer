function fig_peaktime = plotPeakDrift(peaktime_plt, peaktime_err, colorPatch, options)
    arguments
        peaktime_plt; % a k x N array containing the mean of peak time of k days, N groups
        peaktime_err; % corresponding k x N array, or k x 2 x N array containing the negative and positive errors
        colorPatch;
        options.GroupNames = []
        options.FontSize = 16;
        options.Palettes = {'#000000', '#D04848', '#6895D2', '#F3B95F', '#1DA050','#6040D0'};
        options.LineWidth = 2;
        options.MarkerSize = 10;
        options.Marker = 'square';
        options.xtick = 0:4:24;
        options.AxisVisible = 'on';
    end

    [day_N, group_N] = size(peaktime_plt);
    if ismatrix(peaktime_err)
        peaktime_err = cat(3, peaktime_err, peaktime_err);
        peaktime_err = permute(peaktime_err, [1,3,2]);
    end
    if isempty(options.GroupNames)
        group_names = cell(1, group_N);
        for i = 1:group_N
            group_names{i} = "Group "+i;
        end
    end

    % trim color patches to match the available activity data length if necessary
    for i = 1:length(colorPatch)
        if size(colorPatch{i}.HalfMap, 1) > day_N
            colorPatch{i}.HalfMap = colorPatch{i}.HalfMap(1:day_N,:);
        end
    end

    %%% plot peak drifts
    fig_peaktime = figure('Units','normalized','Position',[0.05 0.05 0.45 0.8]);
    hold on;
    set(gca, 'YDir', 'reverse', 'Position', [0.1 0.1 0.8 0.85]);
    set(gca, 'FontSize', options.FontSize);

    % plot color patches
    for i = 1:length(colorPatch)
        patch_map = colorPatch{i}.HalfMap;
        [n_row,n_col] = size(patch_map);
        patch_alphas = unique(patch_map);
        for j = 1:length(patch_alphas)
            if patch_alphas(j) == 0 || isnan(patch_alphas(j))
                continue;
            end
            x = [];
            y = [];
            for rr = 1:n_row
                yy = [rr, rr-1, rr-1, rr]';
                for cc = 1:n_col
                    xx = [cc-1, cc-1, cc, cc]' / 2;
                    if patch_map(rr,cc) == patch_alphas(j)
                        x = cat(2, x, xx);
                        y = cat(2, y, yy);
                    end
                end
            end
            patch(x, y, colorPatch{i}.Color, 'EdgeColor', 'none', 'FaceAlpha', patch_alphas(j));
        end
    end

    % plot peak time
    graphic_objs = gobjects(1, group_N);
    for g = 1:group_N
        graphic_objs(g) = errorbar(peaktime_plt(:,g), 1:day_N, peaktime_err(:,1,g), peaktime_err(:,2,g), 'horizontal', 'Marker',options.Marker, 'Color', options.Palettes{g}, 'MarkerFaceColor', options.Palettes{g}, 'LineWidth', options.LineWidth, 'MarkerSize', options.MarkerSize);
    end
    hold off;
    xlim([0 24]);
    ylim([0 day_N+1]);
    xlabel('Time (h)', 'FontSize', options.FontSize);
    ylabel('Day', 'FontSize', options.FontSize);
    xticks(options.xtick);
    set(gca, 'TickDir', 'out');
    axis(options.AxisVisible);
    legend(graphic_objs, options.GroupNames, 'Location', 'southeastoutside');
end