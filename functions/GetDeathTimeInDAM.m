function death_time = GetDeathTimeInDAM(file_i, t_start, death_threshold)
    arguments
        file_i; % file_i should be a struct with fields: data, textdata, which are the output of importdata function when reading a .txt file.
        t_start = []; % start time in <datetime>
        death_threshold = 24; % hours, if the fly is inactive for [more than] this time, it is considered dead.
    end

    if ~isempty(t_start) % find death start from given time to the file end.
        idx_range = find(file_i.Datetime > t_start);
        data_raw = file_i.data(idx_range, end-31:end); % all activity data in the file
        datetime_range = file_i.Datetime(idx_range);
    else
        data_raw = file_i.data(:, end-31:end);
        datetime_range = file_i.Datetime;
    end

    death_time = NaT(32, 1);
    for fly_i = 1:32
        die_at = NaT;
        for t = 1:size(data_raw, 1) - (death_threshold * 2 - 1)
            if all(data_raw(t:t+(death_threshold*2-1), fly_i) == 0)
                die_at = datetime_range(t);
                break;
            end
        end
        death_time(fly_i) = die_at;
    end
    death_time.Format = 'yyyy-MM-dd HH:mm:ss';
end