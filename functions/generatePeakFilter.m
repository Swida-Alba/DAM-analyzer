function pfs = generatePeakFilter(xlsxfilepath, options)
    arguments
        xlsxfilepath = ''; % the path of the xlsx file
        options.BlockNumPerDay = 48; % number of 30-min blocks per day
        options.LengthPerBlock = 0.5; % in hour
    end
    % read xlsx file sheet by sheet to generate peak filter
    sn = sheetnames(xlsxfilepath);
    sheetN = length(sn);
    pfs = cell(sheetN, 1);
    for i = 1:sheetN
        T = readtable(xlsxfilepath, 'Sheet', sn{i}, 'VariableNamingRule','preserve');
        filter_raw = table2array(T);
        filter_raw = filter_raw(:, 2:end); % remove the first column
        dayN = size(filter_raw, 1);
        pf = zeros(dayN, options.BlockNumPerDay);
        for j = 1:dayN
            t0 = filter_raw(j,1);
            t1 = filter_raw(j,2);
            block0 = floor(t0/options.LengthPerBlock);
            block1 = ceil(t1/options.LengthPerBlock);
            block0 = max(block0, 1);
            block1 = min(block1, options.BlockNumPerDay);
            pf(j, block0:block1) = 1;
            pfs{i} = pf;
        end
    end
end