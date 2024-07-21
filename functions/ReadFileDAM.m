function file_struct = ReadFileDAM(file_path)
    file_i = importdata(file_path);
    for t = 1:length(file_i.textdata)
        date_day = datetime(file_i.textdata{t,2},'InputFormat', 'dd MMM yy', "Format",'yyyy-MM-dd', 'Locale', 'en_US');
        file_i.textdata{t,2} = string(date_day); % change the date format from 'dd MMM yy' to 'yyyy-MM-dd'
    end
    file_i.Datetime = datetime(file_i.textdata(:,2) + " " + file_i.textdata(:,3), 'Format', 'yyyy-MM-dd HH:mm:ss');
    file_struct = file_i;
end