function directory = DirDiagnostics(directory)
    % if the directory is not found, ask the user to select the folder manually
    try
        CheckDir(directory);
    catch

        fprintf(directory+" not found! Please select data folder manually!\n")
        directory = uigetdir(pwd, 'Select the folder containing the data');
        if directory == 0
            error("No folder selected!");
        end
        CheckDir(directory);
    end
end