function CheckDir(directory)
if ~exist(directory,"dir")
    error("Not Found Directory: " + directory);
end
disp("Current processing directory: " + directory);
end