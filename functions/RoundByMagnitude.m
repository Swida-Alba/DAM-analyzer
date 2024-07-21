function num = RoundByMagnitude(varargIn,Method)
arguments
    varargIn;
    Method = 'ceil';
end

mag = fix(log10(abs(varargIn)));
res = varargIn ./ (10.^mag);
switch Method
    case 'ceil'
        res = ceil(res);
    case 'floor'
        res = floor(res);
    case 'fix'
        res = fix(res);
    case 'round'
        res = round(res);
end
num = res.*(10.^mag);
end