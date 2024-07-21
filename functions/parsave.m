function parsave(fname, varargin)
varN = length(varargin);
var_names = cell(1,varN);
for i = 1:varN
    eval(sprintf('%s=varargin{%d};',inputname(i+1),i));
    var_names{i} = inputname(i+1);
end
save(fname,var_names{:});
end