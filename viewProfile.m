function viewProfile(varargin)
% this function views the performance profile saved in resultsDir
% usage: viewProfile -f=<resultsDir>
% --Longfei Li
resultsDir='.'; % default dir is pwd
for i=1:nargin
    line = varargin{i};
    if(strncmp(line,'-f=',3))
        resultsDir=line(4:end); 
    end
end

load(sprintf('%s/profile.mat',resultsDir));
profview(0,pf);

end