function runShell(varargin)
%=========================================================================
% This is the main interface for setting up parameters and running various
% cases for the shell paper
% usage:
% runShell -options
%=========================================================================
infoPrefix = '--runShell--: '; % all info displayed by this function includes this prefix

%--------------- default options -----------------
parameters.savePlot=false;
parameters.caseName='biharmonic'; % supported cases:biharmonic,exPicard,imPicard,newton 
%---------------------------------------------------

% read command line args
for i=1:nargin
    line = varargin{i};
    if(strcmp(line,'-savePlot'))
        parameters.savePlot=true;
    elseif(strncmp(line,'-case=',6))
        parameters.caseName=line(7:end);
    end
    
%     if(strncmp(line,'-beam=',6))
%         beamOption = line(7:end);
%     end
%      if(strncmp(line,'-tf=',3))
%         tf = sscanf(varargin{i},'-tf=%e');
%      end   
end



solve(parameters);



fprintf('%sFinished running caseName=%s successfully.\n',infoPrefix,parameters.caseName);


end