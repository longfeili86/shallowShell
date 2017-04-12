function runShell(varargin)
%=========================================================================
% This is the main interface for running various cases for the shell paper
% usage:
% runShell -options
%=========================================================================

%--------------- default options -----------------
savePlot=false;
caseName='biharmonic'; % supported cases:biharmonic,exPicard,imPicard,newton 
%---------------------------------------------------

% read command line args
for i=1:nargin
    line = varargin{i};
    if(strcmp(line,'-savePlot'))
        savePlot=true;
    elseif(strncmp(line,'-case=',6))
        caseName=line(7:end);
    end
    
%     if(strncmp(line,'-beam=',6))
%         beamOption = line(7:end);
%     end
%      if(strncmp(line,'-tf=',3))
%         tf = sscanf(varargin{i},'-tf=%e');
%      end   
end


infoPrefix = '--runShell--: '; % all info displayed by this function includes this prefix
switch (caseName)
    case 'biharmonic'
        fprintf('%sTest the biharmonic solver\n',infoPrefix);
    case 'exPicard'
        fprintf('%sTest the explicit Picard iterations\n',infoPrefix);
        fprintf('%sFinish me...\n',infoPrefix);
        return
    case 'exPicard'
        fprintf('%sTest the semi-implicit Picard iterations\n',infoPrefix);
        fprintf('%sFinish me...\n',infoPrefix);
        return
    case 'newton'
        fprintf('%sTest the newton iterations (direct solve)\n',infoPrefix);
        fprintf('%sFinish me...\n',infoPrefix);
        return
    otherwise
        fprintf('%sSupported cases: biharmonic,exPicard,imPicard,newton.\n',infoPrefix);
        fprintf('%sTerminated.\n',infoPrefix);
        return
end


figureName='tempFig.eps';
if(savePlot)
    fprintf('%splot saved. Filename=%s\n',infoPrefix,figureName);
    %print('-depsc2',figureName);
end






end