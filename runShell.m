function runShell(varargin)
%=========================================================================
% This is the main interface for setting up parameters and running various
% cases for the shell paper
% usage:
% runShell -options
%=========================================================================
infoPrefix = '--runShell--: '; % all info displayed by this function includes this prefix

%--------------- some default options -----------------
parameters.isPlot=true;
parameters.savePlot=false;
parameters.useLU=false;
parameters.caseName='biharmonic'; % supported cases:biharmonic,exPicard,imPicard,newton 
parameters.bcType=2; % bcTypes:0 periodic; 1 simply supported; 2 clamped edge; 3 free edge; 4 CS; 5 CF
% domain=[xa,xb]x[ya,yb]
parameters.xa=0.;
parameters.xb=1.;
parameters.ya=0.;
parameters.yb=1.;
% resolution nx by ny
parameters.nx=10;
parameters.ny=10;
% Name of the m file that defines the rhs.
% we use function handle to define rhs: rhs.w(x,y,t),rhs.p(x,y,t) for the W
% and phi equations
parameters.rhsFile='rhsFileDefault'; 
% physical parameters
parameters.D=1.;
parameters.nu=1.;
parameters.E=1.;
parameters.h=1.;
%---------------------------------------------------

% read command line args
for i=1:nargin
    line = varargin{i};
    if(strcmp(line,'-savePlot'))
        parameters.savePlot=true;
    elseif(strcmp(line,'-noplot'))
        parameters.isPlot=false;
    elseif(strcmp(line,'-useLU'))
        parameters.useLU=true;          
    elseif(strncmp(line,'-case=',6))
        parameters.caseName=line(7:end);
    elseif(strncmp(line,'-bcType=',8))
        parameters.bcType=sscanf(line,'-bcType=%i');
    elseif(strncmp(line,'-xa=',4))
        parameters.xa=sscanf(line,'-xa=%e');
    elseif(strncmp(line,'-xb=',4))
        parameters.xb=sscanf(line,'-xb=%e');
    elseif(strncmp(line,'-ya=',4))
        parameters.ya=sscanf(line,'-ya=%e');
    elseif(strncmp(line,'-yb=',4))
        parameters.yb=sscanf(line,'-yb=%e'); 
    elseif(strncmp(line,'-nx=',4))
        parameters.nx=sscanf(line,'-nx=%i');  
    elseif(strncmp(line,'-ny=',4))
        parameters.ny=sscanf(line,'-ny=%i');
    elseif(strncmp(line,'-rhsFile=',9))
        parameters.rhsFile=sscanf(line,'-rhsFile=%s');
    elseif(strncmp(line,'-D=',3))
        parameters.D=sscanf(line,'-D=%e');      
    elseif(strncmp(line,'-nu=',4))
        parameters.nu=sscanf(line,'-nu=%e'); 
    elseif(strncmp(line,'-E=',3))
        parameters.E=sscanf(line,'-E=%e'); 
    elseif(strncmp(line,'-h=',3))
        parameters.h=sscanf(line,'-h=%e');      
    end  
end


profile on
solve(parameters);
pf = profile('info');
save('profile.mat','pf');


fprintf('%sFinished running caseName=%s successfully.\n',infoPrefix,parameters.caseName);


end