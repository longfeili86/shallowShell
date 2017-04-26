function runShell(varargin)
%=========================================================================
% This is the main interface for setting up parameters and running various
% cases for the shell paper
% usage:
% runShell -options
% --Longfei Li
%=========================================================================
infoPrefix = '--runShell--: '; % all info displayed by this function includes this prefix

%--------------- some default options -----------------
parameters.resultsDir='.'; % directory to save results. By default it's pwd
parameters.saveDiary=true; % flag to Save Command Window text to file
parameters.isPlot=true;
parameters.savePlot=false;
parameters.useLU=false;
parameters.caseName='biharmonic'; % supported cases:biharmonic,exPicard,imPicard,newton 
parameters.bcType=2; % bcTypes:0 periodic; 1 simply supported; 2 clamped edge; 3 free edge; 4 CS; 5 CF
% parameters for coupled system
parameters.isLinear=true; % flag indicate if the coupled system is linear or not
parameters.solver='fsolve'; % available solvers for the coupled system: 1. fsolve 2. exPicard 3. imPicard 
% domain=[xa,xb]x[ya,yb]
parameters.xa=0.;
parameters.xb=1.;
parameters.ya=0.;
parameters.yb=1.;
% resolution nx by ny
parameters.nx=10;
parameters.ny=10;
% Name of the m file that defines given functions:
% such as the external forcing terms of the W and phi equations: f.phi(x,y), f.w(x,y)
% and the exact solution (if known): exact.phi(x,y), exact.w(x,y)
% and the initial shell shape: w0(x,y).
parameters.funcDefFile='funcDefFileDefault'; 
% flag for exact solution. Exact solution (if any) is defined in rhsFile
% make sure the change this flag when defining an exact solution
parameters.knownExactSolution=false; 

%iteration parameters
parameters.maxIter=500;
parameters.tol=1e-6;

% physical parameters
parameters.nu=.1; % poisson ratio is ranging from 0.0 to 0.5
% these are not needed for dimensionless equations
parameters.E=1.;
parameters.h=1.;
parameters.D=1.; % D = Eh^3/(12(1-nu^2)), we specify it for now
%---------------------------------------------------

% read command line args
for i=1:nargin
    line = varargin{i};
    if(strncmp(line,'-f=',3))
        parameters.resultsDir=line(4:end); 
    elseif(strcmp(line,'-savePlot'))
        parameters.savePlot=true;
    elseif(strcmp(line,'-nodiary'))
        parameters.saveDiary=false;        
    elseif(strcmp(line,'-noplot'))
        parameters.isPlot=false;
    elseif(strcmp(line,'-useLU'))
        parameters.useLU=true;          
    elseif(strncmp(line,'-case=',6))
        parameters.caseName=line(7:end);
    elseif(strncmp(line,'-bcType=',8))
        parameters.bcType=sscanf(line,'-bcType=%i');
    elseif(strcmp(line,'-nonlinear'))
        parameters.isLinear=false;
    elseif(strncmp(line,'-solver=',8))
        parameters.solver=line(9:end);         
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
    elseif(strncmp(line,'-funcDefFile=',12))
        parameters.funcDefFile=sscanf(line,'-funcDefFile=%s');
    elseif(strncmp(line,'-maxIter=',9))
        parameters.maxIter=sscanf(line,'-maxIter=%i');
    elseif(strncmp(line,'-tol=',5))
        parameters.tol=sscanf(line,'-tol=%e'); 
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

% create resultsDir if doesn't exist
if ~exist(parameters.resultsDir, 'dir')
  mkdir(parameters.resultsDir);
  fprintf('%smkdir %s\n',infoPrefix,parameters.resultsDir);
end

% keep a diary for the run
if(parameters.saveDiary)
   diaryFilename=sprintf('%s/diary.txt',parameters.resultsDir);
   diary(diaryFilename); 
   diary on;
   fprintf('%sSave Command Window text to %s\n',infoPrefix,diaryFilename);
end

fprintf('%sRunning: caseName=%s\n',infoPrefix,parameters.caseName);
fprintf('%sResults are saved in directory: %s/\n',infoPrefix,parameters.resultsDir);

% save the profile for the solve
profile on %-history
solve(parameters);
profile off
pf = profile('info');
save(sprintf('%s/profile.mat',parameters.resultsDir),'pf');


fprintf('%sCode exited successfully.\n',infoPrefix);

if(parameters.saveDiary)
   diary off;
end
end