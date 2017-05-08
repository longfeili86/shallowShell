%this script plots results from batchRun.pl
%and move all eps files to a folder named figure

% create figure directory if doesn't exist
destination='figure';
if ~exist(destination, 'dir')
  mkdir(destination);
  fprintf('%mkdir %s\n',destination);
end

% plot convRate start from grid 2. grid 1 is too coarse
convRate -case=biharmonic -test=biharmPolyTest -rStart=2 -rEnd=6;
convRate -case=biharmonic -test=biharmTrigTest -rStart=2 -rEnd=6;
convRate -case=coupledSystem -solver=exPicard -test=coupledSystemTest -rStart=2 -rEnd=6 -nonlinear;
convRate -case=coupledSystem -solver=imPicard -test=coupledSystemTest -rStart=2 -rEnd=6 -nonlinear;
convRate -case=coupledSystem -solver=fsolve -test=coupledSystemTest -rStart=2 -rEnd=6 -nonlinear;
convRate -case=coupledSystem -solver=newton -test=coupledSystemTest -rStart=2 -rEnd=6 -nonlinear;
convRate -case=coupledSystem -solver=exPicard -test=coupledSystemTest -rStart=2 -rEnd=6;
convRate -case=coupledSystem -solver=imPicard -test=coupledSystemTest -rStart=2 -rEnd=6;
convRate -case=coupledSystem -solver=fsolve -test=coupledSystemTest -rStart=2 -rEnd=6;
convRate -case=coupledSystem -solver=newton -test=coupledSystemTest -rStart=2 -rEnd=6;

close all;
movefile('*.eps',destination);


%plot solution and errors from the finest grid

%  coupled system
Cases={
       'biharmTrigTestSupportedG64',...
       'biharmTrigTestClampedG64',...
       'biharmTrigTestFreeG64',...
       'biharmPolyTestSupportedG64',...
       'biharmPolyTestClampedG64',...
       'biharmPolyTestFreeG64',...
       'exPicardLinearcoupledSystemTestSupportedG64',...
       'exPicardLinearcoupledSystemTestClampedG64',...
       'imPicardLinearcoupledSystemTestSupportedG64',...
       'imPicardLinearcoupledSystemTestClampedG64',...
       'fsolveLinearcoupledSystemTestSupportedG64',...
       'fsolveLinearcoupledSystemTestClampedG64',...
       'newtonLinearcoupledSystemTestSupportedG64',...
       'newtonLinearcoupledSystemTestClampedG64',...
       'exPicardNonlinearcoupledSystemTestSupportedG64',...
       'exPicardNonlinearcoupledSystemTestClampedG64',...
       'imPicardNonlinearcoupledSystemTestSupportedG64',...
       'imPicardNonlinearcoupledSystemTestClampedG64',...
       'fsolveNonlinearcoupledSystemTestSupportedG64',...
       'fsolveNonlinearcoupledSystemTestClampedG64',...
       'newtonNonlinearcoupledSystemTestSupportedG64',...
       'newtonNonlinearcoupledSystemTestClampedG64'...
       };
   
   
for i=1:length(Cases)
    source=Cases{i};
    arg1=sprintf('-f=%s',source);
    
    plotSavedResults(arg1);
    plotSavedResults(arg1,'-contour');
    close all;
    moveFigureInSavedResults(source,destination)
end

   
%% difference between solutions of different solvers

contour=true;
if(contour)
    suffix='Contour';
end
resultsDir='.'; % pwd

%clamped 
load exPicardNonlinearcoupledSystemTestClampedG64/results.mat
WplotExPicard=Wplot;
PHIplotExPicard=PHIplot;

load imPicardNonlinearcoupledSystemTestClampedG64/results.mat
WplotImPicard=Wplot;
PHIplotImPicard=PHIplot;

load newtonNonlinearcoupledSystemTestClampedG64/results.mat
WplotNewton=Wplot;
PHIplotNewton=PHIplot;

eWNE = WplotNewton-WplotExPicard;
ePHINE = PHIplotNewton-PHIplotExPicard;

eWNI = WplotNewton-WplotImPicard;
ePHINI = PHIplotNewton-PHIplotImPicard;

figure
mySurf(Xplot,Yplot,eWNE,'$w$ difference (Newton-ExPicard)',contour);
printPlot(sprintf('differenceWNewtonExPicardClampedG64%s',suffix),resultsDir);

mySurf(Xplot,Yplot,ePHINE,'$\phi$ difference (Newton-ExPicard)',contour);
printPlot(sprintf('differencePHINewtonExPicardClampedG64%s',suffix),resultsDir);

mySurf(Xplot,Yplot,eWNI,'$w$ difference (Newton-ImPicard)',contour);
printPlot(sprintf('differenceWNewtonImPicardClampedG64%s',suffix),resultsDir);

mySurf(Xplot,Yplot,ePHINI,'$\phi$ difference (Newton-ImPicard)',contour);
printPlot(sprintf('differencePHINewtonImPicardClampedG64%s',suffix),resultsDir);


close all;
movefile('*.eps',destination);










