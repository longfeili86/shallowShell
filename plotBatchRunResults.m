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
%convRate -case=coupledSystem -solver=fsolve -test=coupledSystemTest -rStart=2 -rEnd=6 -nonlinear;
%convRate -case=coupledSystem -solver=newton -test=coupledSystemTest -rStart=2 -rEnd=6 -nonlinear;
convRate -case=coupledSystem -solver=exPicard -test=coupledSystemTest -rStart=2 -rEnd=6;
convRate -case=coupledSystem -solver=imPicard -test=coupledSystemTest -rStart=2 -rEnd=6;
%convRate -case=coupledSystem -solver=fsolve -test=coupledSystemTest -rStart=2 -rEnd=6;
%convRate -case=coupledSystem -solver=newton -test=coupledSystemTest -rStart=2 -rEnd=6;

close all;
movefile('*.eps',destination);


%plot solution and errors from the finest grid

%  coupled system
Cases={
       'biharmTrigTestSupportedG64',...
       'biharmTrigTestClampedG64',...
       'biharmTrigTestFreeG64',...
       'biharmTrigTestCSG64',...
       'biharmTrigTestCFG64',...
       'biharmPolyTestSupportedG64',...
       'biharmPolyTestClampedG64',...
       'biharmPolyTestFreeG64',...
       'biharmPolyTestCSG64',...
       'biharmPolyTestCFG64',...
       'exPicardLinearcoupledSystemTestSupportedG64',...
       'exPicardLinearcoupledSystemTestClampedG64',...
       'exPicardLinearcoupledSystemTestCSG64',....
       'exPicardLinearcoupledSystemTestCFG64',....
       'imPicardLinearcoupledSystemTestSupportedG64',...
       'imPicardLinearcoupledSystemTestClampedG64',...
       'imPicardLinearcoupledSystemTestCSG64',...
       'imPicardLinearcoupledSystemTestCFG64',...
       'exPicardNonlinearcoupledSystemTestSupportedG64',...
       'exPicardNonlinearcoupledSystemTestClampedG64',...
       'exPicardNonlinearcoupledSystemTestCSG64',...
       'exPicardNonlinearcoupledSystemTestCFG64',...
       'imPicardNonlinearcoupledSystemTestSupportedG64',...
       'imPicardNonlinearcoupledSystemTestClampedG64',...
       'imPicardNonlinearcoupledSystemTestCSG64',...
       'imPicardNonlinearcoupledSystemTestCFG64',...
       };
   
   
for i=1:length(Cases)
    source=Cases{i};
    arg1=sprintf('-f=%s',source);
    
    plotSavedResults(arg1);
    plotSavedResults(arg1,'-contour');
    close all;
    moveFigureInSavedResults(source,destination)
end

   
% plot bifurcation results
setupFigure;
colors='rbkcm';
bcNames={'Supported','Clamped','Free','CS','CF'};
hold on
for b=1:length(bcNames)
    hh(b)=plotBifurcation(strcat('-f=bifurcationExImEx',bcNames{b}),strcat('-color=',colors(b)));
end
hold off
box on
legend(hh,bcNames,'FontSize',figOptions.FS,'Location','best');
printPlot('bifurcationExImEx');
close all
movefile('*.eps',destination);

%plot nonuniform thermal loading case:
solvers={'exPicard','imPicard','newton','fsolve'};
for k=1:length(solvers)
for i=1:length(bcNames)
  source=strcat('nonuniformTL',bcNames{i},'_',solvers{k});
  arg1=strcat('-f=',source);
  plotSavedResults(arg1);
  plotSavedResults(arg1,'-contour');
  close all
  moveFigureInSavedResults(source,destination);
end
end









