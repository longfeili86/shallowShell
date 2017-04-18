%this script does convergence study

grids=2:6; % do conv study on these grids

% simply supported
for i=1:length(grids);
    gn=2^(grids(i)-1);
    numGrids=gn*10;
    resultsDir=sprintf('SimplySupportedG%i',gn);
    cmd=sprintf('runShell -nx=%i -ny=%i -bcType=1 -noplot -f=%s',numGrids,numGrids,resultsDir);
    fprintf('--running--: %s\n',cmd);
    eval(cmd);
    load(sprintf('%s/results.mat',resultsDir));
    bc1ErrorL2(i)=norm(errPlot(:),2);
    bc1ErrorLmax(i)=max(abs(errPlot(:)));
    hh(i)=1./numGrids;
end

% clamped supported
for i=1:length(grids);
    gn=2^(grids(i)-1);
    numGrids=gn*10;
    resultsDir=sprintf('ClampedG%i',gn);
    cmd=sprintf('runShell -nx=%i -ny=%i -bcType=2 -noplot -f=%s',numGrids,numGrids,resultsDir);
    fprintf('--running--: %s\n',cmd);
    eval(cmd);
    load(sprintf('%s/results.mat',resultsDir));
    bc2ErrorL2(i)=norm(errPlot(:),2);
    bc2ErrorLmax(i)=max(abs(errPlot(:)));
    hh(i)=1./numGrids;
end

% free supported
for i=1:length(grids);
    gn=2^(grids(i)-1);
    numGrids=gn*10;
    resultsDir=sprintf('FreeG%i',gn);
    cmd=sprintf('runShell -nx=%i -ny=%i -bcType=3 -noplot -f=%s',numGrids,numGrids,resultsDir);
    fprintf('--running--: %s\n',cmd);
    eval(cmd);
    load(sprintf('%s/results.mat',resultsDir));
    bc3ErrorL2(i)=norm(errPlot(:),2);
    bc3ErrorLmax(i)=max(abs(errPlot(:)));
    hh(i)=1./numGrids;
end


setupFigure;
loglog(hh,(hh).^2,'k','LineWidth',figOptions.LW,'MarkerSize',figOptions.MS);
hold on
    loglog(hh,bc1ErrorLmax,'r-+','LineWidth',figOptions.LW,'MarkerSize',figOptions.MS); 
    loglog(hh,bc2ErrorLmax,'g-+','LineWidth',figOptions.LW,'MarkerSize',figOptions.MS); 
    loglog(hh,bc3ErrorLmax,'b-+','LineWidth',figOptions.LW,'MarkerSize',figOptions.MS); 
hold off

legend({'2nd order',....
   'Supported'...
   'Clamped'...
   'Free'},....
   'Location','NorthWest','FontSize',figOptions.FS);

xlabel('grid size','FontSize',figOptions.FS) ;
ylabel('error','FontSize',figOptions.FS) ;
set(gca,'FontSize',figOptions.FS);
box on
grid on
printPlot('.','convRateBiharmonicSolver.eps');




