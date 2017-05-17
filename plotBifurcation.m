function plotBifurcation(varargin)

% this function plot bifurcation from the saved results in resultsName
% usage: plotBifurcation -f=<resultsName> (optional: -showIC -distBranch -color=<r> -savePlot)
% --Longfei Li
resultsName=''; 
savePlot=false;
showIC=false;
distBranch=false;
color='b';

if(nargin==0)
    fprintf('Usage: plotBifurcation -f=<resultsName> (optional: -showIC -distBranch -color=<r> -savePlot)\n');
    return
end

for i=1:nargin
    line = varargin{i};
    if(strncmp(line,'-f=',3))
        resultsName=line(4:end); 
    elseif(strcmp(line,'-savePlot'))
        savePlot=true;
    elseif(strcmp(line,'-showIC'))
        showIC=true;
    elseif(strcmp(line,'-distBranch'))
        distBranch=true;
    elseif(strncmp(line,'-color=',7))
        color=sscanf(line,'-color=%c');
    end
end


branchShapes='+>*';
branchNames='abc';
setupFigure;

hold on
for b=1:3 % branches
    wc=0;Xi=0;wi=0;
    Dirs= dir(sprintf('%s*%c',resultsName,branchNames(b)));   % number of dirs
    nDirs=length(Dirs);
    for i=1:nDirs
        results=sprintf('%s/results.mat',Dirs(i).name);
        load(results);
        [nx,ny]=size(Wplot);
        wc(i)=Wplot(nx/2,ny/2);
        wi(i)=Wiplot(nx/2,ny/2); % ic    
        Xi(i)=xi;
    end
    brachShape=branchShapes(1);
    if(distBranch)
        brachShape=branchShapes(b);
    end
    plot(Xi,wc,strcat(color,brachShape),'LineWidth',figOptions.LW);
    if(showIC)
        plot(Xi,wi,'k.','LineWidth',figOptions.LW);
    end
    set(gca,'FontSize',figOptions.FS);
end
hold off


% save eps file
if(savePlot)
    figName=sprintf('%s',resultsName);
    printPlot(figName);
end

end