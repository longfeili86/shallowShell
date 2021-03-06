function h=plotBifurcation(varargin)

% this function plot bifurcation from the saved results in resultsName
% usage: plotBifurcation -f=<resultsName> (optional: -showIC -distBranch -color=<r> -savePlot)
% output:
%   handle of this plot
% --Longfei Li
resultsName=''; 
savePlot=false;
showIC=false;
distBranch=false;
color='b';
option='center'; % center point or l2 norm %plot the bifucation for the center point or l2 norm

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
    elseif(strncmp(line,'-option=',8))
        option=sscanf(line,'-option=%s');
    end
end

if(~strcmp(option,'center') && ~strcmp(option,'l2'))
    option='center'; % use center for unknonw options.
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
        Xi(i)=xi; % values for xi
        if(strcmp(option,'center'))
            wc(i)=Wplot(nx/2,ny/2);
            wi(i)=Wiplot(nx/2,ny/2); % ic    
        else
            wc(i)=sqrt(sum(Wplot(:).^2)/(nx*ny)); % l2 norm
            wi(i)=sqrt(sum(Wiplot(:).^2)/(nx*ny)); % ic  l2 norm  
        end
        
    end
    brachShape=branchShapes(1);
    if(distBranch)
        brachShape=branchShapes(b);
    end
    h=plot(Xi,wc,strcat(color,brachShape),'LineWidth',figOptions.LW);
    if(showIC)
        plot(Xi,wi,'k.','LineWidth',figOptions.LW);
    end
    set(gca,'FontSize',figOptions.FS);
end
hold off

% save eps file
if strcmp(option,'l2')
    resultsName=strcat(resultsName,'_l2');
end
if(savePlot)
    figName=sprintf('%s',resultsName);
    printPlot(figName);
end

end