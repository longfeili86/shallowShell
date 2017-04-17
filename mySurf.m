function mySurf(Xplot,Yplot,Zplot,titleMessage)
    
    infoPrefix = '--mySurf--: '; % all info displayed by this function includes this prefix

    useColormapRainbow; % use rainbow colormap
    setupFigure; % setup figure options, linewidth,fontsize ect.
    
    surf(Xplot,Yplot,Zplot);
    shading(figOptions.SD);
    if(figOptions.CT)
        hold on
        [~,hh]=contour3(Xplot,Yplot,Zplot,figOptions.NL);
        for i=1:numel(hh)
            set(hh(i),'EdgeColor','k','LineWidth',figOptions.LW)
        end
        hold off
    end   
    xlabel('x','FontSize',figOptions.FS);
    ylabel('y','FontSize',figOptions.FS)
    title(titleMessage,'FontSize',figOptions.FS);
    set(gca,'FontSize',figOptions.FS);
    view(figOptions.VW);
end