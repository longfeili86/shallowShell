
% exactOption= 1 for trig, 2 for poly
exactOption=2; 
isPlot=true;
contour=~true;
savePlot=~true;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(exactOption==1)
    fprintf('trig exact solution\n');
    forcingFigName='exactForcingtrig';
    exactFigName='exactSolutionTrig';  
elseif(exactOption==2)
    fprintf('poly exact solution\n');
    forcingFigName='exactForcingpoly';  
    exactFigName='exactSolutionPoly';  
else
   fprintf('exactOption = 1 for trig; 2 for poly. No other values are valid\n');
   return
end

if(contour)
   forcingFigName=sprintf('%sContour',forcingFigName);
   exactFigName=sprintf('%sContour',exactFigName); 
end

close all;

syms x y xa ya xb yb;

fx=2*pi/(xb-xa);
fy=2*pi/(yb-ya);

xc=(xa+xb)/2.;
yc=(ya+yb)/2.;
lx=(xb-xa)/3.;
ly=(yb-ya)/3.;

% the exact solution used for convergence test
Utrig = (sin(fx*(x-xa)).*sin(fy*(y-ya))).^4;
Upoly = ((x-xc).*(x-xa).*(x-xb).*(y-yc).*(y-ya).*(y-yb)/(lx^3*ly^3)).^7/100;
if(exactOption==1)
    UU=Utrig;
elseif(exactOption==2)
    UU=Upoly;
else
    fprintf('error\n');
    return;
end
% symbolic computation
Biharm_UU = diff(UU,x,4)+diff(UU,y,4)+2*diff(diff(UU,x,2),y,2);
UU_x = diff(UU,x,1);
UU_y = diff(UU,y,1);
UU_xy = diff(diff(UU,y,1),x,1);
UU_xx = diff(UU,x,2);
UU_yy = diff(UU,y,2);
UU_xxx = diff(UU,x,3);
UU_yyy = diff(UU,y,3);
UU_xxy = diff(UU_xx,y,1);
UU_yyx = diff(UU_yy,x,1);

Ue = matlabFunction(UU,'vars',[x,y,xa,ya,xb,yb]);
Ue_x = matlabFunction(UU_x,'vars',[x,y,xa,ya,xb,yb]);
Ue_y = matlabFunction(UU_y,'vars',[x,y,xa,ya,xb,yb]);
Ue_xy = matlabFunction(UU_y,'vars',[x,y,xa,ya,xb,yb]);
Ue_xx = matlabFunction(UU_xx,'vars',[x,y,xa,ya,xb,yb]);
Ue_yy = matlabFunction(UU_yy,'vars',[x,y,xa,ya,xb,yb]);
Ue_xxx = matlabFunction(UU_xxx,'vars',[x,y,xa,ya,xb,yb]);
Ue_yyy = matlabFunction(UU_yyy,'vars',[x,y,xa,ya,xb,yb]);
Ue_xxy = matlabFunction(UU_xxy,'vars',[x,y,xa,ya,xb,yb]);
Ue_yyx = matlabFunction(UU_yyx,'vars',[x,y,xa,ya,xb,yb]);
Biharm_Ue = matlabFunction(Biharm_UU,'vars',[x,y,xa,ya,xb,yb]);

clear x y xa ya xb yb;

% define grid
xa=0;xb=1;ya=0;yb=1;
domain=[xa,xb,ya,yb];
n=10*32;
nx=n;ny=n;
myGrid = buildGrid(domain,nx,ny);
Xvec = myGrid.XX(:);%column vector
Yvec = myGrid.YY(:);%column vector
hx=myGrid.hx;
hy=myGrid.hy;

% define index
Index=getIndex(nx,ny);
Xplot = reshape(Xvec(Index.interiorBoundary),ny,nx);
Yplot = reshape(Yvec(Index.interiorBoundary),ny,nx);


rhsf=@(x,y) Biharm_Ue(x,y,xa,ya,xb,yb);
 fprintf('check compatibility conditions: integral(rhs)=%e\n',integral2(rhsf,xa,xb,ya,yb));

if(isPlot)
figure
Uplot=Biharm_Ue(Xplot,Yplot,xa,ya,xb,yb);
mySurf(Xplot,Yplot,Uplot,'Forcing',contour);
if(savePlot)
    printPlot(forcingFigName);
end

figure
Uplot=Ue(Xplot,Yplot,xa,ya,xb,yb);
mySurf(Xplot,Yplot,Uplot,'Exact Solution',contour);
if(savePlot)
    printPlot(exactFigName);
end

figure
Uplot=Ue_x(Xplot,Yplot,xa,ya,xb,yb);
mySurf(Xplot,Yplot,Uplot,'Ue_x',contour);

figure
Uplot=Ue_y(Xplot,Yplot,xa,ya,xb,yb);
mySurf(Xplot,Yplot,Uplot,'Ue_y',contour);

figure
Uplot=Ue_xy(Xplot,Yplot,xa,ya,xb,yb);
mySurf(Xplot,Yplot,Uplot,'Ue_{xy}',contour);

figure
Uplot=Ue_xx(Xplot,Yplot,xa,ya,xb,yb);
mySurf(Xplot,Yplot,Uplot,'Ue_{xx}',contour);

figure
Uplot=Ue_yy(Xplot,Yplot,xa,ya,xb,yb);
mySurf(Xplot,Yplot,Uplot,'Ue_{yy}',contour);

figure
Uplot=Ue_xxx(Xplot,Yplot,xa,ya,xb,yb);
mySurf(Xplot,Yplot,Uplot,'Ue_{xxx}',contour);

figure
Uplot=Ue_yyy(Xplot,Yplot,xa,ya,xb,yb);
mySurf(Xplot,Yplot,Uplot,'Ue_{yyy}',contour);    
    

figure
Uplot=Ue_xxy(Xplot,Yplot,xa,ya,xb,yb);
mySurf(Xplot,Yplot,Uplot,'Ue_{xxy}',contour); 

figure
Uplot=Ue_yyx(Xplot,Yplot,xa,ya,xb,yb);
mySurf(Xplot,Yplot,Uplot,'Ue_{yyx}',contour); 
end


  