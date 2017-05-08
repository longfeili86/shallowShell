
isPlot=true;
contour=true;
savePlot=~true;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% if(contour)
%    forcingFigName=sprintf('%sContour',forcingFigName);
%    exactFigName=sprintf('%sContour',exactFigName); 
% end

close all;

syms x y xa ya xb yb lambda;

fx=2*pi/(xb-xa);
fy=2*pi/(yb-ya);


% the exact solution used for convergence test
w = (sin(fx*(x-xa)).*sin(fy*(y-ya))).^4;
phi =  (sin(fx*(x-xa)).*sin(fy*(y-ya))).^5;
w0= (sin(fx*(x-xa)).*sin(fy*(y-ya)));

L=@(u,v) diff(u,x,2).*diff(v,y,2)+ diff(u,y,2).*diff(v,x,2)-2.*diff(diff(u,y,1),x,1).*diff(diff(v,y,1),x,1);
% symbolic computation
Biharm_w = diff(w,x,4)+diff(w,y,4)+2*diff(diff(w,x,2),y,2);
w_x = diff(w,x,1);
w_y = diff(w,y,1);
w_xy = diff(diff(w,y,1),x,1);
w_xx = diff(w,x,2);
w_yy = diff(w,y,2);
w_xxx = diff(w,x,3);
w_yyy = diff(w,y,3);
w_xxy = diff(w_xx,y,1);
w_yyx = diff(w_yy,x,1);
w_f =Biharm_w- L(w,phi)*lambda-L(w0,phi);


Biharm_phi = diff(phi,x,4)+diff(phi,y,4)+2*diff(diff(phi,x,2),y,2);
phi_x = diff(phi,x,1);
phi_y = diff(phi,y,1);
phi_xy = diff(diff(phi,y,1),x,1);
phi_xx = diff(phi,x,2);
phi_yy = diff(phi,y,2);
phi_xxx = diff(phi,x,3);
phi_yyy = diff(phi,y,3);
phi_xxy = diff(phi_xx,y,1);
phi_yyx = diff(phi_yy,x,1);
phi_f =  -(Biharm_phi+ 0.5*L(w,w)*lambda+L(w0,w));


we = matlabFunction(w,'vars',[x,y,xa,ya,xb,yb]);
we_x = matlabFunction(w_x,'vars',[x,y,xa,ya,xb,yb]);
we_y = matlabFunction(w_y,'vars',[x,y,xa,ya,xb,yb]);
we_xy = matlabFunction(w_y,'vars',[x,y,xa,ya,xb,yb]);
we_xx = matlabFunction(w_xx,'vars',[x,y,xa,ya,xb,yb]);
we_yy = matlabFunction(w_yy,'vars',[x,y,xa,ya,xb,yb]);
we_xxx = matlabFunction(w_xxx,'vars',[x,y,xa,ya,xb,yb]);
we_yyy = matlabFunction(w_yyy,'vars',[x,y,xa,ya,xb,yb]);
we_xxy = matlabFunction(w_xxy,'vars',[x,y,xa,ya,xb,yb]);
we_yyx = matlabFunction(w_yyx,'vars',[x,y,xa,ya,xb,yb]);
Biharm_we = matlabFunction(Biharm_w,'vars',[x,y,xa,ya,xb,yb]);
f_we =  matlabFunction(w_f,'vars',[x,y,xa,ya,xb,yb,lambda]);

phie = matlabFunction(phi,'vars',[x,y,xa,ya,xb,yb]);
phie_x = matlabFunction(phi_x,'vars',[x,y,xa,ya,xb,yb]);
phie_y = matlabFunction(phi_y,'vars',[x,y,xa,ya,xb,yb]);
phie_xy = matlabFunction(phi_y,'vars',[x,y,xa,ya,xb,yb]);
phie_xx = matlabFunction(phi_xx,'vars',[x,y,xa,ya,xb,yb]);
phie_yy = matlabFunction(phi_yy,'vars',[x,y,xa,ya,xb,yb]);
phie_xxx = matlabFunction(phi_xxx,'vars',[x,y,xa,ya,xb,yb]);
phie_yyy = matlabFunction(phi_yyy,'vars',[x,y,xa,ya,xb,yb]);
phie_xxy = matlabFunction(phi_xxy,'vars',[x,y,xa,ya,xb,yb]);
phie_yyx = matlabFunction(phi_yyx,'vars',[x,y,xa,ya,xb,yb]);
Biharm_phie = matlabFunction(Biharm_phi,'vars',[x,y,xa,ya,xb,yb]);
f_phie =  matlabFunction(phi_f,'vars',[x,y,xa,ya,xb,yb,lambda]);




clear x y xa ya xb yb lambda;

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

if(isPlot)
figure
Uplot=we(Xplot,Yplot,xa,ya,xb,yb);
mySurf(Xplot,Yplot,Uplot,'$w_e$',contour);

figure
Uplot=we_x(Xplot,Yplot,xa,ya,xb,yb);
mySurf(Xplot,Yplot,Uplot,'${w_e}_{x}$',contour);

figure
Uplot=we_y(Xplot,Yplot,xa,ya,xb,yb);
mySurf(Xplot,Yplot,Uplot,'${w_e}_{y}$',contour);

figure
Uplot=we_xy(Xplot,Yplot,xa,ya,xb,yb);
mySurf(Xplot,Yplot,Uplot,'${w_e}_{xy}$',contour);

figure
Uplot=we_xx(Xplot,Yplot,xa,ya,xb,yb);
mySurf(Xplot,Yplot,Uplot,'${w_e}_{xx}$',contour);

figure
Uplot=we_yy(Xplot,Yplot,xa,ya,xb,yb);
mySurf(Xplot,Yplot,Uplot,'${w_e}_{yy}$',contour);

figure
Uplot=we_xxx(Xplot,Yplot,xa,ya,xb,yb);
mySurf(Xplot,Yplot,Uplot,'${w_e}_{xxx}$',contour);

figure
Uplot=we_yyy(Xplot,Yplot,xa,ya,xb,yb);
mySurf(Xplot,Yplot,Uplot,'${w_e}_{yyy}$',contour);

figure
Uplot=we_xxy(Xplot,Yplot,xa,ya,xb,yb);
mySurf(Xplot,Yplot,Uplot,'${w_e}_{xxy}$',contour);

figure
Uplot=we_yyx(Xplot,Yplot,xa,ya,xb,yb);
mySurf(Xplot,Yplot,Uplot,'${w_e}_{yyx}$',contour);

figure
Uplot=Biharm_we(Xplot,Yplot,xa,ya,xb,yb);
mySurf(Xplot,Yplot,Uplot,'$\nabla^4w_e$',contour);

figure
Uplot=f_we(Xplot,Yplot,xa,ya,xb,yb,0);
mySurf(Xplot,Yplot,Uplot,'$f_w$ (Linear case)',contour);

figure
Uplot=f_we(Xplot,Yplot,xa,ya,xb,yb,1);
mySurf(Xplot,Yplot,Uplot,'$f_w$ (Nonlinear case)',contour);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% phi
figure
Uplot=phie(Xplot,Yplot,xa,ya,xb,yb);
mySurf(Xplot,Yplot,Uplot,'$\phi_e$',contour);

figure
Uplot=phie_x(Xplot,Yplot,xa,ya,xb,yb);
mySurf(Xplot,Yplot,Uplot,'${\phi_e}_{x}$',contour);

figure
Uplot=phie_y(Xplot,Yplot,xa,ya,xb,yb);
mySurf(Xplot,Yplot,Uplot,'${\phi_e}_{y}$',contour);

figure
Uplot=phie_xy(Xplot,Yplot,xa,ya,xb,yb);
mySurf(Xplot,Yplot,Uplot,'${\phi_e}_{xy}$',contour);

figure
Uplot=phie_xx(Xplot,Yplot,xa,ya,xb,yb);
mySurf(Xplot,Yplot,Uplot,'${\phi_e}_{xx}$',contour);

figure
Uplot=phie_yy(Xplot,Yplot,xa,ya,xb,yb);
mySurf(Xplot,Yplot,Uplot,'${\phi_e}_{yy}$',contour);

figure
Uplot=phie_xxx(Xplot,Yplot,xa,ya,xb,yb);
mySurf(Xplot,Yplot,Uplot,'${\phi_e}_{xxx}$',contour);

figure
Uplot=phie_yyy(Xplot,Yplot,xa,ya,xb,yb);
mySurf(Xplot,Yplot,Uplot,'${\phi_e}_{yyy}$',contour);

figure
Uplot=phie_xxy(Xplot,Yplot,xa,ya,xb,yb);
mySurf(Xplot,Yplot,Uplot,'${\phi_e}_{xxy}$',contour);

figure
Uplot=phie_yyx(Xplot,Yplot,xa,ya,xb,yb);
mySurf(Xplot,Yplot,Uplot,'${\phi_e}_{yyx}$',contour);

figure
Uplot=Biharm_phie(Xplot,Yplot,xa,ya,xb,yb);
mySurf(Xplot,Yplot,Uplot,'$\nabla^4\phi_e$',contour);

figure
Uplot=f_phie(Xplot,Yplot,xa,ya,xb,yb,0);
mySurf(Xplot,Yplot,Uplot,'$f_\phi$ (Linear case)',contour);

figure
Uplot=f_phie(Xplot,Yplot,xa,ya,xb,yb,1);
mySurf(Xplot,Yplot,Uplot,'$f_\phi$ (Nonlinear case)',contour);

end


  