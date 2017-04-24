% This is the default script that defines given functions:
% such as the external forcing terms of the W and phi equations: f.phi(x,y), f.w(x,y)
% and the exact solution (if known): exact.phi(x,y), exact.w(x,y)
% and the initial shell shape: w0(x,y).
% 
% we can write other .m script to define these function. To use a different
% funcDefFile, simply pass the file name to the runShell -funcDefFile=<filename>


fx=2*pi/(xb-xa);
fy=2*pi/(yb-ya);
xc=(xa+xb)/2.;
yc=(ya+yb)/2.;

% forcing term
f.w=@(x,y) 256*fx^4*sin(fx*x - fx*xa).^4.*sin(fy*y - fy*ya).^4 - 240*fx^4*sin(fx*x - fx*xa).^2.*sin(fy*y - fy*ya).^4 + 24*fx^4*sin(fy*y - fy*ya).^4 + 512*fx^2*fy^2*sin(fx*x - fx*xa).^4.*sin(fy*y - fy*ya).^4 - 384*fx^2*fy^2*sin(fx*x - fx*xa).^4.*sin(fy*y - fy*ya).^2 - 384*fx^2*fy^2*sin(fx*x - fx*xa).^2.*sin(fy*y - fy*ya).^4 + 288*fx^2*fy^2*sin(fx*x - fx*xa).^2.*sin(fy*y - fy*ya).^2 + 256*fy^4*sin(fx*x - fx*xa).^4.*sin(fy*y - fy*ya).^4 - 240*fy^4*sin(fx*x - fx*xa).^4.*sin(fy*y - fy*ya).^2 + 24*fy^4*sin(fx*x - fx*xa).^4;
f.phi=@(x,y) 0.*x+0.*y+0.;



%exact solution is known for this test
knownExactSolution=true;
exact.w=@(x,y) sin(fx*(x-xa)).^4.*sin(fy*(y-ya)).^4;
exact.phi=@(x,y) 0.*x+0.*y;

% initial shell shape: w0(x,y).
w0=@(x,y) 0.*x+0.*y;