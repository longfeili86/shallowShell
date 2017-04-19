% the default file that defines the rhs for the W and phi eqns
% we can use other files to define rhs, simply pass the filename to
% using runShell -rhsFile='Filename'
% exact solution if any is defined here as well

fx=2*pi/(xb-xa);
fy=2*pi/(yb-ya);
xc=(xa+xb)/2.;
yc=(ya+yb)/2.;

rhs.w=@(x,y) 256*fx^4*sin(fx*x - fx*xa).^4.*sin(fy*y - fy*ya).^4 - 240*fx^4*sin(fx*x - fx*xa).^2.*sin(fy*y - fy*ya).^4 + 24*fx^4*sin(fy*y - fy*ya).^4 + 512*fx^2*fy^2*sin(fx*x - fx*xa).^4.*sin(fy*y - fy*ya).^4 - 384*fx^2*fy^2*sin(fx*x - fx*xa).^4.*sin(fy*y - fy*ya).^2 - 384*fx^2*fy^2*sin(fx*x - fx*xa).^2.*sin(fy*y - fy*ya).^4 + 288*fx^2*fy^2*sin(fx*x - fx*xa).^2.*sin(fy*y - fy*ya).^2 + 256*fy^4*sin(fx*x - fx*xa).^4.*sin(fy*y - fy*ya).^4 - 240*fy^4*sin(fx*x - fx*xa).^4.*sin(fy*y - fy*ya).^2 + 24*fy^4*sin(fx*x - fx*xa).^4;
rhs.p=@(x,y) 0.*x+0.*y+0.;



%exact solution is known for this test
knownExactSolution=true;
exact=@(x,y) sin(fx*(x-xa)).^4.*sin(fy*(y-ya)).^4;
