% the default file that defines the rhs for the W and phi eqns
% we can use other files to define rhs, simply pass the filename to
% using runShell -rhsFile='Filename'

rhs.w=@(x,y,t) 0.*x+0.*y-10.;
rhs.p=@(x,y,t) 0.*x+0.*y+0.;