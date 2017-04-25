% This is the default script that defines given functions:
% such as the external forcing terms of the W and phi equations: f.phi(x,y), f.w(x,y)
% and the exact solution (if known): exact.phi(x,y), exact.w(x,y)
% and the initial shell shape: w0(x,y).
% 
% we can write other .m script to define these function. To use a different
% funcDefFile, simply pass the file name to the runShell -funcDefFile=<filename>



% forcing term
f.w=@(x,y) 32634.2*max(-100.0*((x - 0.75)*2 + (y - 0.25)*2) + 1,0);
f.phi=@(x,y) 0;



%exact solution is known for this test
%knownExactSolution=true;
%exact.w=@(x,y) sin(fx*(x-xa)).^4.*sin(fy*(y-ya)).^4;
%exact.phi=@(x,y) 0.*x+0.*y;

% initial shell shape: w0(x,y).
w0=@(x,y) 0.1 - 0.4*(y - 0.5)*2;