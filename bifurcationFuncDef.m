% this function sets up the given functions for the bifurcation application


xi=-7100;
% forcing terms:
f.phi=@(x,y) 0.*x+0.*y+xi; % this is the thermal loading
f.w=@(x,y) 0.*x+0.*y;

% initial shell shape:
w0=@(x,y) 0.3*(1.-(x - 0.5).^2 - (y -0.5).^2);

% initial guess:
% wi=@(x,y) -w0(x,y).*100;