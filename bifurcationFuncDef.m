% this is the default funcDefFile for bifurcation runs
% it sets up the given functions for the bifurcation application


xi=parameters.xi;

fprintf('thermal loading xi=%f\n',xi);

% forcing terms:
f.phi=@(x,y) 0.*x+0.*y+xi; % this is the thermal loading
f.w=@(x,y) 0.*x+0.*y;

% initial shell shape:
w0=@(x,y) 0.3*(1.-(x - 0.5).^2 - (y -0.5).^2);

% initial guess:
if(~isReadIC)
    wi=@(x,y) 0.*w0(x,y);
end