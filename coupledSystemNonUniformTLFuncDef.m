% this function sets up the given functions for solution of nonuniform
% thermal loading and CF bc



% forcing terms:
f.phi=@(x,y) 32634.2*max(-100.0*((x - 0.75).^2 + (y - 0.25).^2 ) + 1.,0); % this is the thermal loading
f.w=@(x,y) 0.*x+0.*y;

% initial shell shape:
w0=@(x,y) 0.1-0.4*(y-0.5).^2;

% initial guess:
if(~isReadIC)
wi=@(x,y) w0(x,y);
end