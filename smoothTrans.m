function omega=smoothTrans(x,xc,rc)
% this function build a smooth transition function from clamped region
% to other regions. The clamped region is defined by center xc and radius
% rc

e=0.01;
omega = 1-0.5*(tanh((abs(x-xc)-rc)/e)+1);

end