function coeff = ClampedSimplysupportedCoefficients(dx)
% Compute finite difference coefficient for singularity where the boundary
% condition changes from clamped to free
% modify the coefficient for inhomogeneous biharmonic equation


r = [sqrt(2)*dx, sqrt(2)*dx,2*dx,sqrt(5)*dx,sqrt(5)*dx,sqrt(5)*dx,sqrt(5)*dx]';
theta = [3*pi/4,pi/4,pi/2,pi-atan(0.5),atan(0.5),pi-atan(2),atan(2)]';


F = zeros(7,7);


F(1,:) = r.^1.5.*(cos(1.5*theta)-cos(0.5*theta));
F(2,:) = r.^2.5.*(cos(2.5*theta)-cos(0.5*theta));
F(3,:) = r.^3.5.*(cos(3.5*theta)-cos(1.5*theta));
F(4,:) = r.^3.*(sin(3*theta)-3*sin(theta));
F(5,:) = r.^4.*(2.0*sin(4*theta)-4.0*sin(2.0*theta));
F(6,:) = r.^2.*log(r);
F(7,:) = r.^2;


r0 = dx;
theta0 = pi/2;
v0 = zeros(7,1);
v0(1) = r0^1.5*(cos(1.5*theta0)-cos(0.5*theta0));
v0(2) = r0^2.5*(cos(2.5*theta0)-cos(0.5*theta0));
v0(3) = r0^3.5*(cos(3.5*theta0)-cos(1.5*theta0));
v0(4) = r0^3*(sin(3*theta0)-3*sin(theta0));
v0(5) = r0^4*(2.0*sin(4*theta0)-4.0*sin(2.0*theta0));
v0(6) = r0^2*log(r0);
v0(7) = r0^2;


coeff = F\v0;
coeff = [-1;coeff];



end

