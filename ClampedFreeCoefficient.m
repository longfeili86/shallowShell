function coeff = ClampedFreeCoefficient(dx,nu)
% Compute finite difference coefficient for singularity where the boundary
% condition changes from clamped to free

%% 6 points stencil
k = atanh((1+nu)/2)/pi;
r = [sqrt(2)*dx, sqrt(2)*dx,2*dx,sqrt(5)*dx,sqrt(5)*dx]';
theta = [3*pi/4,pi/4,pi/2,pi-atan(0.5),atan(0.5)]';

F = zeros(5,5);

q = complex(0.5,k);
coeffB = -2*(q-1)/(1+nu)*cot(q*pi);
F(1,:) = r.^(q+1).*(cos((q+1)*theta)-cos((q-1)*theta)+...
    coeffB*(sin((q+1)*theta)/(q+1)-sin((q-1)*theta)/(q-1)));

q = complex(0.5,-k);
coeffB = -2*(q-1)/(1+nu)*cot(q*pi);
F(2,:) = r.^(q+1).*(cos((q+1)*theta)-cos((q-1)*theta)+...
    coeffB*(sin((q+1)*theta)/(q+1)-sin((q-1)*theta)/(q-1)));


F(3,:) = r.^2;
F(4,:) = r.^2.*log(r);
F(5,:) = r.^2.*theta;

r0 = dx;
theta0 = pi/2;
v0 = zeros(5,1);

q = complex(0.5,k);
coeffB = -2*(q-1)/(1+nu)*cot(q*pi);
v0(1) = r0^(q+1)*(cos((q+1)*theta0)-cos((q-1)*theta0)+coeffB*(sin((q+1)*theta0)/(q+1)-sin((q-1)*theta0)/(q-1)));

q = complex(0.5,-k);
coeffB = -2*(q-1)/(1+nu)*cot(q*pi);
v0(2) = r0^(q+1)*(cos((q+1)*theta0)-cos((q-1)*theta0)+coeffB*(sin((q+1)*theta0)/(q+1)-sin((q-1)*theta0)/(q-1)));

v0(3) = r0^2;
v0(4) = r0^2*log(r0);
v0(5) = r0^2*theta0;

coeff = F\v0;
coeff = [-1;real(coeff)];

end

