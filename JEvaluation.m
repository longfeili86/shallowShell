function J=JEvaluation(x,n,W0,Aphi,Aw,mtx,Index,parameters,dx)
% this function computes the jacobian matrix of F(x) at x
%  Fp = Aphi*phi+0.5*L[w,w]+L[w0,w]         phi equation
%  Fw = Aw*w-L[w,phi]-L[w0,phi]             w equation
%  
%  Jpp = Aphi;  
%  Jpw = 0.5 dLww+Lw0;
%  Jwp = -Lw-Lw0;
%  Jww = Aw-Lphi;
%
%  bc has already implemented in Aw and Aphi
%  we need to implement bc for L matrices
% 
%    |Jpp Jpw|
% J= |       |
%    |Jwp Jww|
%
% -- Longfei Li


Nphi=1:n; % phi solutions
Nw=n+1:2*n;  % w solutions
Nlambda=2*n+1:2*n+3; % lambda solutions

phi=x(Nphi);
w=x(Nw);

Lw0= getLMatrix(mtx,W0);
Lw0=assignBoundaryConditionsLMatrix(Lw0,Index,parameters); 
% first add all linear components:
Jpp = Aphi;

nadd=0; % number of additional varibles in w equations
if(parameters.bcType==3)
    nadd=3; % add 3 additional varibles to remove singularity
end
Jpw = sparse(n,n+nadd); % n phi eqns and n+nadd w variables
Jpw(1:n,1:n) = Lw0;
Jwp = sparse(n+nadd,n); % n+nadd w eqns and n phi variables
Jwp(1:n,1:n) = -Lw0;
Jww = Aw; 
% add nonlinear components
if(~parameters.isLinear)
    Lw = getLMatrix(mtx,w);
    Lw=assignBoundaryConditionsLMatrix(Lw,Index,parameters); 
    Lphi = getLMatrix(mtx,phi);
    Lphi=assignBoundaryConditionsLMatrix(Lphi,Index,parameters); 

    dLww = 2.*Lw;
    Jpw(1:n,1:n) = Jpw(1:n,1:n)+ 0.5*dLww; 
    Jwp(1:n,1:n) = Jwp(1:n,1:n)-Lw;
    Jww(1:n,1:n) = Jww(1:n,1:n)-Lphi;
end

J = [Jpp,Jpw;
     Jwp,Jww];
 
% augment J for PAC method 
if(parameters.usePAC) 
    Lpx = ones(n,1); % xi derivatives of phi eqns
    Lpx=assignBoundaryConditionsLMatrix(Lpx,Index,parameters); 
    J = [J,[Lpx;zeros(n+nadd,1)];
         dx'];
end

 

end