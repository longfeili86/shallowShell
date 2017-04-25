function coupledSystem(parameters)
% given parameters, this function solves the coupled system:
%   pde:  \nabla^4 phi = -1/2 L[w,w]*lambda - L[w0,w] -f.phi
%         \nabla^4 w   = L[w,phi]*lambda + L[w0,phi] + f.w
%   (Note: for linear case: lambda=0, for nonlinear case: lambda=1)
%   bcTypes: 0 periodic; 1 simply supported; 2 clamped edge; 3 free edge; 4 CS; 5 CF
%         0: periodic bc  to do FINISH ME ......
%         1: w=0, d^2wdn^2=0, phi=0, d^2phidn^2=0,
%         2: w=0, dwdn=0,  phi=0, dphidn=0    
%         3: d^2wdn^2+nu*d^2ds^2=0, d^3wdn^3+(2-nu)d^3wdnds^2=0, phi=0, dphidn=0 
%         4: mixed CS to do FINISH ME ......
%         5: mixed CF to do FINISH ME ......
% --Longfei Li

infoPrefix = '--coupledSystem--: '; % all info displayed by this function includes this prefix


% parse parameters
resultsDir=parameters.resultsDir;

xa=parameters.xa;xb=parameters.xb;ya=parameters.ya;yb=parameters.yb;
domain=[xa,xb,ya,yb]; % rectangle [xa,xb,ya,yb]
nu=parameters.nu; % physical parameters (poisson ratio)
bcType=parameters.bcType;% bcType: 2 clamped edge, 3 free edge
nx=parameters.nx; % number of grid points in x direction
ny=parameters.ny; % number of grid points in y direction

isPlot=parameters.isPlot;
savePlot=parameters.savePlot;
useLU=parameters.useLU;

isLinear=parameters.isLinear;
solver=parameters.solver;

funcDefFile=parameters.funcDefFile;
knownExactSolution=parameters.knownExactSolution;


% print some information
sysInfo='nonlinear';
if(isLinear)
    sysInfo='linear';
end
fprintf('%sThe coupled system is %s\n',infoPrefix,sysInfo);
fprintf('%sUsing solver: %s\n',infoPrefix,solver);


% preprocess

% define grid
myGrid = buildGrid(domain,nx,ny);
Xvec = myGrid.XX(:);%column vector
Yvec = myGrid.YY(:);%column vector

% define diff matrix
hx = myGrid.hx;
hy = myGrid.hy;
mtx = getDiffMatrix(nx,ny,hx,hy);

% define index
Index=getIndex(nx,ny);

% define given functions
fprintf('%sGetting definitions of all the given functions from file: %s.m\n',infoPrefix,funcDefFile);
run(funcDefFile); % f.w(x,y), f.phi(x,y) and w0(x,y) are defined here
F.phi=f.phi(Xvec,Yvec);
F.w=f.w(Xvec,Yvec);
W0=w0(Xvec,Yvec);

%initial guess is (W0,PHI0)
PHI0=0.*W0;

% setup equations
Aphi=mtx.BiDh;
Aw = mtx.BiDh;
%%%%%NOTE NOT WORKING FOR FREE BC FOR NOW!!! FINISH ME %%%%%%%%%%%%%%%
Aphi = assignBoundaryConditionsCoefficient(Aphi,Index,mtx,parameters);
Aw = assignBoundaryConditionsCoefficient(Aw,Index,mtx,parameters); 
%%%%%NOTE NOT WORKING FOR FREE BC FOR NOW!!! FINISH ME %%%%%%%%%%%%%%%

% bc for RHSs are already implemented in side of getRHS functions
% the RHSs are only for the Index.UsedPoints
RHSphi=@(w,phi) getRHS_phiEqn(w,phi,F.phi,W0,mtx,parameters,Index);
RHSw=@(w,phi)   getRHS_wEqn(w,phi,F.w,W0,mtx,parameters,Index);


% build problem for fsolve
nUsed=(nx+4)*(ny+4);
problem.options = optimoptions('fsolve','Display','iter','Algorithm', 'trust-region-dogleg');
problem.objective = @(x) [Aphi*x(1:nUsed)-RHSphi(x(1:nUsed),x(nUsed+1:end));Aw*x(nUsed+1:end)-RHSw(x(1:nUsed),x(nUsed+1:end))];
problem.x0 = [PHI0;W0];
problem.solver = 'fsolve';
[x,fval,exitflag,output] = fsolve(problem);

PHI =x(1:nUsed); 
W =  x(nUsed+1:end); 


% postprocess results
Xplot = reshape(Xvec(Index.interiorBoundary),ny,nx);
Yplot = reshape(Yvec(Index.interiorBoundary),ny,nx);
Wplot = reshape(W(Index.interiorBoundary),ny,nx);
PHIplot = reshape(PHI(Index.interiorBoundary),ny,nx);

figure
mySurf(Xplot,Yplot,Wplot,'w');

figure
mySurf(Xplot,Yplot,PHIplot,'\phi');


end