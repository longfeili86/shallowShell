function biharmonic(parameters)
% this function do all the tests for biharmonic solver
% We solve:
%   pde: \nabla^4 w = rhs
%   bcTypes: 1 simply supported; 2 clamped edge; 3 free edge
%         1: w=0, d^2wdn^2=0
%         2: w=0, dwdn=0
%         3: d^2wdn^2+nu*d^2ds^2=0, d^3wdn^3+(2-nu)d^3wdnds^2=0
%         4: mixed CS to do FINISH ME ......
%         5: mixed CF to do FINISH ME ......

infoPrefix = '--biharmonic--: '; % all info displayed by this function includes this prefix


% parse parameters
domain=parameters.domain; % rectangle for now [xa,xb,ya,yb]
nu=parameters.nu; % physical parameters (poisson ratio)
bcType=parameters.bcType;% bcType: 2 clamped edge, 3 free edge
resolution=parameters.resolution; % resolution of the grid [nx,ny]
nx=resolution(1);ny=resolution(2);
rhs=parameters.rhs; % this rhs is defined as a function handle rhs(x,y)

isPlot=parameters.isPlot;
savePlot=parameters.savePlot;

% preprossess
myGrid = buildGrid(domain,nx,ny);
Xvec = myGrid.XX(:);%column vector
Yvec = myGrid.YY(:);%column vector
RHS = rhs(Xvec,Yvec);%vectorized rhs

hx = myGrid.hx;
hy = myGrid.hy;
mtx = getDiffMatrix(nx,ny,hx,hy);
A = mtx.BiDh;

Index=getIndex(nx,ny);
A = assignBoundaryConditionsCoefficient(A,Index,mtx,parameter);
RHS = assignBoundaryConditionsRHS(RHS,Index,parameter);

% we need a cornor condition to remove the singulariy for free bc
if(bcType==2) 
    [A,RHS]=assignCornerConditions(A,RHS,Index);
end

% solve

% We might need condition number, eig values ect. of the matrix.
% There are unused ghost points, we need to remove those points from the A.
Aused=A(Index.UsedPoints,Index.UsedPoints);

% we use LU decomposition of Aused. This doesn't do much for a single
% biharmonic solve, but for iterative methods, we can reuse L, U to 
% to solve the system faster than simply backslash Aused.
[L,U]=lu(Aused);
y = L\RHS;
x = U\y;
W = 0.*RHS; % zero out stuff to store solution
W(Index.UsedPoints) = x;


% postprossess results

setupFigure; % setup figure options, linewidth,fontsize ect.
useColormapRainbow; % use rainbow colormap

%FontSize
if (isPlot)
    figure
    Xplot = reshape(Xvec(Index.interiorBoundary),nx,ny);
    Yplot = reshape(Yvec(Index.interiorBoundary),nx,ny);
    surf(Xplot,Yplot,reshape(U(Index.interiorBoundary),nx,ny));
    title('Solution','FontSize',figOptions.FS)
    set(gca,'FontSize',figOptions.FS)
    shading(figOptions.SD);
end



end