function biharmonic(parameters)
% given parameters, this function solves:
%   pde: D*s\nabla^4 w = rhs
%   bcTypes: 0 periodic; 1 simply supported; 2 clamped edge; 3 free edge; 4 CS; 5 CF
%         0: periodic bc  to do FINISH ME ......
%         1: w=0, d^2wdn^2=0
%         2: w=0, dwdn=0
%         3: d^2wdn^2+nu*d^2ds^2=0, d^3wdn^3+(2-nu)d^3wdnds^2=0
%         4: mixed CS to do FINISH ME ......
%         5: mixed CF to do FINISH ME ......

infoPrefix = '--biharmonic--: '; % all info displayed by this function includes this prefix


% parse parameters
domain=[parameters.xa,parameters.xb,parameters.ya,parameters.yb]; % rectangle [xa,xb,ya,yb]
nu=parameters.nu; % physical parameters (poisson ratio)
bcType=parameters.bcType;% bcType: 2 clamped edge, 3 free edge
nx=parameters.nx; % number of grid points in x direction
ny=parameters.ny; % number of grid points in y direction

isPlot=parameters.isPlot;
savePlot=parameters.savePlot;
useLU=parameters.useLU;

rhsFile=parameters.rhsFile;

% physical parameters
D=parameters.D;
nu=parameters.nu;
E=parameters.E;
h=parameters.h;

% preprossess
fprintf('%sGetting rhs definition from file: %s.m\n',infoPrefix,rhsFile);
run(rhsFile);

myGrid = buildGrid(domain,nx,ny);
Xvec = myGrid.XX(:);%column vector
Yvec = myGrid.YY(:);%column vector
RHS = rhs.w(Xvec,Yvec);%vectorized rhs

hx = myGrid.hx;
hy = myGrid.hy;
mtx = getDiffMatrix(nx,ny,hx,hy);
A = D*mtx.BiDh;

Index=getIndex(nx,ny);
A = assignBoundaryConditionsCoefficient(A,Index,mtx,parameters);
RHS = assignBoundaryConditionsRHS(RHS,Index,parameters);

% we need a cornor condition to remove the singulariy for free bc
if(bcType==3) 
    [A,RHS]=assignCornerConditions(A,RHS,Index);
end

% solve

% We might need condition number, eig values ect. of the matrix.
% There are unused ghost points, we need to remove those points from the A.
Aused=A(Index.UsedPoints,Index.UsedPoints);

% we can use LU decomposition of Aused. This is slower for a single
% biharmonic solve, but for iterative methods, we can reuse L, U to 
% to solve the biharmonic system faster than simply backslash Aused every time.
if(useLU)
    fprintf('%suse LU decomposition\n',infoPrefix);
    %tic;
    [L,U]=lu(Aused,0.);
    %toc;
    %tic;
    y = L\RHS(Index.UsedPoints);
    x = U\y;
    %toc;
else
    %tic;
    x=Aused\RHS(Index.UsedPoints);
    %toc;
end
W = 0.*RHS; % zero out stuff to store solution
W(Index.UsedPoints) = x;

% postprossess results


if (isPlot)
    figure
    useColormapRainbow; % use rainbow colormap
    setupFigure; % setup figure options, linewidth,fontsize ect.

    Xplot = reshape(Xvec(Index.interiorBoundary),nx,ny);
    Yplot = reshape(Yvec(Index.interiorBoundary),nx,ny);
    Wplot = reshape(W(Index.interiorBoundary),nx,ny);
    surf(Xplot,Yplot,Wplot);
    shading(figOptions.SD);
    if(figOptions.CT)
        hold on
        [~,hh]=contour3(Xplot,Yplot,Wplot);
        for i=1:numel(hh)
            set(hh(i),'EdgeColor','k','LineWidth',figOptions.LW)
        end
        hold off
    end   
    title('Solution','FontSize',figOptions.FS)
    set(gca,'FontSize',figOptions.FS)
end



end