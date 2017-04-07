%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This is the main program to solve biharmonic equation
%  D \nabla^4 w(x,y) = p(x,y)
% 
% Clamped BC : w = w_n = 0
% or simply supported BC: w = w_nn = 0
% imposed on C
% 
% Free boundary conditions imposed on AllEdges\C
% Left & Right: w_xx+nu*w_yy=0, w_xxx+(2-nu)*w_xyy=0
% Top & Bottom: w_yy+nu*w_xx=0, w_yyy+(2-nu)*w_xxy=0
%
% corner conditions: w_{xy}=0
%
% Schematic boundary conditions:
% 
%  --------**********--------
% |                         |		  
% |                         |		  
% |                         |		  
% |                         |		  
%  --------**********--------
%
% by H. Ji 4/30/2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% use Richardson extrapolation result as exact solution
% second order accuracy with SS+Clamped, Free+Clamped BCs

close all
clear
clc

out = fopen('accuracy_biharmonic.dat','w');

% implement mesh refinement study?
flag_refinement = 1; 

% setup parameters
parameter.D = 1; 
parameter.nu = 0.22; 

parameter.bcType = 3;% bcType: 1 simply supported, 2 clamped, 3 free
parameter.mixedbcType =2 ; % mixedbcType: 1 simply supported, 2 clamped

% specify temperature type for T
% T_Type: 
% 1 uniform heating P = -1
% 2 uniform heating P = -10
% 3 parabolic localized forcing at (L*0.5,L*0.5)
% 4 parabolic localized forcing at (L*0.25,L2*0.25)

T_Type = 1;

domain = [0,2,0,1]; %rectangular domain
nx_refined=769;
nx_coarse = 13;
n_iter = 6;

% domain = [0,4,0,1]; %rectangular domain
% nx_refined=1281; 
% n_iter = 5;
% nx_coarse = 41;

domain_l = domain(2)-domain(1);
domain_h = domain(4)-domain(3);

% build clamped on part of edges
mixed_seg = [0.25*domain_l,0.75*domain_l];
mixed_seg_h = [0.25*domain_h,0.75*domain_h];
nomixed_seg = [0,-2];
nomixed_seg_h = [0,-2];

% mixed_seg = partially clamped on that edge
% nomixed_seg = BC does not change on that edge

mixedBoundary.domainL = nomixed_seg_h;
mixedBoundary.domainR = nomixed_seg_h;
mixedBoundary.domainT = mixed_seg;
mixedBoundary.domainB = mixed_seg;

%% obtain a relative accurate solution by refined grid

ny_refined = floor((nx_refined-1)*domain_h/domain_l)+1;

myGrid_refined = buildGrid(domain,nx_refined,ny_refined,mixedBoundary);
Xvec_refined = myGrid_refined.XX(:);%column vector
Yvec_refined = myGrid_refined.YY(:);%column vector

parameter.hx = myGrid_refined.hx;
parameter.hy = myGrid_refined.hy;

Index_refined=getIndex(nx_refined,ny_refined,myGrid_refined,mixedBoundary);

LapT = forcing(Xvec_refined,Yvec_refined,domain,T_Type);

RHS = LapT;

hx_refined = myGrid_refined.hx;
hy_refined = myGrid_refined.hy;
mtx = getDiffMatrix(nx_refined,ny_refined,hx_refined,hy_refined);
A = mtx.BiDh;
A=assignBoundaryConditionsCoefficient(A,Index_refined,parameter);
A = (parameter.D)*A;

RHS = assignBoundaryConditionsRHS(RHS,Index_refined,parameter);
RHS = RHS*hx_refined^2*hy_refined^2;

Uexact =  A\RHS;

figure(1)
FS = 20;
Xplot = reshape(Xvec_refined,ny_refined,nx_refined);
Yplot = reshape(Yvec_refined,ny_refined,nx_refined);
contour(Xplot,Yplot,reshape(Uexact,ny_refined,nx_refined));hold on;
%surf(Xplot,Yplot,reshape(Uexact,ny_refined,nx_refined));
set(gca,'FontSize',20)
title('Exact Solution','FontSize',FS)
colorbar('FontSize',FS)
%shading interp
colorbar('FontSize',13);
if(mixedBoundary.domainB(2)>mixedBoundary.domainB(1))
    plot(mixed_seg,[domain(3),domain(3)],'k-','LineWidth',8), hold on;
end
if(mixedBoundary.domainT(2)>mixedBoundary.domainT(1))
    plot(mixed_seg,[domain(4),domain(4)],'k-','LineWidth',8), hold on;
end
if(mixedBoundary.domainL(2)>mixedBoundary.domainL(1))
    plot([domain(1),domain(1)],mixed_seg_h,'k-','LineWidth',8), hold on;
end
if(mixedBoundary.domainR(2)>mixedBoundary.domainR(1))
    plot([domain(2),domain(2)],mixed_seg_h,'k-','LineWidth',8), hold on;
end


figure(20)
surf(Xplot,Yplot,reshape(LapT,ny_refined,nx_refined));
title('Laplace(T)','FontSize',FS)

if flag_refinement
    disp('Start mesh refinement study:');
    %% obtain a relative accurate solution by refined grid/2
    nx = (nx_refined+1)/2;
    ny = (ny_refined+1)/2;

    myGrid = buildGrid(domain,nx,ny,mixedBoundary);
    Xvec = myGrid.XX(:);%column vector
    Yvec = myGrid.YY(:);%column vector

    parameter.hx = myGrid.hx;
    parameter.hy = myGrid.hy;

    Index=getIndex(nx,ny,myGrid,mixedBoundary);

    LapT = forcing(Xvec,Yvec,domain,T_Type);

    RHS = LapT;

    hx = myGrid.hx;
    hy = myGrid.hy;
    mtx = getDiffMatrix(nx,ny,hx,hy);
    A = mtx.BiDh;
    A=assignBoundaryConditionsCoefficient(A,Index,parameter);
    A = (parameter.D)*A;

    RHS = assignBoundaryConditionsRHS(RHS,Index,parameter);
    RHS = RHS*hx^2*hy^2;
    U =  A\RHS;

    Index_refined.compare=Index_refined.all(linspace(1,ny_refined,ny),linspace(1,nx_refined,nx));
    Index_refined.compare=Index_refined.compare(:);

    Uexact = (Uexact(Index_refined.compare)*4 - U)/3;
    Index_refined = Index;
    ny_refined = ny;
    nx_refined = nx;

    %% start mesh refinement study
    nx = nx_coarse;
    ny = floor((nx-1)*domain_h/domain_l)+1;

    dx_vec = zeros(n_iter,1);
    err_vec = zeros(n_iter,1);
    errL2_vec = zeros(n_iter,1);
    
    for i_iter = 1:n_iter

        myGrid = buildGrid(domain,nx,ny,mixedBoundary);
        Xvec = myGrid.XX(:);%column vector
        Yvec = myGrid.YY(:);%column vector

        Index=getIndex(nx,ny,myGrid,mixedBoundary);

        LapT = forcing(Xvec,Yvec,domain,T_Type);

        RHS = LapT;
        RHS = RHS*myGrid.hx^2*myGrid.hy^2;


        hx = myGrid.hx;
        hy = myGrid.hy;
        parameter.hx = hx;
        parameter.hy = hy;

        mtx = getDiffMatrix(nx,ny,hx,hy);
        A = (parameter.D)*mtx.BiDh;

        A=assignBoundaryConditionsCoefficient(A,Index,parameter);
        RHS = assignBoundaryConditionsRHS(RHS,Index,parameter);

        U = A\RHS;

        Index_refined.compare=Index_refined.all(linspace(1,ny_refined,ny),linspace(1,nx_refined,nx));
        Index_refined.compare=Index_refined.compare(:);

        UUexact = Uexact(Index_refined.compare);
        err = U-UUexact;

        maxErr = max(abs(err));
        Err_L2 = norm(abs(err),2)*sqrt(hx*hy);

        dx_vec(i_iter) = max(hx,hy);
        err_vec(i_iter) = maxErr;    
        errL2_vec(i_iter) = Err_L2;

        nx = nx*2-1;
        ny = ny*2-1;

    end

    nx = (nx+1)/2;
    ny = (ny+1)/2;

    %FontSize
    FS = 20;
    figure(200)
    Xplot = reshape(Xvec,ny,nx);
    Yplot = reshape(Yvec,ny,nx);
%     contour(Xplot,Yplot,reshape(U,ny,nx));
    surf(Xplot,Yplot,reshape(U,ny,nx));
    colorbar('FontSize',FS)
    title('w','FontSize',FS)
    set(gca,'FontSize',20)
    shading interp


    figure(300)
    contour(Xplot,Yplot,reshape(err,ny,nx));hold on;
    %surf(Xplot,Yplot,reshape(err,ny,nx));
    colorbar('FontSize',FS)
    title('Error','FontSize',FS)
    set(gca,'FontSize',20)
    % shading interp

    figure
    y1 = 0.01*dx_vec;
    y2 = 0.2*dx_vec.^2;
    y3 = 0.2*dx_vec.^1.5;

    %loglog(dx_vec,err_vec,'b-*');hold on;
    %loglog(dx_vec,errL1_vec,'g-o');hold on;
    loglog(dx_vec,errL2_vec,'r-*','DisplayName','norm(err,2)'),
    legend('-DynamicLegend'),hold on;

    % loglog(dx_vec,y1,'r--','DisplayName','O(x)'),hold on;
    loglog(dx_vec,y2,'b--','DisplayName','O(x^2)'),hold on;
    %loglog(dx_vec,y3,'k--','DisplayName','O(x^{1.5})');
    set(gca,'FontSize',15)
    xlabel('dx','FontSize',FS);
    ylabel('|w - w_{exact}|_2','FontSize',FS);

    for i = 1:length(dx_vec)
        fprintf(out,'%12.10f %12.10f \n',dx_vec(i),errL2_vec(i));
    end
end



