%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script is for snap-through bifurcation study of the following coupled
% nonlinear system (DIMENSIONLESS) with critical thermal forcing 
% \nabla^2(T) = \xi
%
%  \nabla^4 w(x,y) = [phi,w0] + [phi,w]
%  \nabla^4 phi(x,y) = -0.5[w,w]-[w0,w] - \nabla^2(T)
% 
% Algorithm:
% Step 1:
%  Solve the following coupled linear system to get initial guess (phi,w)
%  \nabla^4 w(x,y) = [phi,w0]
%  \nabla^4 phi(x,y) = -[w0,w] - \nabla^2(T)
%  is solved first to get initial guess for the nonlinear system solutions
%
% Step 2:
% Solve the nonlinear system with Picard-type iterative method 
% by solving the first and second equation iteratively.
% 
% Clamped BC : w = w_n = 0
% or simply supported BC: w = w_nn = 0
% imposed on C for w
%
% Free boundary conditions imposed on AllEdges\C for w
% Left & Right: w_xx+nu*w_yy=0, w_xxx+(2-nu)*w_xyy=0
% Top & Bottom: w_yy+nu*w_xx=0, w_yyy+(2-nu)*w_xxy=0
%
% Clamped or free BC:  phi = phi_n = 0
%   simply supported BC: phi = phi_nn = 0
%
% corner conditions: w_{xy}=0
%
% Used continuation to find buckling (multiple solns) with elliptic
% paraboloid w0 and mixed boundary conditions
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
% by H.Ji 4/14/2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Specify unstressed shell shape in getInitialShape.m.
% Specify Temperature function in getTempDist.m

%% Start calculation
close all
clear
clc

out = fopen('out_bifurcation.dat','w');
out2 = fopen('out_solution.dat','w');

%============================ Physical parameters=============================
% For nondimensional calculation, only Poisson's ratio \nu is needed

% parameter.E = 7e4; % Young's Modulus
% parameter.h = 0.5; % thickness
parameter.nu=0.22; % Poisson's ratio
% parameter.alpha = 5.0e-6; % Coefficient of thermal expansion
% 
% parameter.D = parameter.E*parameter.h^3/(1.0-parameter.nu^2);
% domain size = 1000*1000 (mm^2)

err_tol = 1e-6;


% specify bcType for phi and w here
% clamped bc for phi, do not change it.
% bcType: 1 simply supported, 2 clamped, 3 free
% mixedbcType: 1 simply supported, 2 clamped, 0 no mixed bc
phi_bcType = 1;
phi_mixedbcType = 2;
w_bcType = 1;
w_mixedbcType =2;

% specify temperature type for T
% T_Type: 
% 1 uniform heating P = -1
% 2 uniform heating P = -10
% 3 parabolic localized forcing at (L*0.5,L*0.5)
% 4 parabolic localized forcing at (L*0.75,L2*0.25)

T_Type = 1;

% specify w0 type for phi and w here
% w0_Type:  
%           0 uniform W0;
%           1 elliptic paraboloid W0; 
%           2 hyperbolic paraboloid W0;
%           3 cylindrical shell W0;
%           4 Gaussian shape W0;

w0_Type = 1;

isPlot_w = 1;
isPlot_phi = 1;

isPlot_w0 = 1; %imposed w0.
isPlot_LapT = 1; %imposed nabla^2(T);
isPlot_T = false; %imposed T;

% buildGrid

domain = [0,1,0,1]; % square domain

domain_l = domain(2)-domain(1);
domain_h = domain(4)-domain(3);

% build clamped bc on part of edges

mixed_seg = [0.4*domain_l,0.6*domain_l];
mixed_seg_h = [0.25*domain_h,0.75*domain_h];
nomixed_seg = [0,-2];
nomixed_seg_h = [0,-2];

mixedBoundary.domainL = nomixed_seg_h;
mixedBoundary.domainR = nomixed_seg_h;
mixedBoundary.domainT = mixed_seg;
mixedBoundary.domainB = mixed_seg;


flag_afterbif = 0;

flag_unstablebranch = 0; % use semi-implicit iterative scheme to get the unstable branch

LapT_target = 3000; % target \nabla^2 T
 
if(~flag_unstablebranch)
 dLapT = 50; 
 LapTInit = LapT_target; 

else
 LapTInit =LapT_target;
 dLapT = 20;
%  dLapT = -20;
end

%==========================================================================

%% problem setup
% buildGrid

nx=81;
%nx=161;

ny = floor((nx-1)*domain_h/domain_l)+1;

myGrid = buildGrid(domain,nx,ny,mixedBoundary);
Xvec = myGrid.XX(:);%column vector
Yvec = myGrid.YY(:);%column vector

getInitialShape;

LapT = forcing(Xvec,Yvec,domain,T_Type);

LapT = LapTInit*LapT; % rescale for appropriate domain
LapT_rec = LapT;

w0 = Ue(Xvec,Yvec); % initial guess for deflection W, fixed in this model
% Phi0 = Phi0e(Xvec,Yvec); % initial guess for the unknown Phi.

%w0 = w0*0.3; % shell with small curvature;
%w0 = w0*3; % for shell with large curvature precast shape

w0 = w0*0.1;

if (w0_Type == 0)
    w0 = ones(size(Xvec))*w0;
end

Index=getIndex(nx,ny,myGrid,mixedBoundary);

if (isPlot_LapT)
    FS = 20;
    figure(99);
    Xplot = reshape(Xvec,nx,ny);
    Yplot = reshape(Yvec,nx,ny);
    Zplot = reshape(LapT,nx,ny);

    surf(Xplot,Yplot,Zplot);
    title('imposed Laplace(T)','FontSize',FS)
    set(gca,'FontSize',20)
    shading interp
end

if (isPlot_w0)
    FS = 20;
    figure(101);
    Xplot = reshape(Xvec,nx,ny);
    Yplot = reshape(Yvec,nx,ny);
    Zplot = reshape(w0,nx,ny);

    surf(Xplot,Yplot,Zplot);
    title('imposed w0 (unstressed shell shape)','FontSize',FS)
    set(gca,'FontSize',20)
    colorbar('FontSize',20)
    shading interp
end

hx = myGrid.hx;
hy = myGrid.hy;
parameter.hx = myGrid.hx;
parameter.hy = myGrid.hy;
mtx = getDiffMatrix(nx,ny,hx,hy);


%% start iteratively solving the coupled linear system
% use its solution as the initial guess of the nonlinear system

flag_iterate = 1;

if (w0_Type~=0)
    
    w = w0; % initial guess
    

    phi = 0*Xvec;
    disp('Start solving the coupled linear system');
    
    phi_old =phi;

    for ii = 1:100
        
        %plot
        %FontSize
        w_old = w;

        if (isPlot_w)
            FS = 20;
            figure(11);
            Xplot = reshape(Xvec,ny,nx);
            Yplot = reshape(Yvec,ny,nx);
            Zplot = reshape(w,ny,nx);

            surf(Xplot,Yplot,Zplot);view(2);
            str = sprintf('w initial guess from linear solver');
            title(str,'FontSize',FS)
            set(gca,'FontSize',20)
            shading interp
             pause(.1)
        end  
       
        % set clamped edge for stress function phi
        % bcType: 1 simply supported, 2 clamped, 3 free
        % mixedbcType: 1 simply supported, 2 clamped, 0 no mixed bc

        parameter.bcType = phi_bcType; % bcType: 2 clamped edge, 3 free edge
        parameter.mixedbcType = phi_mixedbcType;

        Lw0w = getLOperator(mtx,w0,w);
        A_phi = mtx.BiDh;
        RHS_phi = - Lw0w/hx^2/hy^2 - LapT;

        A_phi = assignBoundaryConditionsCoefficient(A_phi,Index,parameter);
        RHS_phi = assignBoundaryConditionsRHS(RHS_phi,Index,parameter);

        RHS_phi = RHS_phi*hx^2*hy^2;
        phi = A_phi\RHS_phi;  
        
        err_phi = norm(phi_old-phi,2)/norm(phi,2); 
        phi_old = phi;

        if (isPlot_phi)
            FS = 20;
            figure(23);
            Xplot = reshape(Xvec,ny,nx);
            Yplot = reshape(Yvec,ny,nx);
            Zplot = reshape(phi,ny,nx);
            surf(Xplot,Yplot,Zplot);view(2);
            str = sprintf('phi, after %d iteration(s)',ii);
            title(str,'FontSize',FS)
            set(gca,'FontSize',20)
            shading interp
            pause(.1)
        end

        % set bc for deflection w
        % bcType: 1 simply supported, 2 clamped, 3 free
        % mixedbcType: 1 simply supported, 2 clamped, 0 no mixed bc
        parameter.bcType = w_bcType; 
        parameter.mixedbcType = w_mixedbcType;

        RHS_w = getLOperator(mtx,phi,w0)/hx^2/hy^2;
        A_w = mtx.BiDh;
        A_w = assignBoundaryConditionsCoefficient(A_w,Index,parameter);
        RHS_w = assignBoundaryConditionsRHS(RHS_w,Index,parameter);

        RHS_w = RHS_w*hx^2*hy^2;

        w = A_w\RHS_w;    
        err_w = norm(w_old-w,2)/norm(w,2);

         fprintf('After %d iteration(s): |w - w_old|_2/|w| = %.4g, |phi - phi_old|_2/|phi| = %.4g\n',ii,err_w,err_phi); 

        if (abs(err_w)<err_tol*0.1 && abs(err_phi)<err_tol*0.1)
            break;
        end


    end

    if (abs(err_w)<err_tol && abs(err_phi)<err_tol)
        w_linear = w;
        phi_linear = phi;
    else
        disp('The iterative solver does not converge for the linear system.');
        flag_iterate = 0;
        return;
    end
    
    fprintf('max(w0)=%g\n',max(w0)); 
end


%% start iteratively solving the nonlinear system

w = w_linear; % Initial guess
phi = phi_linear;

n_iter = 100;

w_rec = w;
phi_rec=phi;

center_vec = zeros(1,1);
j=1;


BAD=0;

while (BAD<20)

    %% start iteratively solving the nonlinear system
    disp(' ');
    disp('==============Start solving the coupled nonlinear system==============');
    fprintf('max(LapT) = %g, dLapT = %g\n\n',max(LapT),dLapT); 

    for ii = 1:n_iter
        %plot
        %FontSize
        w_old = w;
        phi_old = phi;

        if (isPlot_w)
            FS = 20;
            figure(12);
            Xplot = reshape(Xvec,ny,nx);
            Yplot = reshape(Yvec,ny,nx);
            Zplot = reshape(w,ny,nx);

            surf(Xplot,Yplot,Zplot);view(2);
            str = sprintf('w, after %d iteration(s)',ii);
            title(str,'FontSize',FS)
            colorbar('FontSize',FS)
            set(gca,'FontSize',20)
             shading interp
             pause(.1)
        end

        % set clamped edge for stress function phi
        % bcType: 1 simply supported, 2 clamped, 3 free
        % mixedbcType: 1 simply supported, 2 clamped, 0 no mixed bc

        parameter.bcType = phi_bcType; % bcType: 2 clamped edge, 3 free edge
        parameter.mixedbcType = phi_mixedbcType;

        Lw0w = getLOperator(mtx,w0,w);
        Lww = getLOperator(mtx,w,w);
        A_phi = mtx.BiDh;
        
        RHS_phi = -0.5*Lww/hx^2/hy^2-Lw0w/hx^2/hy^2 - LapT;

        A_phi = assignBoundaryConditionsCoefficient(A_phi,Index,parameter);
        RHS_phi = assignBoundaryConditionsRHS(RHS_phi,Index,parameter);

        RHS_phi = RHS_phi*hx^2*hy^2;   

        phi = A_phi\RHS_phi;

        err_phi = norm(phi_old-phi,2)/norm(phi,2);
        phi_old = phi;

       if (isPlot_phi)
            FS = 20;
            figure(23);
            Xplot = reshape(Xvec,ny,nx);
            Yplot = reshape(Yvec,ny,nx);
            Zplot = reshape(phi,ny,nx);
            surf(Xplot,Yplot,Zplot);view(2);
            colorbar('FontSize',FS)
            str = sprintf('phi, after %d iteration(s)',ii);
            title(str,'FontSize',FS)
            set(gca,'FontSize',20)
            shading interp
            pause(.1)
        end
        

        % set bc for deflection w
        % bcType: 1 simply supported, 2 clamped, 3 free
        % mixedbcType: 1 simply supported, 2 clamped, 0 no mixed bc
        parameter.bcType = w_bcType; 
        parameter.mixedbcType = w_mixedbcType;

         
         if(flag_unstablebranch)   

             L0 = sparse(length(phi),length(phi));

            L_coeff = spdiags(mtx.Dxx*phi,0,L0)*mtx.Dyy+spdiags(mtx.Dyy*phi,0,L0)*mtx.Dxx...
            -2.0*spdiags(mtx.Dxy*phi,0,L0)*mtx.Dxy;

            RHS_w = getLOperator(mtx,phi,w0)/hx^2/hy^2;
            
            A_w = mtx.BiDh - L_coeff;
            A_w = assignBoundaryConditionsCoefficient_implicit(A_w,Index,parameter);
            
            RHS_w = assignBoundaryConditionsRHS(RHS_w,Index,parameter);    
            RHS_w = RHS_w*hx^2*hy^2;
            
         else
             RHS_w = getLOperator(mtx,phi,w0)+getLOperator(mtx,phi,w);
             A_w = mtx.BiDh;
             A_w = assignBoundaryConditionsCoefficient(A_w,Index,parameter);
             RHS_w = assignBoundaryConditionsRHS(RHS_w,Index,parameter);   
         end
             
          w = A_w\RHS_w;

        err_w = norm(w_old-w,2)/norm(w,2);

         fprintf('After %d iteration(s): |w - w_old|_2/|w| = %.4g, |phi - phi_old|_2/|phi| = %.4g\n',ii,err_w,err_phi); 

        if (abs(err_w)<err_tol*0.1 && abs(err_phi)<err_tol*0.1)
            break;
        end
        
         if (isnan(err_w) || isnan(err_phi))
             break;
         end
    end

    if (abs(err_w)<err_tol && abs(err_phi)<err_tol)

         if(abs(max(LapT)+LapT_target)<10 && flag_afterbif)
             return;
         end
        
        center_vec(j)=w((nx*ny+1)/2);
        
        w_rec_rec = w_rec;
        phi_rec_rec = phi_rec;
        w_rec = w;
        phi_rec = phi;
        w_nonlinear = w;
        phi_nonlinear = phi;
        LapT_rec = LapT;
        BAD = 0;
        
        fprintf(out,'%12.10f %12.10f %d %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f\n',...
            max(LapT),center_vec(j),ii,sum(phi)*hx*hy, sum(w)*hx*hy, norm(phi,2)*sqrt(hx*hy), ...
            norm(w,2)*sqrt(hx*hy),w((nx*ny+1)/2)+w0((nx*ny+1)/2), norm(w+w0,2)*sqrt(hx*hy));

        if (center_vec(j)<0  && flag_unstablebranch==0)
            dLapT_old = dLapT;
            %% uncomment these two to get the lower branch of the second branch in the bifurcation diagram
            dLapT = -1.0*abs(dLapT); 
            flag_afterbif = 1;
 
            save w_phi.mat w_rec phi_rec LapT_rec dLapT_old w0 w_rec_rec phi_rec_rec;
        else
            dLapT_old = dLapT;
        end
        LapT = LapT-dLapT;
        
        j=j+1;
        
        
    else
        disp('The iterative solver does not converge for the nonlinear system.');
        BAD = BAD+1;
        fprintf('BAD = %d\n',BAD); 
        fprintf('flag_afterbif = %d\n',flag_afterbif);
        fprintf('w((nx*ny+1)/2) = %.8f\n',w((nx*ny+1)/2));
        fprintf('w_rec((nx*ny+1)/2) = %.8f\n',w_rec((nx*ny+1)/2));
    end

    if(isPlot_w)
        FS = 20;
        figure(12);
        Xplot = reshape(Xvec,ny,nx);
        Yplot = reshape(Yvec,ny,nx);
        Zplot = reshape(w,ny,nx);

        surf(Xplot,Yplot,Zplot);view(2);
        str = sprintf('nonlinear w, LapT=%.4f',max(LapT));
        title(str,'FontSize',FS)
        colorbar('FontSize',FS)
        set(gca,'FontSize',20)
        shading interp
        pause(.1)
    end
    
    
      
    if (BAD>0) 
        
       if (w_rec((nx*ny+1)/2)<0 && (flag_afterbif || flag_unstablebranch))
           dLapT = dLapT*0.8;
           w = w_rec+(w_rec-w_rec_rec)/(dLapT_old)*dLapT;
           phi = phi_rec+(phi_rec-phi_rec_rec)/dLapT_old*dLapT;
       else 
            dLapT = dLapT*0.5;   
            w = w_rec;
            phi = phi_rec;
            w = -1.0*w;
       end 
       
       LapT = LapT_rec-dLapT;
        
    end
    
   
    
    
end


fclose(out);

out2 = fopen('out_solution.dat','w');
Xplot = reshape(Xvec,ny,nx);
Yplot = reshape(Yvec,ny,nx);
w_plot = reshape(w,ny,nx);
w_plot = w_plot(floor((ny-1)/2),:);

for i = 1:length(w_plot)

fprintf(out2,'%12.10f %12.10f %12.10f\n',...
    Xplot(floor((ny-1)/2),i),Yplot(floor((ny-1)/2),i), w_plot(1,i));

end
fclose(out2);


%% Plot of solutions (w,phi)
FS = 20;
figure;
Xplot = reshape(Xvec,ny,nx);
Yplot = reshape(Yvec,ny,nx);
Zplot = reshape(w,ny,nx);
surf(Xplot,Yplot,Zplot),hold on;
%contour(Xplot,Yplot,Zplot),hold on;view(2);
str = sprintf('w(x,y)',ii);
title(str,'FontSize',FS)
colorbar('FontSize',FS)
set(gca,'FontSize',20)
shading interp
pause(.1)

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

FS = 20;
figure;
Zplot = reshape(phi,ny,nx);
surf(Xplot,Yplot,Zplot),hold on;
%contour(Xplot,Yplot,Zplot),hold on;view(2);
str = sprintf('phi(x,y)',ii);
title(str,'FontSize',FS)
colorbar('FontSize',FS)
set(gca,'FontSize',20)
shading interp

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

FS = 20;
figure(12);
Xplot = reshape(Xvec,nx,ny);
Yplot = reshape(Yvec,nx,ny);
Zplot = reshape(w+w0,nx,ny);

surf(Xplot,Yplot,Zplot),hold on;
%contour(Xplot,Yplot,Zplot),hold on;view(2);
str = sprintf('w+w0',max(LapT));
title(str,'FontSize',FS)
colorbar('FontSize',FS)
set(gca,'FontSize',20)
shading interp

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

