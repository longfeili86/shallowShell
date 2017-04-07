%get the unstressed shape of the shallow shell
syms x y;% delete after use

L = domain(2)-domain(1);
S = domain(4)-domain(3);

% UU is W0, the precast unstressed shell shape for deflection W

switch w0_Type
    % w0_Type:  
    %           0 uniform W0;
    %           1 elliptic paraboloid W0; 
    %           2 hyperbolic paraboloid W0;
    %           3 cylindrical shell W0;

    case 0 % uniform W0
        UU = 1.0;

    case 1 % elliptic paraboloid W0;
         UU = 1.0*(1.0-((x-L*0.5)/S)^2-((y-L*0.5)/S)^2);
        
%         UU = 0.1*(1.0-((x-L*0.5)/S)^2-((y-L*0.5)/S)^2);
        
    case 2      % hyperbolic paraboloid W0
        UU = 10.0*(1.0-((x-L*0.5)/S)^2+((y-L*0.5)/S)^2);
     
        
    case 3 % cylindrical shell
        UU = 2.5-10.0*((y-L*0.5)/S)^2;
        
    case 4      % Gaussian shape W0
        x0 = L*0.5;
        y0 = L*0.5;
        UU = 1.0*(exp(-5.0*(x-x0)^2-5.0*(y-y0)^2));
        
        x0 = L*0.75;
        y0 = L*0.5;
        UU = 1.0*(exp(-10.0*(x-x0)^2-10.0*(y-y0)^2));
        
    otherwise
        disp('default trig surface for w0');
        % UU = sin(pi*x/(domain(2)-domain(1)))*sin(pi*y*(domain(4)-domain(3)));
        %UU = 0.1*sin(pi*x/(domain(2)-domain(1)));

        UU = sin(0.1*pi*x/(domain(2)-domain(1)))*sin(0.1*pi*y*(domain(4)-domain(3)));
end


Phi0 = 0.1*sin(x*pi/L)*sin(y*pi/L);


Biharm_UU = diff(UU,x,4)+diff(UU,y,4)+2*diff(diff(UU,x,2),y,2);
UU_x = diff(UU,x,1);
UU_y = diff(UU,y,1);
UU_xy = diff(UU_x,y,1);
UU_xx = diff(UU,x,2);
UU_yy = diff(UU,y,2);
UU_xxx = diff(UU,x,3);
UU_yyy = diff(UU,y,3);
UU_xxy = diff(UU_xx,y,1);
UU_yyx = diff(UU_yy,x,1);

Phi0_x = diff(Phi0,x,1);
Phi0_y = diff(Phi0,y,1);
Phi0_xy = diff(Phi0_x,y,1);
Phi0_xx = diff(Phi0,x,2);
Phi0_yy = diff(Phi0,y,2);
Phi0_xxx = diff(Phi0,x,3);
Phi0_yyy = diff(Phi0,y,3);
Phi0_xxy = diff(Phi0_xx,y,1);
Phi0_yyx = diff(Phi0_yy,x,1);

Loperator = UU_xx*Phi0_yy + UU_yy*Phi0_xx - 2.0*UU_xy*Phi0_xy;

Phi0e = matlabFunction(Phi0,'vars',[x,y]);

Ue = matlabFunction(UU,'vars',[x,y]);
Ue_x = matlabFunction(UU_x,'vars',[x,y]);
Ue_y = matlabFunction(UU_y,'vars',[x,y]);
Ue_xx = matlabFunction(UU_xx,'vars',[x,y]);
Ue_yy = matlabFunction(UU_yy,'vars',[x,y]);
Ue_xxx = matlabFunction(UU_xxx,'vars',[x,y]);
Ue_yyy = matlabFunction(UU_yyy,'vars',[x,y]);
Ue_xxy = matlabFunction(UU_xxy,'vars',[x,y]);
Ue_yyx = matlabFunction(UU_yyx,'vars',[x,y]);
Biharm_Ue = matlabFunction(Biharm_UU,'vars',[x,y]);

Le_operator = matlabFunction(Loperator,'vars',[x,y]);

clear x y;














    
    
  