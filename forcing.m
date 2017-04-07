function laplace_T = forcing(Xvec,Yvec,domain,T_Type)
% Forcing term in the biharmonic equation

L = domain(2)-domain(1);
L2 = domain(4)-domain(3);

% specify temperature type for T here
% T_Type: 
% 1 uniform heating P = -1
% 2 uniform heating P = -10
% 3 parabolic localized forcing at (L*0.5,L*0.5)
% 4 parabolic localized forcing at (L*0.75,L2*0.25)

switch T_Type
        
    case 1 % uniform heating P = -1
         laplace_T = -1.0*ones(size(Xvec)); 
         
    case 2 % uniform heating P = -1
         laplace_T = -10.0*ones(size(Xvec)); 
        
    case 3 % parabolic localized forcing at (L*0.5,L2*0.5)
        x0 = L*0.5;
        y0 = L2*0.5;
        laplace_T = -1000.0/L^2*((Xvec-x0).^2+(Yvec-y0).^2)+10.0; 
        laplace_T = max(laplace_T,0);
    case 4 % parabolic localized forcing at (L*0.75,L2*0.25)
        x0 = L*0.75;
        y0 = L2*0.25;
        laplace_T = -1000.0/L^2*((Xvec-x0).^2+(Yvec-y0).^2)+10.0; 
        laplace_T = max(laplace_T,0);
    otherwise
        disp('Incorrect thermal forcing type. T_Type supported: 1,2');
end





end

