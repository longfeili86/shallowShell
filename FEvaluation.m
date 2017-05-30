function F=FEvaluation(x,n,Aphi,Aw,RHSphi,RHSw,R,Index,parameters,xOld,dx)
% evaluate the F function
%
% optional input: nearby solutions and its derivative xOld and dx, needed by the normalization eqn
% for the PAC method
%
% -- Longfei Li


Nphi=1:n; % phi solutions
Nw=n+1:2*n;  % w solutions
Nlambda=2*n+1:2*n+3; % lambda solutions



phi=x(Nphi);
w=x(Nw);

F(Nphi,1)= Aphi*phi-RHSphi(w,phi); % phi eqn
if(parameters.bcType==3) % w eqn is singular, we need extra equations for lambda
    lambda=x(Nlambda); %  extra unknowns
    F([Nw,Nlambda],1) = Aw*[w;lambda]-[RHSw(w,phi);R]; % w + lambda eqn
else
    F(Nw,1) = Aw*w-RHSw(w,phi); % w eqn
end

% 20170523: new for pseudo-arclength continuation (PAC) method:
% if use pseudo-arclength continuation method, we add an additial unknown vxi
% to the phi eqn and an additional normalization equation
if(parameters.usePAC)
    ds=parameters.ds;
    vxi(1:length(Nphi),1)=x(end); % we need add to vxi to each of phi eqautions
    bcTypeSaved=parameters.bcType; % save bcType
    if(parameters.bcType==3 || parameters.bcType==5)
    % for free bc and CF bc, the phi eqn uses the clamped bc: phi=dphidn=0
    % so we overwrite the bcType value here to reuse
    % the assignBoundaryConditionsRHS fucntion for the phi eqn as well
        parameters.bcType=2; 
    end
    vxi=assignBoundaryConditionsRHS(vxi,Index,parameters);
    parameters.bcType=bcTypeSaved; % reset bcType to saved value

    F(Nphi,1)= F(Nphi,1)+vxi;
    nE=length(F); % number of equations
    F(nE+1,1)= dx'*(x-xOld)-ds;
end


end