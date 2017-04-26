function RHS=getRHS_phiEqn(W,PHI,F,W0,mtx,parameters,Index)
% get the RHS for the phi eqn: \nabla^4 phi = -1/2 L[w,w]*lambda - L[w0,w] -f.phi
% Input:
%       W: vectorized w
%     PHI: vectorized phi
%       F: external forcing (vectorized) 
%      W0: initla shell shape (vectorized) 
%     mtx: diff matrix
%   parameters:
% Output:
%    RHS: rhs of the phi eqn
%
% -- Longfei Li

RHS  = -getLOperator(mtx,W0,W)-F;

if(~parameters.isLinear) % for nonlinear problems
    RHS = RHS -0.5*getLOperator(mtx,W,W);
end


bcTypeSaved=parameters.bcType; % save bcType

if(parameters.bcType==3)
    % for free bc, the phi eqn uses the clamped bc: phi=dphidn=0
    % so we overwrite the bcType value here to reuse
    % the assignBoundaryConditionsRHS fucntion for the phi eqn as well
    parameters.bcType=2; 
end

RHS=assignBoundaryConditionsRHS(RHS,Index,parameters);
parameters.bcType=bcTypeSaved; % reset bcType to saved value


end