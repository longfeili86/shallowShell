function RHS=getRHS_wEqn(W,PHI,F,W0,mtx,parameters,Index)
% get the RHS for the w eqn: \nabla^4 w   = L[w,phi]*lambda + L[w0,phi] + f.w
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


RHS = getLOperator(mtx,W0,PHI)+F;

addNonlinearTerm = true;
nonlinearFactor=1.; % the portion of the nonlinear term that will be added to the RHS
if(strcmp(parameters.solver,'imPicard'))
    implicitFactor=parameters.implicitFactor;
    nonlinearFactor=1.-implicitFactor; % for imPicard solver,we alread added  the implicitFactor*NL term to the LHS
end
if(~parameters.isLinear) % for nonlinear problems
    RHS = RHS + nonlinearFactor*getLOperator(mtx,W,PHI);
end

RHS=assignBoundaryConditionsRHS(RHS,Index,parameters);

end