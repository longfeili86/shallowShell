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
if(strcmp(parameters.solver,'imPicard'))
    addNonlinearTerm=false; % for imPicard solver,we add the NL term to the LHS
end
if(~parameters.isLinear && addNonlinearTerm) % for nonlinear problems
    RHS = RHS + getLOperator(mtx,W,PHI);
end

RHS=assignBoundaryConditionsRHS(RHS,Index,parameters);

end