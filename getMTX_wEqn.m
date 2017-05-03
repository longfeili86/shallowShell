function  A = getMTX_wEqn(Index,mtx,parameters,PHI,quiet)
% get the matrix for w equation according to differerent solver and bc:
%
% -- Longfei Li

if(nargin==4)
    quiet=false;
end

infoPrefix = '--getMTX_wEqn--: '; % all info displayed by this function includes this prefix

A=mtx.BiDh;


% we add -Lmatrix to A for imPicard solver and nonlinear problem
if(strcmp(parameters.solver,'imPicard') && ~parameters.isLinear)
    A = A-getLMatrix(mtx,PHI);
    if(~quiet)
        fprintf('%sModify matrix of the w equation for implicit picard iteration\n',infoPrefix);
    end
end

A = assignBoundaryConditionsCoefficient(A,Index,mtx,parameters);

end