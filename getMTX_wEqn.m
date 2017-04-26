function  A = getMTX_wEqn(Index,mtx,parameters,PHI)
% get the matrix for w equation according to differerent solver and bc:
%
% -- Longfei Li

A=mtx.BiDh;

% we add -Lmatrix to A for imPicard solver and nonlinear problem
if(strcmp(parameters.solver,'imPicard') && ~parameters.isLinear)
    A = A-getLMatrix(mtx,PHI);
end

A = assignBoundaryConditionsCoefficient(A,Index,mtx,parameters);

end