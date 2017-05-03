function L=assignBoundaryConditionsLMatrix(L,Index,parameters)
%assign boundary conditions for the L matrix
%    L matrices are associated with the derivatives of the rhs. 
%    So the bc for L should be 0 zeros for the boundary nodes
%Input:
% L: L matrices that are components of the Jacobean matrix. 
% Index: 
% parameters:
%Output:
% L: with bc implemented


%simplify names
bcType = parameters.bcType;
L(Index.GhostAll,:) = 0.; 
if (bcType==1) % simply support
    L(Index.Boundary,:) = 0.;   % w = 0 on boundary, so rhs=0 here
elseif (bcType==2) % clamped edge
    L(Index.Boundary,:) = 0.;   % w = 0 on boundary, so rhs=0 here
elseif (bcType==3) % free edge
    % all the ghost are 0. already set!
end




end