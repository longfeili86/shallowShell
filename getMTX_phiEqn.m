function  A = getMTX_phiEqn(Index,mtx,parameters)
% get the matrix for phi equation according to differerent solver and bc:
%
% -- Longfei Li

A=mtx.BiDh;

bcTypeSaved=parameters.bcType; % save bcType
if(parameters.bcType==3)
    % for free bc, the phi eqn uses the clamped bc: phi=dphidn=0
    % so we overwrite the bcType value here to reuse
    % the assignBoundaryConditionsCoefficient fucntion for the phi eqn as well
    parameters.bcType=2; 
end

A = assignBoundaryConditionsCoefficient(A,Index,mtx,parameters);

parameters.bcType=bcTypeSaved; % reset bcType to saved value

end