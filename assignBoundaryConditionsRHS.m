function RHS=assignBoundaryConditionsRHS(RHS,Index,parameter)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this function assigns the rhs for corresponding boundary conditions
% Input:
%   RHS: the RHS vector of the system
%   Index: index object
%   parameter: parameter object
% Output:
% Output:
%   RHS: the RHS vector with BC implemented
% 
% by L. Li 07/01/2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    RHS(Index.GhostAll) = 0;

    if(parameter.bcType==1) % simply supported 
        RHS(Index.Boundary) = 0;
        RHS(Index.GhostL1) = 0; 
        RHS(Index.GhostR1) =  0; 
        RHS(Index.GhostB1) = 0;    
        RHS(Index.GhostT1) =  0;
        RHS(Index.Ghost2) = 0; %extrapolation rhs=0
    elseif(parameter.bcType==2) % clamped edge 
        RHS(Index.Boundary) = 0;
        RHS(Index.GhostL1) = 0; 
        RHS(Index.GhostR1) =  0; 
        RHS(Index.GhostB1) = 0;    
        RHS(Index.GhostT1) =  0;
        RHS(Index.Ghost2) = 0; %extrapolation rhs=0
    elseif(parameter.bcType==3) % free edge 
        %RHS(Index.GhostAll)=100;
        %left boundary
        RHS(Index.GhostL1)=0;
        RHS(Index.GhostL2)=0;
        %right boundary   
        RHS(Index.GhostR1)=0;
        RHS(Index.GhostR2)=0;
        %top boundary   
        RHS(Index.GhostT1)=0;
        RHS(Index.GhostT2)=0;
        %bottom boundary   
        RHS(Index.GhostB1)=0;
        RHS(Index.GhostB2)=0;
        %corner conditions for free bc
        RHS(Index.UsedGhostCorners)=0.;
    end


end