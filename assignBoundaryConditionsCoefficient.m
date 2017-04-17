function A=assignBoundaryConditionsCoefficient(A,Index,diffMtx,parameter)
%assign boundary conditions for the coefficient matrix
%Input:
% A: differentiation matrix for biharmonic operator
% Index: 
% diffMtx: differentiation matrices
% bcType: boundary condition types: Supported=1, clamped=2, free=3


%simplify names
bcType = parameter.bcType;
nu = parameter.nu;  
D = parameter.D;
I = speye(size(A));
A(Index.GhostAll,:) = I(Index.GhostAll,:); 
if (bcType==1) % simply support
    A(Index.Boundary,:) = I(Index.Boundary,:);   % w = 0 on boundary
    A(Index.Ghost1,:) = 0; % zero out existing coefficient
    A(Index.GhostL1,:)=  diffMtx.Dxx(Index.BoundaryL,:) ; %l : Dxx
    A(Index.GhostR1,:)=  diffMtx.Dxx(Index.BoundaryR,:) ; %r : Dxx
    A(Index.GhostB1,:)=  diffMtx.Dyy(Index.BoundaryB,:) ; %b : Dyy  
    A(Index.GhostT1,:)=  diffMtx.Dyy(Index.BoundaryT,:) ; %t : Dyy
    %extrapolate the second ghostline
    A(Index.Ghost2,:) = 0; % zero out existed coefficient
    for i = 1:length(Index.Ghost2)
        A(Index.Ghost2(i),[Index.Ghost2(i),Index.Ghost1(i),...
            Index.Boundary(i),Index.Ghostin1(i)])=[1,-3,3,-1]; %third order extrapolation
    end
    
elseif (bcType==2) % clamped edge
    A(Index.Boundary,:) = I(Index.Boundary,:);   % w = 0 on boundary
    A(Index.Ghost1,:) = 0; % zero out existing coefficient
    A(Index.GhostL1,:)= -diffMtx.D0x(Index.BoundaryL,:) ; %l :n*grad=(-1,0)*(D0x,D0y)
    A(Index.GhostR1,:)=  diffMtx.D0x(Index.BoundaryR,:) ; %r :n*grad=(1,0)*(D0x,D0y)
    A(Index.GhostB1,:)= -diffMtx.D0y(Index.BoundaryB,:) ; %b :n*grad=(0,-1)*(D0x,D0y)  
    A(Index.GhostT1,:)=  diffMtx.D0y(Index.BoundaryT,:) ; %t :n*grad=(0,1)*(D0x,D0y)
    %extrapolate the second ghostline
    A(Index.Ghost2,:) = 0; % zero out existed coefficient
    for i = 1:length(Index.Ghost2)
        A(Index.Ghost2(i),[Index.Ghost2(i),Index.Ghost1(i),...
            Index.Boundary(i),Index.Ghostin1(i)])=[1,-3,3,-1]; %third order extrapolation
    end
elseif (bcType==3) % free edge
    A(Index.Ghost1,:) = 0; % zero out existing coefficient
    A(Index.Ghost2,:) = 0; % zero out existing coefficient
    % left boundary
    A(Index.GhostL1,:)= diffMtx.Dxx(Index.BoundaryL,:)+nu*diffMtx.Dyy(Index.BoundaryL,:) ; 
    A(Index.GhostL2,:)= diffMtx.Dxxx(Index.BoundaryL,:)+(2-nu)*diffMtx.Dyyx(Index.BoundaryL,:) ;
    % right boundary
    A(Index.GhostR1,:)= diffMtx.Dxx(Index.BoundaryR,:)+nu*diffMtx.Dyy(Index.BoundaryR,:) ; 
    A(Index.GhostR2,:)= diffMtx.Dxxx(Index.BoundaryR,:)+(2-nu)*diffMtx.Dyyx(Index.BoundaryR,:) ; 
    % top bboundary
    A(Index.GhostT1,:)= diffMtx.Dyy(Index.BoundaryT,:)+nu*diffMtx.Dxx(Index.BoundaryT,:) ; 
    A(Index.GhostT2,:)= diffMtx.Dyyy(Index.BoundaryT,:)+(2-nu)*diffMtx.Dxxy(Index.BoundaryT,:) ;
    % bottom bboundary
    A(Index.GhostB1,:)= diffMtx.Dyy(Index.BoundaryB,:)+nu*diffMtx.Dxx(Index.BoundaryB,:) ; 
    A(Index.GhostB2,:)= diffMtx.Dyyy(Index.BoundaryB,:)+(2-nu)*diffMtx.Dxxy(Index.BoundaryB,:) ;   
end




end