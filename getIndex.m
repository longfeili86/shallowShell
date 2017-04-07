function Index = getIndex(nx,ny,grid,mixedBoundary)
% This function get the indices for the boudanry and ghost line ect.
% Input:
%    nx: number of grids in x direction       
%    ny: number of grids in y direction
% Output:
%    Index: an object that holds index info
% added used ghost corner point to the index list

% First version by L.Li
% modified by H.Ji 12/25/15

    Index.all = reshape(1:nx*ny,ny,nx);
    Index.nx = nx;
    Index.ny = ny;

    Index.interior=Index.all(2:ny-1,2:nx-1);
    Index.interior=Index.interior(:);
    
    % Indices of left,right, top and bottom boundary of a rectangular
    % domain (not including corners)
    Index.BoundaryB = Index.all(ny,2:nx-1);   Index.BoundaryB=Index.BoundaryB(:); 
    Index.BoundaryT = Index.all(1,2:nx-1);Index.BoundaryT=Index.BoundaryT(:);
    Index.BoundaryR = Index.all(2:ny-1,nx);Index.BoundaryR=Index.BoundaryR(:);
    Index.BoundaryL = Index.all(2:ny-1,1);   Index.BoundaryL=Index.BoundaryL(:);
    Index.Boundary = [Index.BoundaryL;Index.BoundaryR;Index.BoundaryB;Index.BoundaryT]; % column vector
    
    % Indices of left,right, top and bottom boundary of a rectangular
    % domain (not including corners and the point adjacent to corner on the edge)
    % This is needed for free corner conditions
    Index.BoundaryB2 = Index.all(ny,3:nx-2);   Index.BoundaryB2=Index.BoundaryB2(:); 
    Index.BoundaryT2 = Index.all(1,3:nx-2);Index.BoundaryT2=Index.BoundaryT2(:);
    Index.BoundaryR2 = Index.all(3:ny-2,nx);Index.BoundaryR2=Index.BoundaryR2(:);
    Index.BoundaryL2 = Index.all(3:ny-2,1);   Index.BoundaryL2=Index.BoundaryL2(:);
    Index.Boundary2 = [Index.BoundaryL2;Index.BoundaryR2;Index.BoundaryB2;Index.BoundaryT2]; % column vector
    
    % four corners in the domain 
    Index.CornerBL = Index.all(ny,1);
    Index.CornerBR = Index.all(ny,nx);
    Index.CornerTL = Index.all(1,1);
    Index.CornerTR = Index.all(1,nx);
    Index.Corners = [Index.CornerBL(:);Index.CornerBR(:);Index.CornerTL(:);Index.CornerTR(:)];
    
    % Eight points adjacent to corners on the edge
    Index.CornerBLL = Index.all(ny-1,1);
    Index.CornerBLB = Index.all(ny,2);
    Index.CornerBRR = Index.all(ny-1,nx);
    Index.CornerBRB = Index.all(ny,nx-1);
    Index.CornerTLL = Index.all(2,1);
    Index.CornerTLT = Index.all(1,2);
    Index.CornerTRR = Index.all(2,nx);
    Index.CornerTRT = Index.all(1,nx-1);
    Index.CornersAdjacent = [Index.CornerBLL(:),1;Index.CornerBLB(:),2;Index.CornerBRR(:),3;Index.CornerBRB(:),4;...
        Index.CornerTLL(:),5;Index.CornerTLT(:),6;Index.CornerTRR(:),7;Index.CornerTRT(:),8];
    
    Index.BoundaryCorners = [Index.BoundaryL2,ones(ny-4,1);Index.BoundaryR2,2*ones(ny-4,1);...
                        Index.BoundaryB2,3*ones(nx-4,1);Index.BoundaryT2,4*ones(nx-4,1);...
                        Index.CornerBL,5;Index.CornerBR,6;...
                        Index.CornerTL,7;Index.CornerTR,8;...
                        Index.CornerBLL,9;Index.CornerBLB,10;...
                        Index.CornerBRR,11;Index.CornerBRB,12;...
                        Index.CornerTLL,13;Index.CornerTLT,14;...
                        Index.CornerTRR,15;Index.CornerTRT,16];
    

    
    % First line in the domain adjacent to the edge (not including corners)
    Index.InB1 = Index.all(ny-1,3:nx-2);   Index.InB1=Index.InB1(:); 
    Index.InT1 = Index.all(2,3:nx-2);Index.InT1=Index.InT1(:);
    Index.InR1 = Index.all(3:ny-2,nx-1);Index.InR1=Index.InR1(:);
    Index.InL1 = Index.all(3:ny-2,2);   Index.InL1=Index.InL1(:);
    
    % four inner corners inside the corner
    Index.InCornerBL = Index.all(ny-1,2);
    Index.InCornerBR = Index.all(ny-1,nx-1);
    Index.InCornerTL = Index.all(2,2);
    Index.InCornerTR = Index.all(2,nx-1);
    Index.InCorners = [Index.InCornerBL(:);Index.InCornerBR(:);Index.InCornerTL(:);Index.InCornerTR(:)];
    
    Index.InBoundaryCorner = zeros(2*nx+2*ny-12,2);
    Index.InBoundaryCorner = [Index.InL1,ones(ny-4,1);Index.InR1,2*ones(ny-4,1);...
                        Index.InB1,3*ones(nx-4,1);Index.InT1,4*ones(nx-4,1);...
                        Index.InCornerBL,5;Index.InCornerBR,6;...
                        Index.InCornerTL,7;Index.InCornerTR,8];
                    
    Index.InCorners = [Index.InCornerBL;Index.InCornerBR;Index.InCornerTL;Index.InCornerTR];                
    Index.InBoundary = [Index.InL1;Index.InR1;Index.InB1;Index.InT1];

    
    % index for part of edges where mixed bc (clamped or simply supported) are imposed. 
    La = mixedBoundary.domainL(1); Lb = mixedBoundary.domainL(2);
    Ra = mixedBoundary.domainR(1); Rb = mixedBoundary.domainR(2);
    Ta = mixedBoundary.domainT(1); Tb = mixedBoundary.domainT(2);
    Ba = mixedBoundary.domainB(1); Bb = mixedBoundary.domainB(2);


    Index.mixedB = Index.all(ny,grid.x<=Bb & grid.x>=Ba);
    Index.mixedB = Index.mixedB(:);
    Index.mixedB1 = -nx*ny;
    Index.mixedB2 = -nx*ny;
    if ~isempty(Index.mixedB)
        Index.mixedB1 = min(Index.mixedB);
        Index.mixedB2 = max(Index.mixedB);
    end
    
    
    Index.mixedT = Index.all(1,grid.x<=Tb & grid.x>=Ta);
    Index.mixedT = Index.mixedT(:);
    Index.mixedT1 = -nx*ny;
    Index.mixedT2 = -nx*ny;
    if ~isempty(Index.mixedT)
        Index.mixedT1 = min(Index.mixedT);
        Index.mixedT2 = max(Index.mixedT);
    end
    
    Index.mixedR = Index.all(grid.y<=Rb & grid.y>=Ra,nx);
    Index.mixedR = Index.mixedR(:);
    Index.mixedR1 = -nx*ny;
    Index.mixedR2 = -nx*ny;
    if ~isempty(Index.mixedR)
        Index.mixedR1 = min(Index.mixedR);
        Index.mixedR2 = max(Index.mixedR);
    end
    
    Index.mixedL = Index.all(grid.y<=Lb & grid.y>=La,1); 
    Index.mixedL= Index.mixedL(:); 
    Index.mixedL1 = -nx*ny;
    Index.mixedL2 = -nx*ny;
    if ~isempty(Index.mixedL)
        Index.mixedL1 = min(Index.mixedL);
        Index.mixedL2 = max(Index.mixedL);
    end
    
    Index.mixedBoundary = [Index.mixedB;Index.mixedT;Index.mixedR;Index.mixedL];
    
    Index.mixedBoundaryEnd = [Index.mixedB1,1;Index.mixedB2,2;Index.mixedT1,3;Index.mixedT2,4;...
        Index.mixedR1,5;Index.mixedR2,6;Index.mixedL1,7;Index.mixedL2,8];
    
    Index.mixedBoundaryEndAdjacent = [Index.mixedB1-ny,1;Index.mixedB2+ny,2;Index.mixedT1-ny,3;Index.mixedT2+ny,4;...
        Index.mixedR1-1,5;Index.mixedR2+1,6;Index.mixedL1-1,7;Index.mixedL2+1,8];

    % First line in the domain adjacent to the edge with mixed bc (not including corners)
    Index.mixedInB1 = Index.all(ny-1,grid.x<=Bb & grid.x>=Ba);   Index.mixedInB1=Index.mixedInB1(:); 
    Index.mixedInT1 = Index.all(2,grid.x<=Tb & grid.x>=Ta);Index.mixedInT1=Index.mixedInT1(:);
    Index.mixedInR1 = Index.all(grid.y<=Rb & grid.y>=Ra,nx-1);Index.mixedInR1=Index.mixedInR1(:);
    Index.mixedInL1 = Index.all(grid.y<=Lb & grid.y>=La,2);   Index.mixedInL1=Index.mixedInL1(:);
    
    Index.mixedInBoundary = [Index.mixedInB1;Index.mixedInT1;Index.mixedInR1;Index.mixedInL1];
    
    Index.mixedInBoundaryIdx = [Index.mixedInB1,ones(size(Index.mixedInB1));...
        Index.mixedInT1,2*ones(size(Index.mixedInT1));...
        Index.mixedInR1,3*ones(size(Index.mixedInR1));...
        Index.mixedInL1,4*ones(size(Index.mixedInL1))];

end