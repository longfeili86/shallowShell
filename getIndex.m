function Index = getIndex(nx,ny)
% This function get the indices for the boudanry and ghost line ect.
% Input:
%    nx: number of grids in x direction       
%    ny: number of grids in y direction
% Output:
%    Index: an object that holds index info

    Index.all = reshape(1:(nx+4)*(ny+4),nx+4,ny+4);

    Index.interiorBoundary=Index.all(3:nx+2,3:ny+2);
    Index.interiorBoundary=Index.interiorBoundary(:);
    
    % Indices of left,right, top and bottom boundary of a rectangular domain
    Index.BoundaryB = Index.all(3,3:ny+2);   Index.BoundaryB=Index.BoundaryB(:); 
    Index.BoundaryT = Index.all(nx+2,3:ny+2);Index.BoundaryT=Index.BoundaryT(:);
    Index.BoundaryR = Index.all(3:nx+2,ny+2);Index.BoundaryR=Index.BoundaryR(:);
    Index.BoundaryL = Index.all(3:nx+2,3);   Index.BoundaryL=Index.BoundaryL(:);
    Index.Boundary = [Index.BoundaryL;Index.BoundaryR;Index.BoundaryB;Index.BoundaryT]; % column vector
    
    % First ghost line
    Index.GhostB1 = Index.all(2,3:ny+2);   Index.GhostB1=Index.GhostB1(:); 
    Index.GhostT1 = Index.all(nx+3,3:ny+2);Index.GhostT1=Index.GhostT1(:);
    Index.GhostR1 = Index.all(3:nx+2,ny+3);Index.GhostR1=Index.GhostR1(:);
    Index.GhostL1 = Index.all(3:nx+2,2);   Index.GhostL1=Index.GhostL1(:);
    Index.Ghost1 = [Index.GhostL1;Index.GhostR1;Index.GhostB1;Index.GhostT1];
    
    
    % Second ghost line
    Index.GhostB2 = Index.all(1,3:ny+2);   Index.GhostB2=Index.GhostB2(:); 
    Index.GhostT2 = Index.all(nx+4,3:ny+2);Index.GhostT2=Index.GhostT2(:);
    Index.GhostR2 = Index.all(3:nx+2,ny+4);Index.GhostR2=Index.GhostR2(:);
    Index.GhostL2 = Index.all(3:nx+2,1);   Index.GhostL2=Index.GhostL2(:);
    Index.Ghost2 = [Index.GhostL2;Index.GhostR2;Index.GhostB2;Index.GhostT2];

    %First line in that are needed to do finite difference on the boundary
    Index.GhostBin1 = Index.all(4,3:ny+2);   Index.GhostBin1=Index.GhostBin1(:); 
    Index.GhostTin1 = Index.all(nx+1,3:ny+2);Index.GhostTin1=Index.GhostTin1(:);
    Index.GhostRin1 = Index.all(3:nx+2,ny+1);Index.GhostRin1=Index.GhostRin1(:);
    Index.GhostLin1 = Index.all(3:nx+2,4);   Index.GhostLin1=Index.GhostLin1(:);
    Index.Ghostin1 = [Index.GhostLin1;Index.GhostRin1;Index.GhostBin1;Index.GhostTin1];

    
    
    Index.GhostCornersLL = Index.all(1:2,1:2);            % Lower Left Ghost Corners 
    Index.GhostCornersUL = Index.all(1:2,ny+3:ny+4);      % Upper Left Ghost Corners 
    Index.GhostCornersLR = Index.all(nx+3:nx+4,1:2);      % Lower Right Ghost Corners
    Index.GhostCornersUR = Index.all(nx+3:nx+4,ny+3:ny+4);% Upper Right Ghost Corners
    
    
    % these are the unused points to be excluded from our system
    Index.GhostCorners = [Index.GhostCornersLL(:);Index.GhostCornersUL(:);...
        Index.GhostCornersLR(:);Index.GhostCornersUR(:)]; 
    
    Index.GhostAll = setdiff(Index.all(:),Index.interiorBoundary);
    Index.UsedPoints = setdiff(Index.all(:),Index.GhostCorners);

     
end