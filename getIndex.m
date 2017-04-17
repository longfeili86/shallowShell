function Index = getIndex(nx,ny)
% This function get the indices for the boudanry and ghost line ect.
% Input:
%    nx: number of grids in x direction       
%    ny: number of grids in y direction
% Output:
%    Index: an object that holds index info

    Index.all = reshape(1:(nx+4)*(ny+4),ny+4,nx+4);

    Index.interiorBoundary=Index.all(3:ny+2,3:nx+2);
    Index.interiorBoundary=Index.interiorBoundary(:);
    
    % Indices of left,right, top and bottom boundary of a rectangular domain
    Index.BoundaryB = Index.all(3,3:nx+2);   Index.BoundaryB=Index.BoundaryB(:); 
    Index.BoundaryT = Index.all(ny+2,3:nx+2);Index.BoundaryT=Index.BoundaryT(:);
    Index.BoundaryR = Index.all(3:ny+2,nx+2);Index.BoundaryR=Index.BoundaryR(:);
    Index.BoundaryL = Index.all(3:ny+2,3);   Index.BoundaryL=Index.BoundaryL(:);
    Index.Boundary = [Index.BoundaryL;Index.BoundaryR;Index.BoundaryB;Index.BoundaryT]; % column vector
    
    % First ghost line
    Index.GhostB1 = Index.all(2,3:nx+2);   Index.GhostB1=Index.GhostB1(:); 
    Index.GhostT1 = Index.all(ny+3,3:nx+2);Index.GhostT1=Index.GhostT1(:);
    Index.GhostR1 = Index.all(3:ny+2,nx+3);Index.GhostR1=Index.GhostR1(:);
    Index.GhostL1 = Index.all(3:ny+2,2);   Index.GhostL1=Index.GhostL1(:);
    Index.Ghost1 = [Index.GhostL1;Index.GhostR1;Index.GhostB1;Index.GhostT1];
    
    
    % Second ghost line
    Index.GhostB2 = Index.all(1,3:nx+2);   Index.GhostB2=Index.GhostB2(:); 
    Index.GhostT2 = Index.all(ny+4,3:nx+2);Index.GhostT2=Index.GhostT2(:);
    Index.GhostR2 = Index.all(3:ny+2,nx+4);Index.GhostR2=Index.GhostR2(:);
    Index.GhostL2 = Index.all(3:ny+2,1);   Index.GhostL2=Index.GhostL2(:);
    Index.Ghost2 = [Index.GhostL2;Index.GhostR2;Index.GhostB2;Index.GhostT2];

    %First line in that are needed to do finite difference on the boundary
    Index.GhostBin1 = Index.all(4,3:nx+2);   Index.GhostBin1=Index.GhostBin1(:); 
    Index.GhostTin1 = Index.all(ny+1,3:nx+2);Index.GhostTin1=Index.GhostTin1(:);
    Index.GhostRin1 = Index.all(3:ny+2,nx+1);Index.GhostRin1=Index.GhostRin1(:);
    Index.GhostLin1 = Index.all(3:ny+2,4);   Index.GhostLin1=Index.GhostLin1(:);
    Index.Ghostin1 = [Index.GhostLin1;Index.GhostRin1;Index.GhostBin1;Index.GhostTin1];

    
    
    Index.GhostCornersLL = Index.all(1:2,1:2);            % Lower Left Ghost Corners 
    Index.GhostCornersLR = Index.all(1:2,nx+3:nx+4);      % Lower Right Ghost Corners 
    Index.GhostCornersUL = Index.all(ny+3:ny+4,1:2);      % Upper Right Ghost Corners
    Index.GhostCornersUR = Index.all(ny+3:ny+4,nx+3:nx+4);% Upper Right Ghost Corners
    % all the ghost corners
    Index.UnusedGhostCorners = [Index.GhostCornersLL(:);Index.GhostCornersUL(:);...
        Index.GhostCornersLR(:);Index.GhostCornersUR(:)];  
    
    % these are the ghost corners needed by the corner bc W_xy=0 and the
    % stencil for \nabla^4 at the corners
    Index.GhostCornersLL1 = Index.all(2,2);           
    Index.GhostCornersLR1 = Index.all(2,nx+3);  
    Index.GhostCornersUL1 = Index.all(ny+3,2);  
    Index.GhostCornersUR1 = Index.all(ny+3,nx+3);
    Index.UsedGhostCorners = [Index.GhostCornersLL1;Index.GhostCornersUL1;...
        Index.GhostCornersLR1;Index.GhostCornersUR1];
    
    % these are the unused points   
    Index.UnusedGhostCorners = setdiff(Index.UnusedGhostCorners,Index.UsedGhostCorners);
    
    Index.CornersLL = Index.all(3,3);           
    Index.CornersLR = Index.all(3,nx+2);  
    Index.CornersUL = Index.all(ny+2,3);  
    Index.CornersUR = Index.all(ny+2,nx+2);
    Index.Corners = [Index.CornersLL;Index.CornersUL;...
        Index.CornersLR;Index.CornersUR]; 
    
   
    
    Index.GhostAll = setdiff(Index.all(:),Index.interiorBoundary);
    Index.UsedPoints = setdiff(Index.all(:),Index.UnusedGhostCorners);

     
end