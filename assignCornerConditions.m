function [A,RHS]=assignCornerConditions(A,RHS,Index,diffMtx,bcType)
   
    if(bcType==3) % assign the corner conditions W_xy=0  for free bc
        A(Index.UsedGhostCorners,:)=0.; % zero out existing coefficients
        A(Index.UsedGhostCorners,:)=diffMtx.Dxy(Index.Corners,:) ;   
        RHS(Index.UsedGhostCorners)=0.;
    else
        return
    end
    
end