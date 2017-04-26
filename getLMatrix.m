function Lmtx = getLMatrix(diffMtx,u)
% Calculate the discretizated matrix for L[u,~]
% L[u,v] = u_xx * v_yy + u_yy * v_xx - 2 * u_xy * v_xy;
% Boundary points are included
%
% -- Longfei Li

[r,c]=size(diffMtx.Dxx);        
Lmtx = spdiags(diffMtx.Dxx*u,0,r,c)*(diffMtx.Dyy )+...
    spdiags(diffMtx.Dyy * u,0,r,c)*(diffMtx.Dxx )-...
    2.0*spdiags(diffMtx.Dxy*u,0,r,c)*(diffMtx.Dxy);

end
