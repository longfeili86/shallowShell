function L = getLOperator(diffMtx,u,v)
% Calculate the discretizated bilinear operator
% L[u,v] = u_xx * v_yy + u_yy * v_xx - 2 * u_xy * v_xy;
% Boundary points are included

L = (diffMtx.Dxx * u) .* (diffMtx.Dyy * v) + ...
    (diffMtx.Dyy * u) .* (diffMtx.Dxx * v) - ...
    2.0* (diffMtx.Dxy*u) .*(diffMtx.Dxy*v);

end
