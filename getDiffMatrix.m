function Mtx = getDiffMatrix(nx,ny,hx,hy)
% Calculate the differentiation matrices using finite difference schemes
    % Define differetiation Matrices
    % (nx+4)*(ny+4) by (nx+4)*(ny+4) matrices that includes 
    % all Boundary and ghost points
    
    % 1d differentiation matrices
    D0x = spdiags(ones(nx+4,1),1,nx+4,nx+4)- spdiags(ones(nx+4,1),-1,nx+4,nx+4); % centered fd 
    Dpx = spdiags(ones(nx+4,1),1,nx+4,nx+4)- spdiags(ones(nx+4,1),0,nx+4,nx+4);  % foward fd
    Dmx = spdiags(ones(nx+4,1),0,nx+4,nx+4)- spdiags(ones(nx+4,1),-1,nx+4,nx+4); % backward fd
    D0y = spdiags(ones(ny+4,1),1,ny+4,ny+4)- spdiags(ones(ny+4,1),-1,ny+4,ny+4); % centered fd 
    Dpy = spdiags(ones(ny+4,1),1,ny+4,ny+4)- spdiags(ones(ny+4,1),0,ny+4,ny+4);  % forward fd
    Dmy = spdiags(ones(ny+4,1),0,ny+4,ny+4)- spdiags(ones(ny+4,1),-1,ny+4,ny+4); % backward fd

    % 2d differentiation matrices
    Mtx.D0x = kron(D0x,speye(ny+4))/(2*hx);
    Mtx.Dpx = kron(Dpx,speye(ny+4))/hx;
    Mtx.Dmx = kron(Dmx,speye(ny+4))/hx;
    Mtx.D0y = kron(speye(nx+4),D0y)/(2*hy);
    Mtx.Dpy = kron(speye(nx+4),Dpy)/hy;
    Mtx.Dmy = kron(speye(nx+4),Dmy)/hy;
    Mtx.Dxx = Mtx.Dpx*Mtx.Dmx;
    Mtx.Dyy = Mtx.Dpy*Mtx.Dmy;
    Mtx.Dxy = Mtx.D0x*Mtx.D0y; % this is needed for corner bc  
    Mtx.Dxxx= Mtx.D0x*Mtx.Dxx;
    Mtx.Dyyy=Mtx.D0y*Mtx.Dyy;
    Mtx.Dyyx=Mtx.D0x*Mtx.Dyy;
    Mtx.Dxxy=Mtx.D0y*Mtx.Dxx;
    Mtx.Dh = Mtx.Dpx*Mtx.Dmx+Mtx.Dpy*Mtx.Dmy; % laplace operator
    Mtx.BiDh = Mtx.Dxx*Mtx.Dxx+2*Mtx.Dxx*Mtx.Dyy+Mtx.Dyy*Mtx.Dyy; % Biharmonic operator
    Mtx.I = speye((nx+4)*(ny+4));
end
