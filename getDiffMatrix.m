function Mtx = getDiffMatrix(nx,ny,hx,hy)
% Calculate the differentiation matrices using finite difference schemes
    % Define differetiation Matrices
    % nx*ny by nx*ny matrix
    % Boundary points are included
    D0x = spdiags(ones(nx,1),1,nx,nx)- spdiags(ones(nx,1),-1,nx,nx); % centered fd 
    Dpx = spdiags(ones(nx,1),1,nx,nx)- spdiags(ones(nx,1),0,nx,nx);  % foward fd
    Dmx = spdiags(ones(nx,1),0,nx,nx)- spdiags(ones(nx,1),-1,nx,nx); % backward fd
    D0y = spdiags(ones(ny,1),1,ny,ny)- spdiags(ones(ny,1),-1,ny,ny); % centered fd 
    Dpy = spdiags(ones(ny,1),1,ny,ny)- spdiags(ones(ny,1),0,ny,ny);  % backward fd
    Dmy = spdiags(ones(ny,1),0,ny,ny)- spdiags(ones(ny,1),-1,ny,ny); % centered fd
    D2 = spdiags(ones(nx,1),1,nx,nx)-2*spdiags(ones(nx,1),0,nx,nx)...
        + spdiags(ones(nx,1),-1,nx,nx); % centered fd for second x derivative
    D1 = spdiags(ones(nx,1),1,nx,nx)-spdiags(ones(nx,1),-1,nx,nx);% centered fd 
    % vectorize the grids by columns
%     Mtx.D2 = D2/(hx^2);
    
    Mtx.D2 = D2;
    
    %Mtx.D1 = D1/(2*hx);
    
    Mtx.D1 = D1/(2);
    
%     Mtx.D0x = kron(D0x,speye(ny))/(2*hx);
    Mtx.D0x = kron(D0x,speye(ny))/2;
%     Mtx.Dpx = kron(Dpx,speye(ny))/hx;
    Mtx.Dpx = kron(Dpx,speye(ny));
%     Mtx.Dmx = kron(Dmx,speye(ny))/hx;
    Mtx.Dmx = kron(Dmx,speye(ny));
%     Mtx.D0y = kron(speye(nx),D0y)/(2*hy);
    Mtx.D0y = kron(speye(nx),D0y)/2;
%     Mtx.Dpy = kron(speye(nx),Dpy)/hy;
    Mtx.Dpy = kron(speye(nx),Dpy);
%     Mtx.Dmy = kron(speye(nx),Dmy)/hy;
    Mtx.Dmy = kron(speye(nx),Dmy);
    Mtx.Dxx = Mtx.Dpx*Mtx.Dmx;
    Mtx.Dyy = Mtx.Dpy*Mtx.Dmy;
    
    % ?? correct??
    Mtx.Dpxx = Mtx.Dpx*Mtx.Dpx;
    Mtx.Dmxx = Mtx.Dmx*Mtx.Dmx;
    Mtx.Dpyy = Mtx.Dpy*Mtx.Dpy;
    Mtx.Dmyy = Mtx.Dmy*Mtx.Dmy;
    
    Mtx.Dpxxx = Mtx.Dpx*Mtx.Dpx*Mtx.Dpx;
    Mtx.Dmxxx = Mtx.Dmx*Mtx.Dmx*Mtx.Dpx;
    Mtx.Dpyyy = Mtx.Dpy*Mtx.Dpy*Mtx.Dpy;
    Mtx.Dmyyy = Mtx.Dmy*Mtx.Dmy*Mtx.Dmy;
    

    Mtx.Dpxpy = Mtx.Dpx*Mtx.Dpy;
    Mtx.Dpxmy = Mtx.Dpx*Mtx.Dmy;
    Mtx.Dmxpy = Mtx.Dmx*Mtx.Dpy;
    Mtx.Dmxmy = Mtx.Dmx*Mtx.Dmy;
    
    Mtx.Dxy = Mtx.D0x*Mtx.D0y;     

    Mtx.Dxxx= Mtx.D0x*Mtx.Dxx;
    Mtx.Dyyy=Mtx.D0y*Mtx.Dyy;
    Mtx.Dyyx=Mtx.D0x*Mtx.Dyy;
    Mtx.Dxxy=Mtx.D0y*Mtx.Dxx;
    Mtx.Dxxxx = Mtx.Dxx*Mtx.Dxx;
    Mtx.Dyyyy = Mtx.Dyy*Mtx.Dyy;
    Mtx.Dxxyy = Mtx.Dxx*Mtx.Dyy;
    Mtx.Dh = Mtx.Dpx*Mtx.Dmx+Mtx.Dpy*Mtx.Dmy; % laplace operator
    Mtx.BiDh = Mtx.Dxx*Mtx.Dxx+2*Mtx.Dxx*Mtx.Dyy+Mtx.Dyy*Mtx.Dyy; % Biharmonic operator
    
end
