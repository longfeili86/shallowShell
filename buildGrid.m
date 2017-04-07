function grid = buildGrid(domain, nx,ny,mixedBoundary)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input:
%  domain = [xa,xb,ya,yb]
%  nx: number of grids in x direction (exclude ghost)
%  ny: number of grids in y direction (exclude ghost)
%Output:
%  grid: object that holds grid information
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xa = domain(1);
xb = domain(2);
ya = domain(3);
yb = domain(4);

grid.x = linspace(xa,xb,nx);
grid.y = linspace(ya,yb,ny);

grid.hx = grid.x(2)-grid.x(1);
grid.hy = grid.y(2)-grid.y(1);
% Add two ghost lines
% grid.x = [xa-2*grid.hx,xa-grid.hx,grid.x,xb+grid.hx,xb+2*grid.hx];
% grid.y = [ya-2*grid.hy,ya-grid.hy,grid.y,yb+grid.hy,yb+2*grid.hy];

[grid.XX, grid.YY] = meshgrid(grid.x,grid.y);

% [grid.YY, grid.XX] = meshgrid(grid.y,grid.x);

% generate grid to impose mixed boundary conditions
La = mixedBoundary.domainL(1); Lb = mixedBoundary.domainL(2);
Ra = mixedBoundary.domainR(1); Rb = mixedBoundary.domainR(2);
Ta = mixedBoundary.domainT(1); Tb = mixedBoundary.domainT(2);
Ba = mixedBoundary.domainB(1); Bb = mixedBoundary.domainB(2);

grid.mixedLx = xa;
grid.mixedLy = grid.y(grid.y<Lb & grid.y>La);

grid.mixedRx = xb;
grid.mixedRy = grid.y(grid.y<Rb & grid.y>Ra);

grid.mixedBx = grid.x(grid.x<Bb & grid.x>Ba);
grid.mixedBy = ya;

grid.mixedTx = grid.x(grid.x<Tb & grid.x>Ta);
grid.mixedTy = yb;


end