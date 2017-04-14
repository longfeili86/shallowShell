function grid = buildGrid(domain, nx,ny)
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


assert((xa<xb)&&(ya<yb),'Error: invalid domain. make sure xa<xb and ya<yb\n');

grid.x = linspace(xa,xb,nx);
grid.y = linspace(ya,yb,ny);

grid.hx = grid.x(2)-grid.x(1);
grid.hy = grid.y(2)-grid.y(1);
% Add two ghost lines
grid.x = [xa-2*grid.hx,xa-grid.hx,grid.x,xb+grid.hx,xb+2*grid.hx];
grid.y = [ya-2*grid.hy,ya-grid.hy,grid.y,yb+grid.hy,yb+2*grid.hy];

[grid.XX, grid.YY] = meshgrid(grid.x,grid.y);

end