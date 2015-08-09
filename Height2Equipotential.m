function H_eq = Height2Equipotential(filename,lat_grid,lon_grid,mu,Rref,lmcosi_g,T)

w = 2*pi./(3600*T);
step = (lon_grid(end) - lon_grid(end-1))*180/pi;

%% Getting random shape
Npts = 1000;
[x_rand,y_rand,z_rand]=ReadSPC(filename,Npts,'rand');
x_rand = x_rand * 1000;
y_rand = y_rand * 1000;
z_rand = z_rand * 1000;
[lon_rand, lat_rand, r_rand] = cart2sph(x_rand,y_rand,z_rand);
r_mean = mean(r_rand);
g_approx = mu/(r_mean^2);
U_g = GravityPotential(mu,Rref,lmcosi_g,x_rand,y_rand,z_rand);
U_r=RotPot(x_rand,y_rand,w);
U_tot = U_g + U_r;
U_mean = mean(U_tot);

%% Grid shape

grid = zeros(2,numel(lat_grid));
s    = size(lat_grid);

grid(1,:) = lon_grid(:);
grid(2,:) = lat_grid(:);

% get points from dsk
handle = cspice_dasopr( filename );
[dladsc, found] = cspice_dlabfs( handle );
[spoints, ~] = cspice_llgrid_pl02( handle, dladsc, grid );

x_grid = spoints(1,:);
y_grid = spoints(2,:);
z_grid = spoints(3,:);

x_grid = x_grid * 1000;
y_grid = y_grid * 1000;
z_grid = z_grid * 1000;

x_grid = reshape(x_grid,s);
y_grid = reshape(y_grid,s);
z_grid = reshape(z_grid,s);

r_grid = sqrt(x_grid.*x_grid+y_grid.*y_grid+z_grid.*z_grid);

U_g = GravityPotential(mu,Rref,lmcosi_g,x_grid,y_grid,z_grid);
U_r=RotPot(x_grid,y_grid,w);
U_tot = U_g + U_r;

dr   = (U_tot - U_mean)/g_approx;
r_eq = r_grid + dr;
[x_eq,y_eq,z_eq] = sph2cart(lon_grid,lat_grid,r_eq);

tol = 1;
eps = tol + 1;

while (eps > tol)
    U_g = GravityPotential(mu,Rref,lmcosi_g,x_eq,y_eq,z_eq);
    U_r=RotPot(x_eq,y_eq,w);
    U_tot = U_g + U_r;
    dr   = (U_tot - U_mean)/g_approx;
    r_eq = r_eq + dr;
    [x_eq,y_eq,z_eq] = sph2cart(lon_grid, lat_grid,r_eq);
    eps = max(abs(dr(:)));
end

H_eq = r_grid - r_eq;


