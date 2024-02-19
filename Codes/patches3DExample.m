clear; clc

vidName = 'video3D.mp4';
% Parámetros del dominio y condiciones iniciales
domain = [-200, -100,-100,200,100,100,100,50,50]; % [x_min, y_min, z_min, x_max, y_max, z_max, Nx, Ny, Nz]
Lx = domain(4) - domain(1);
Ly = domain(5) - domain(2);
Lz = domain(6) - domain(3);

tSpan = [0, 0.2, 185]; % [t_inicial, dt, t_final]

params.IN   = 0.00225;
params.mu   = 0.5;
params.k    = 0.01;
params.mN   = 0.015;
params.fP   = 1.8;
params.HN   = 0.005;
params.HP   = 4.0;
params.d    = 0.04;
params.vmax = 0.3;
params.A    = 0.5;

nEddy = 100;
eddy.x = domain(4)*rand(nEddy, 1) - domain(1);  % x_i
eddy.y = domain(5)*rand(nEddy, 1) - domain(2);  % y_i
eddy.z = domain(6)*rand(nEddy, 1) - domain(3);  % z_i
eddy.s = 2*(round(rand(nEddy, 1))-0.5);         % sigma_i
eddy.r = abs(normrnd(20,5,[nEddy,1]));          % r_i

% Condiciones iniciales
[x,y,z] = meshgrid( linspace(domain(1), domain(4), domain(7)), ...
                    linspace(domain(2), domain(5), domain(8)), ...
                    linspace(domain(3), domain(6), domain(9)) );

iCondN = 1 + params.A * sin(pi/domain(4)*2*x) .* sin(pi/domain(5)*2*y) .* sin(pi/domain(6)*2*z);
iCondP = 1 + params.A * sin(pi/domain(4)*2*x) .* sin(pi/domain(5)*2*y) .* sin(pi/domain(6)*2*z);

iCondN = 1 + params.A*sin(pi/Lx*2*x);
iCondP = 1 + params.A*sin(pi/Lx*2*x);

% Simulación
[N, P, T] = patches3D(domain, tSpan, params, eddy, iCondN, iCondP,'adjust');

%==========================================================================
% Plots
%==========================================================================

% Color map
cmapC = [0.023529411764706   0.070588235294118   0.388235294117647;
         0.509803921568627   1.000000000000000   0.509803921568627;
         0.933333333333333   1.000000000000000   0.349019607843137];
cmapP = [0,0.9,1];
cmap = interp1(cmapP,cmapC,(1:255)/255);


