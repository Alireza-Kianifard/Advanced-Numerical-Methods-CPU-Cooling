clear; clc; close all;

%% 1. Parameters Definition
% Geometry
Lx          = 0.02;             % Width (m)
Ly          = 0.01;             % Height (m)

% Material Properties (Aluminum)
k           = 237;              % Thermal conductivity (W/m.K)
rho         = 2700;             % Density (kg/m^3)
cp          = 900;              % Specific heat (J/kg.K)
alpha       = k / (rho * cp);   % Thermal diffusivity (m^2/s)

% Boundary & Initial Conditions
q_flux      = 200000;           % Bottom constant heat flux (W/m^2)
h_side      = 64;               % Side natural convection (W/m^2.K)
T_inf       = 25;               % Ambient temperature (C)
T_init      = 25;               % Initial temperature (C)

% Cooling Parameters (Smart Fan)
h_base      = 1000;             % Base convection coefficient (fan idle)
T_target    = 55;               % Temperature threshold to trigger cooling (C)
Kp          = 100;              % Proportional gain for controller
h_max_limit = 4000;             % Maximum capacity of the cooling fan

% Mesh Grid
Nx          = 81;
Ny          = 41;
dx          = Lx / (Nx - 1);
dy          = Ly / (Ny - 1);

% Time Settings (Explicit FTCS Stability)
dt_max_cond = 1 / (2 * alpha * (1/dx^2 + 1/dy^2));
dt_max_conv = (rho * cp * dy) / (2 * h_max_limit + 2 * k / dy); 
dt_max      = min(dt_max_cond, dt_max_conv);

dt          = 0.9 * dt_max;     % Safe time step
t_final     = 50;               % Total simulation time (s)
Nt          = ceil(t_final / dt);

% Pack parameters into struct
params.Nx = Nx; params.Ny = Ny; params.dx = dx; params.dy = dy;
params.alpha = alpha; params.dt = dt; params.Nt = Nt;
params.k = k; params.q_flux = q_flux; params.h_side = h_side; 
params.T_inf = T_inf; params.T_init = T_init;
params.h_base = h_base; params.T_target = T_target;
params.Kp = Kp; params.h_max_limit = h_max_limit;

%% 2. Solve using Explicit FTCS
disp('Running Explicit FTCS Solver');
tic;
[T_explicit, X, Y, hist_time, hist_T_max, hist_h_top] = solve_Explicit(params);
time_explicit = toc;
disp(['Explicit completed in ', num2str(time_explicit, '%.2f'), ' seconds.']);

%% 3. Post-Processing
disp('Generating Plots');

% Figure 1: Max Temperature vs Time
figure('Name', 'Max Temperature (Explicit)', 'Color', 'w', 'Position', [100, 100, 800, 500]);
plot(hist_time, hist_T_max, 'r-', 'LineWidth', 2);
hold on;
yline(90, 'k--', 'Safe Limit (90^\circC)', 'LineWidth', 1.5, 'LabelHorizontalAlignment', 'left');
title('Maximum Temperature over Time (Explicit FTCS)');
xlabel('Time (s)');
ylabel('Max Temperature (^\circC)');
grid on;

% Figure 2: Smart Fan Behavior vs Time
figure('Name', 'Smart Fan Behavior (Explicit)', 'Color', 'w', 'Position', [150, 150, 800, 500]);
plot(hist_time, hist_h_top, 'b-', 'LineWidth', 2);
title('Smart Fan Operation (h_{top} vs Time) - Explicit FTCS');
xlabel('Time (s)');
ylabel('h_{top} (W/m^2.K)');
grid on;

% Figure 3: Final Temperature Distribution
figure('Name', 'Temperature Distribution (Explicit)', 'Color', 'w', 'Position', [200, 200, 600, 400]);
contourf(X', Y', T_explicit, 30, 'LineColor', 'none');
colormap('jet');
c = colorbar;
c.Label.String = 'Temperature (^\circC)';
title(sprintf('Temperature Dist. at t = %.1f s (Explicit FTCS)', t_final));
xlabel('x (m)');
ylabel('y (m)');
axis equal; axis tight;
