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

% Time Settings (ADI is unconditionally stable)
dt          = 0.1;              % Chosen time step
t_final     = 50;               % Total simulation time (s)
time_steps  = 0:dt:t_final;
Nt          = length(time_steps) - 1;

% Grid for plotting
x           = linspace(0, Lx, Nx);
y           = linspace(0, Ly, Ny);
[X, Y]      = meshgrid(x, y);

%% 2. Initialization & Pre-allocation
T = T_init * ones(Nx, Ny);

% History arrays
hist_time  = zeros(1, Nt+1);
hist_T_max = zeros(1, Nt+1);
hist_h_top = zeros(1, Nt+1);

hist_time(1)  = 0;
hist_T_max(1) = T_init;
hist_h_top(1) = h_base;

%% 3. Solve using ADI Method
disp('Running ADI Solver');
tic;
for n = 1:Nt
    % Execute one step
    T = solve_ADI(T, Nx, Ny, dx, dy, dt, alpha, k, q_flux, h_side, h_base, T_inf, Kp, T_target, h_max_limit);
    
    % Update top convection
    T_top_avg = mean(T(:, Ny));
    h_top_current = min(h_base + Kp * max(0, T_top_avg - T_target), h_max_limit);
    
    % Record history
    hist_time(n+1)  = n * dt;
    hist_T_max(n+1) = max(T(:));
    hist_h_top(n+1) = h_top_current;
end
time_ADI = toc;
disp(['ADI completed in ', num2str(time_ADI, '%.2f'), ' seconds.']);

%% 4. Post-Processing
disp('Generating Plots');

% Figure 1: Max Temperature vs Time
figure('Name', 'Max Temperature (ADI)', 'Color', 'w', 'Position', [100, 100, 800, 500]);
plot(hist_time, hist_T_max, 'r-', 'LineWidth', 2);
hold on;
yline(90, 'k--', 'Safe Limit (90^\circC)', 'LineWidth', 1.5, 'LabelHorizontalAlignment', 'left');
title('Maximum Temperature over Time (ADI Method)');
xlabel('Time (s)');
ylabel('Max Temperature (^\circC)');
grid on;

% Figure 2: Smart Fan Behavior vs Time
figure('Name', 'Smart Fan Behavior (ADI)', 'Color', 'w', 'Position', [150, 150, 800, 500]);
plot(hist_time, hist_h_top, 'b-', 'LineWidth', 2);
title('Smart Fan Operation (h_{top} vs Time) - ADI Method');
xlabel('Time (s)');
ylabel('h_{top} (W/m^2.K)');
grid on;

% Figure 3: Final Temperature Distribution
figure('Name', 'Temperature Distribution (ADI)', 'Color', 'w', 'Position', [200, 200, 600, 400]);
contourf(X', Y', T, 30, 'LineColor', 'none');
colormap('jet');
c = colorbar;
c.Label.String = 'Temperature (^\circC)';
title(sprintf('Temperature Dist. at t = %.1f s (ADI Method)', t_final));
xlabel('x (m)');
ylabel('y (m)');
axis equal; axis tight;
