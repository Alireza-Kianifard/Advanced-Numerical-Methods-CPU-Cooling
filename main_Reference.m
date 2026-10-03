clear; clc; close all;

%% 1. Parameters Definition
% Geometry
Lx          = 0.02;             % Width (m)
Ly          = 0.01;             % Height (m)

% Material Properties (Aluminum)
k           = 237;              % Thermal conductivity (W/m.K)
rho         = 2700;             % Density (kg/m^3)
cp          = 900;              % Specific heat (J/kg.K)
alpha       = k / (rho * cp);   % Thermal diffusivity (m^2/s) - Unused directly by PDE tool

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

% Time Settings
t_final     = 50;               % Total simulation time (s)
dt          = 0.1;              % Time step for output
tlist       = 0:dt:t_final;

%% 2. Create PDE Model
disp('Setting up PDE Model');
thermalmodel = createpde('thermal', 'transient');

% Create Geometry (Rectangle)
geom = decsg([3; 4; 0; Lx; Lx; 0; 0; 0; Ly; Ly]);
geometryFromEdges(thermalmodel, geom);

% Assign Material Properties
thermalProperties(thermalmodel, 'ThermalConductivity', k, 'MassDensity', rho, 'SpecificHeat', cp);

%% 3. Apply Boundary & Initial Conditions
% Bottom Edge
thermalBC(thermalmodel, 'Edge', 1, 'HeatFlux', q_flux);

% Left & Right Edges
thermalBC(thermalmodel, 'Edge', [2, 4], 'ConvectionCoefficient', h_side, 'AmbientTemperature', T_inf);

% Top Edge (Smart Convection)
thermalBC(thermalmodel, 'Edge', 3, 'ConvectionCoefficient', @(region, state) smartCoolingFunc(region, state, h_base, T_target, Kp, h_max_limit), 'AmbientTemperature', T_inf);

% Initial Condition
thermalIC(thermalmodel, T_init);

%% 4. Mesh Generation & Solving
generateMesh(thermalmodel, 'Hmax', min(Lx/40, Ly/20));

disp('Running MATLAB PDE Toolbox Solver');
tic;
result = solve(thermalmodel, tlist);
time_PDE = toc;
disp(['PDE Toolbox completed in ', num2str(time_PDE, '%.2f'), ' seconds.']);

%% 5. Post-Processing
T_all      = result.Temperature; 
hist_T_max = max(T_all, [], 1);

% Calculate smart fan history
tol        = 1e-9;
topNodes   = find(abs(thermalmodel.Mesh.Nodes(2,:) - Ly) < tol);
T_top_avg  = mean(T_all(topNodes, :), 1); 
hist_h_top = min(h_base + Kp .* max(0, T_top_avg - T_target), h_max_limit); 

%% 6. Plotting Results
disp('Generating Plots');

% Figure 1: Max Temperature vs Time
figure('Name', 'Max Temperature (Reference)', 'Color', 'w', 'Position', [100, 100, 800, 500]);
plot(tlist, hist_T_max, 'r-', 'LineWidth', 2);
hold on;
yline(90, 'k--', 'Safe Limit (90^\circC)', 'LineWidth', 1.5, 'LabelHorizontalAlignment', 'left');
title('Maximum Temperature over Time (Reference PDE)');
xlabel('Time (s)');
ylabel('Max Temperature (^\circC)');
grid on;

% Figure 2: Smart Fan Behavior vs Time
figure('Name', 'Smart Fan Behavior (Reference)', 'Color', 'w', 'Position', [150, 150, 800, 500]);
plot(tlist, hist_h_top, 'b-', 'LineWidth', 2);
title('Smart Fan Operation (h_{top} vs Time) - Reference PDE');
xlabel('Time (s)');
ylabel('h_{top} (W/m^2.K)');
grid on;

% Figure 3: Final Temperature Distribution
figure('Name', 'Temperature Distribution (Reference)', 'Color', 'w', 'Position', [200, 200, 600, 400]);
pdeplot(thermalmodel, 'XYData', T_all(:, end), 'Contour', 'on', 'ColorMap', 'jet');
title(sprintf('Temperature Dist. at t = %.1f s (Reference PDE)', t_final));
xlabel('x (m)');
ylabel('y (m)');
axis equal; axis tight;

%% Local Function for Smart Boundary Condition
function h = smartCoolingFunc(region, state, h_base, T_target, Kp, h_max_limit)
    if any(isnan(state.u))
        h = NaN(1, numel(region.x));
        return;
    end
    h_current = h_base + Kp * max(0, state.u - T_target);
    h = min(h_current, h_max_limit);
end
