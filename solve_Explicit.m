function [T, X, Y, hist_time, hist_T_max, hist_h_top] = solve_Explicit(params)
    % Unpack parameters
    Nx          = params.Nx;        Ny          = params.Ny;
    dx          = params.dx;        dy          = params.dy;
    alpha       = params.alpha;     dt          = params.dt; 
    Nt          = params.Nt;        k           = params.k; 
    q_flux      = params.q_flux;    h_side      = params.h_side;
    T_inf       = params.T_inf;     T_init      = params.T_init;
    h_base      = params.h_base;    T_target    = params.T_target;
    Kp          = params.Kp;        h_max_limit = params.h_max_limit;
    
    % Domain setup
    x = linspace(0, dx*(Nx-1), Nx);
    y = linspace(0, dy*(Ny-1), Ny);
    [X, Y] = meshgrid(x, y);
    
    % Initialize temperature matrix
    T     = T_init * ones(Nx, Ny);
    T_new = T;
    
    % Pre-compute multipliers
    rx = alpha * dt / dx^2;
    ry = alpha * dt / dy^2;
    
    % Allocate history arrays
    hist_time  = zeros(1, Nt);
    hist_T_max = zeros(1, Nt);
    hist_h_top = zeros(1, Nt);
    
    % Index vectors for vectorization
    i = 2:Nx-1;
    j = 2:Ny-1;
    
    % Time integration loop
    for n = 1:Nt
        % Smart cooling controller
        T_top_avg     = mean(T(:, Ny));
        h_top_current = min(h_base + Kp * max(0, T_top_avg - T_target), h_max_limit);
        
        % Record history
        hist_time(n)  = n * dt;
        hist_T_max(n) = max(T(:));
        hist_h_top(n) = h_top_current;
        
        % 1. Interior nodes (Vectorized)
        T_new(i,j) = T(i,j) + rx*(T(i+1,j) - 2*T(i,j) + T(i-1,j)) ...
                            + ry*(T(i,j+1) - 2*T(i,j) + T(i,j-1));
        
        % 2. Boundaries (Vectorized)
        % Bottom (Constant Heat Flux)
        T_new(i,1)  = T(i,1)  + rx*(T(i+1,1) - 2*T(i,1) + T(i-1,1)) ...
                              + ry*(2*T(i,2) - 2*T(i,1) + (2*dy*q_flux)/k);
        % Top (Smart Convection)
        T_new(i,Ny) = T(i,Ny) + rx*(T(i+1,Ny) - 2*T(i,Ny) + T(i-1,Ny)) ...
                              + ry*(2*T(i,Ny-1) - 2*T(i,Ny) - (2*dy*h_top_current/k).*(T(i,Ny) - T_inf));
        
        % Left (Natural Convection)
        T_new(1,j)  = T(1,j)  + rx*(2*T(2,j) - 2*T(1,j) - (2*dx*h_side/k).*(T(1,j) - T_inf)) ...
                              + ry*(T(1,j+1) - 2*T(1,j) + T(1,j-1));
        % Right (Natural Convection)
        T_new(Nx,j) = T(Nx,j) + rx*(2*T(Nx-1,j) - 2*T(Nx,j) - (2*dx*h_side/k).*(T(Nx,j) - T_inf)) ...
                              + ry*(T(Nx,j+1) - 2*T(Nx,j) + T(Nx,j-1));
        
        % 3. Corners
        T_new(1,1)   = T(1,1)   + rx*(2*T(2,1) - 2*T(1,1) - (2*dx*h_side/k)*(T(1,1) - T_inf)) + ry*(2*T(1,2) - 2*T(1,1) + (2*dy*q_flux)/k);
        T_new(Nx,1)  = T(Nx,1)  + rx*(2*T(Nx-1,1) - 2*T(Nx,1) - (2*dx*h_side/k)*(T(Nx,1) - T_inf)) + ry*(2*T(Nx,2) - 2*T(Nx,1) + (2*dy*q_flux)/k);
        T_new(1,Ny)  = T(1,Ny)  + rx*(2*T(2,Ny) - 2*T(1,Ny) - (2*dx*h_side/k)*(T(1,Ny) - T_inf)) + ry*(2*T(1,Ny-1) - 2*T(1,Ny) - (2*dy*h_top_current/k)*(T(1,Ny) - T_inf));
        T_new(Nx,Ny) = T(Nx,Ny) + rx*(2*T(Nx-1,Ny) - 2*T(Nx,Ny) - (2*dx*h_side/k)*(T(Nx,Ny) - T_inf)) + ry*(2*T(Nx,Ny-1) - 2*T(Nx,Ny) - (2*dy*h_top_current/k)*(T(Nx,Ny) - T_inf));
        
        % Update field
        T = T_new;
    end
end
