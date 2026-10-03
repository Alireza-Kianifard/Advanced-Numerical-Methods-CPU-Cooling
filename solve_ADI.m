function T_new = solve_ADI(T, Nx, Ny, dx, dy, dt, alpha, k, q_flux, h_side, h_base, T_inf, Kp, T_target, h_max_limit)
    rx = alpha * (dt / 2) / (dx^2);
    ry = alpha * (dt / 2) / (dy^2);
    
    T_half = zeros(Nx, Ny);
    T_new  = zeros(Nx, Ny);
    
    %% First Half-Step: Implicit in X, Explicit in Y
    T_top_avg_initial = mean(T(:, Ny));
    h_top_initial     = min(h_base + Kp * max(0, (T_top_avg_initial - T_target)), h_max_limit);
    
    % Build A_x ONCE (it is independent of j)
    A_x     = zeros(Nx, Nx);
    C_left  = 2 * dx * h_side / k;
    C_right = 2 * dx * h_side / k;
    
    A_x(1, 1) = 1 + 2*rx + rx*C_left;
    A_x(1, 2) = -2*rx;
    for i = 2:Nx-1
        A_x(i, i-1) = -rx;
        A_x(i, i)   = 1 + 2*rx;
        A_x(i, i+1) = -rx;
    end
    A_x(Nx, Nx-1) = -2*rx;
    A_x(Nx, Nx)   = 1 + 2*rx + rx*C_right;
    
    % Factorize for speed
    dA_x = decomposition(A_x);
    
    for j = 1:Ny
        B_x = zeros(Nx, 1);
        for i = 1:Nx
            if j == 1
                B_x(i) = ry*T(i,2) + (1 - 2*ry)*T(i,j) + ry*(T(i,2) + 2*dy*q_flux/k);
            elseif j == Ny
                C_top_exp  = 2 * dy * h_top_initial / k;
                T_up_ghost = T(i,Ny-1) + C_top_exp*(T_inf - T(i,Ny));
                B_x(i)     = ry*T(i,Ny-1) + (1 - 2*ry)*T(i,j) + ry*T_up_ghost;
            else
                B_x(i) = ry*T(i,j-1) + (1 - 2*ry)*T(i,j) + ry*T(i,j+1);
            end
        end
        B_x(1)  = B_x(1)  + rx * C_left * T_inf;
        B_x(Nx) = B_x(Nx) + rx * C_right * T_inf;
        
        T_half(:, j) = dA_x \ B_x;
    end
    
    %% Second Half-Step: Implicit in Y, Explicit in X
    T_top_avg_half = mean(T_half(:, Ny));
    h_top_half     = min(h_base + Kp * max(0, (T_top_avg_half - T_target)), h_max_limit);
    
    % Build A_y ONCE (it is independent of i)
    A_y = zeros(Ny, Ny);
    A_y(1, 1) = 1 + 2*ry;
    A_y(1, 2) = -2*ry;
    for j = 2:Ny-1
        A_y(j, j-1) = -ry;
        A_y(j, j)   = 1 + 2*ry;
        A_y(j, j+1) = -ry;
    end
    C_top_imp     = 2 * dy * h_top_half / k;
    A_y(Ny, Ny-1) = -2*ry;
    A_y(Ny, Ny)   = 1 + 2*ry + ry*C_top_imp;
    
    dA_y = decomposition(A_y);
    
    for i = 1:Nx
        B_y = zeros(Ny, 1);
        for j = 1:Ny
            if i == 1
                T_left_ghost = T_half(2,j) + C_left*(T_inf - T_half(1,j));
                B_y(j)       = rx*T_left_ghost + (1 - 2*rx)*T_half(i,j) + rx*T_half(i+1,j);
            elseif i == Nx
                T_right_ghost = T_half(Nx-1,j) + C_right*(T_inf - T_half(Nx,j));
                B_y(j)        = rx*T_half(i-1,j) + (1 - 2*rx)*T_half(i,j) + rx*T_right_ghost;
            else
                B_y(j) = rx*T_half(i-1,j) + (1 - 2*rx)*T_half(i,j) + rx*T_half(i+1,j);
            end
        end
        B_y(1)  = B_y(1)  + 2*ry * (dy * q_flux / k);
        B_y(Ny) = B_y(Ny) + ry * C_top_imp * T_inf;
        
        T_new(i, :) = (dA_y \ B_y)';
    end
end
