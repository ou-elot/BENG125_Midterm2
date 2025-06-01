%%% CONFIGURATION %%%

% Define the system using symbolic variables for full modularity
syms x y
fx_sym = x * (3 - 2*x - 2*y);   % dx/dt
fy_sym = y * (2 - y - x);     % dy/dt

% Axis limits
x_min = 0; x_max = 4;
y_min = 0; y_max = 4;

% Resolution for vector field and nullclines
grid_step = 0.1;
nullcline_resolution = 400;

%%%%%%%%%%%%%%%%%%%%%%

% Convert symbolic expressions to numeric function handles
fx = matlabFunction(fx_sym, 'Vars', [x y]);
fy = matlabFunction(fy_sym, 'Vars', [x y]);

% Solve for nullclines symbolically (dx/dt = 0 and dy/dt = 0)
nullcline_fx = solve(fx_sym == 0, y);  % dx/dt = 0
nullcline_gy = solve(fy_sym == 0, y);  % dy/dt = 0

% Convert symbolic nullclines to function handles
nc_fx_fun = arrayfun(@(expr) matlabFunction(expr, 'Vars', x), nullcline_fx, 'UniformOutput', false);
nc_gy_fun = arrayfun(@(expr) matlabFunction(expr, 'Vars', x), nullcline_gy, 'UniformOutput', false);

% Create a grid for the vector field
[xg, yg] = meshgrid(x_min:grid_step:x_max, y_min:grid_step:y_max);
u = fx(xg, yg);
v = fy(xg, yg);

% Plot the vector field
figure;
streamslice(xg, yg, u, v);  % Correct 2D usage
hold on;

% Prepare x-range for nullcline plotting
x_vals = linspace(x_min, x_max, nullcline_resolution);
y_vals = linspace(y_min, y_max, nullcline_resolution);

% Plot dx/dt = 0 nullclines (including vertical x = 0 if applicable)
r1 = plot(zeros(size(y_vals)), y_vals, 'r--', 'LineWidth', 2);  % x = 0
for k = 1:length(nc_fx_fun)
    y_nc = arrayfun(nc_fx_fun{k}, x_vals);
    plot(x_vals, y_nc, 'r--', 'LineWidth', 2);
end

% Plot dy/dt = 0 nullclines (including y = 0)
b1 = plot(x_vals, zeros(size(x_vals)), 'b--', 'LineWidth', 2);  % y = 0
for k = 1:length(nc_gy_fun)
    y_nc = arrayfun(nc_gy_fun{k}, x_vals);
    plot(x_vals, y_nc, 'b--', 'LineWidth', 2);
end

% Final plot formatting
axis([x_min x_max y_min y_max]);
axis equal;
xlabel('x');
ylabel('y');
title('Streamlines with Nullclines');
legend([r1, b1], {'dx/dt = 0', 'dy/dt = 0'}, 'Location', 'best');
grid on;
