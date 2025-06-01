clear; clc;

% Create grid in (a, b) parameter space
a_vals = linspace(0, 0.14, 40);
b_vals = linspace(0, 1.2, 40);

% Prepare figure
figure; hold on;
xlabel('a'); ylabel('b');
title('Q4 Hopf Bifurcation Curve');
axis([0 0.14 0 1.2]);
axis square;
grid on;

% Sweep through all combinations of (a, b)
for i = 1:length(a_vals)
    for j = 1:length(b_vals)
        a = a_vals(i);
        b = b_vals(j);

        % System of equations
        f = @(v) [-v(1) + a*v(2) + v(1)^2 * v(2);
                   b - a*v(2) - v(1)^2 * v(2)];

        % Use fsolve to find fixed point
        [sol, ~, exitflag] = fsolve(f, [0.5 0.5], optimset('Display','off'));

        if exitflag > 0  % if solution converged
            x = sol(1);
            y = sol(2);

            % Jacobian matrix at (x, y)
            J = [ -1 + 2*x*y, a + x^2;
                 -2*x*y,     -a - x^2];

            eigenvalues = eig(J);

            if any(imag(eigenvalues) ~= 0)  % if complex eigenvalues (spiral)
                if all(real(eigenvalues) < 0)
                    h1 = plot(a, b, 'kx', 'MarkerSize', 6);  % Stable spiral
                elseif any(real(eigenvalues) > 0)
                    h2 = plot(a, b, 'ko', 'MarkerSize', 6);  % Unstable spiral
                end
            end
        end
    end
end

%%
a_vals = linspace(0, 0.14, 100);
b_vals = linspace(0, 1.2, 100);

% Create empty matrices for spiral classification
spiral_type = zeros(length(b_vals), length(a_vals)); % 0: none, 1: stable, 2: unstable

% Loop through (a, b)
for i = 1:length(a_vals)
    for j = 1:length(b_vals)
        a = a_vals(i);
        b = b_vals(j);

        % System at (x, y)
        f = @(v) [-v(1) + a*v(2) + v(1)^2 * v(2);
                   b - a*v(2) - v(1)^2 * v(2)];
        [sol, ~, exitflag] = fsolve(f, [0.5 0.5], optimset('Display','off'));

        if exitflag > 0
            x = sol(1); y = sol(2);
            J = [ -1 + 2*x*y, a + x^2;
                 -2*x*y,     -a - x^2];
            eigs = eig(J);

            if any(imag(eigs) ~= 0)  % Spiral check
                if all(real(eigs) < 0)
                    spiral_type(j,i) = 1;  % stable spiral
                elseif any(real(eigs) > 0)
                    spiral_type(j,i) = 2;  % unstable spiral
                end
            end
        end
    end
end


% === Estimate and plot full Hopf boundary (stable <-> unstable)
hopf_lower = [];
hopf_upper = [];

for i = 1:length(a_vals)
    for j = 2:length(b_vals)
        t1 = spiral_type(j-1,i);
        t2 = spiral_type(j,i);

        b_mid = (b_vals(j-1) + b_vals(j)) / 2;
        a_now = a_vals(i);

        if a_now > 0.005
            if t1 == 1 && t2 == 2
                hopf_lower(end+1,:) = [a_now, b_mid]; %#ok<AGROW>
            elseif t1 == 2 && t2 == 1
                hopf_upper(end+1,:) = [a_now, b_mid]; %#ok<AGROW>
            end
        end
    end
end

% Sort each half by increasing a
hopf_lower = sortrows(hopf_lower, 1);
hopf_upper = sortrows(hopf_upper, 1);
hopf_closed = [hopf_lower; flipud(hopf_upper)];

N = 3;  % change this to 4–8 depending on how smooth you want it
hopf = hopf_closed(1:N:end, :);
h3 = plot(hopf(:,1), hopf(:,2), 'r-', 'LineWidth', 2);

legend([h1 h2 h3], {'Stable Spiral (x)', 'Unstable Spiral (o)', 'Hopf Bifurcation Curve'});

%%
clear; clc;

% Parameters
a = 0.06;
b = 0.2;

% System
dxdt = @(x, y) -x + a*y + x.^2 .* y;
dydt = @(x, y) b - a*y - x.^2 .* y;

% Grid
[x, y] = meshgrid(linspace(0, 3, 40), linspace(0, 3, 40));
u = dxdt(x, y);
v = dydt(x, y);

% Normalize vectors for direction field
L = sqrt(u.^2 + v.^2);
u = u ./ L;
v = v ./ L;

% Plot vector field
figure; hold on;
quiver(x, y, u, v, 0.5, 'k');

% Plot nullclines
fimplicit(@(x,y) dxdt(x,y), [0 3 0 3], 'r', 'LineWidth', 1.5);  % dx/dt = 0
fimplicit(@(x,y) dydt(x,y), [0 3 0 3], 'b', 'LineWidth', 1.5);  % dy/dt = 0

% Format
xlabel('x'); ylabel('y');
title('Phase Portrait at a = 0.06, b = 0.2 (Unstable Spiral)');
legend('Vector Field', 'dx/dt = 0', 'dy/dt = 0');
axis tight; axis square;
