r = 1;

[x, y] = meshgrid(linspace(0, 4.5, 25), linspace(0, 4, 25));
z = 1 - x - y;

dx = r .* x .* z;
dy = r .* y .* z;

figure; hold on;
streams = streamslice(x, y, dx, dy);
set(streams, 'Color', [0.4 0.7 1]);  % light blue arrows

h1 = plot(zeros(1, 100), linspace(0, 4, 100), 'r--', 'LineWidth', 2);       % x = 0 → dx/dt = 0
h2 = plot(linspace(0, 4.5, 100), zeros(1, 100), 'b--', 'LineWidth', 2);     % y = 0 → dy/dt = 0
h3 = fimplicit(@(x, y) x + y - 1, [0 4.5 0 4], 'k--', 'LineWidth', 2);      % x + y = 1

xlabel('x');
ylabel('y');
title('Streamlines with Nullclines');

legend([h1 h2 h3], {'\it dx/dt = 0', '\it dy/dt = 0', '\it x + y = 1'}, 'Location', 'northeast');

xlim([-0.5 4.5]);
ylim([-0.5 4]);
axis equal;
box on;
