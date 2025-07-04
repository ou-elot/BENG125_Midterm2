a = 10;
b = 5;
x_fp = a / 5;
y_fp = 1 + x_fp^2;

[x, y] = meshgrid(0:2:50, 0:2:50);

dx = a - x - (4.*x.*y) ./ (1 + x.^2);
dy = b.*x .* (1 - y ./ (1 + x.^2));

figure;
quiver(x, y, dx, dy, 'k', 'AutoScale', 'on', 'AutoScaleFactor', 1.5); hold on;
contour(x, y, dx, [0 0], 'r', 'LineWidth', 2); % dx/dt = 0
contour(x, y, dy, [0 0], 'b', 'LineWidth', 2); % dy/dt = 0
plot(x_fp, y_fp, 'ko', 'MarkerSize', 8, 'LineWidth', 2); % fixed point

xlabel('x'); ylabel('y');
title('Stable Spiral: a = 10, b = 5');
legend('Vector field', 'dx/dt = 0', 'dy/dt = 0', 'Fixed point');
axis([0 50 0 50]); grid on;
