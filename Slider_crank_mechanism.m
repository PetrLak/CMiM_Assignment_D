a = 0.1; b = 0.2; fi_d = 30; omega =  -1;
fi_r = pi * fi_d/180; fi_r_diff = omega;

F = @(x) [a*cos(fi_r)+b*cos(x(1))-x(2);...
          a*sin(fi_r)-b*sin(x(1))];

J = @(x) [-b*sin(x(1)), -1;...
          -b*cos(x(1)),  0];

init_est = [0.5; 0.3];
tol = 1e-12;

[x, n] = NR_method(F, J, init_est, tol);

theta_diff = (a*cos(fi_r)*fi_r_diff)/(0.2*cos(x(1)));
d_diff = -a*sin(fi_r)*fi_r_diff-b*sin(x(1))*theta_diff;

fprintf('theta = %.5g° (%.5g rad), d = %.5g m\n', 180*x(1)/pi, x(1), x(2))
fprintf('theta_velocity = %.5g rad/s, d_velocity = %.5g m/s\n', theta_diff, d_diff)
