a = 0.1; b = 0.2; omega =  -1;
fi_r_diff = omega; fi_r = @(t) pi/6+omega*t;

t_init = 0; % Initial time position.
t_end = 10; % Time at the end.
step = 0.02; % Step time for computing and plotting.
tol = 1e-8; % Tolerance for computing.

init_est = [0.5; 1]; % Initial estimate.

t_plot = t_init:step:t_end; % time array.
n_s = length(t_plot); % Number of solution.


x = zeros(2,n_s); % Matrix for storing solution.

ii = 1;

while t_init <= t_end
    F = @(x) [a*cos(fi_r(t_init))+b*cos(x(1))-x(2);...
              a*sin(fi_r(t_init))-b*sin(x(1))];

    J = @(x) [-b*sin(x(1)), -1;...
              -b*cos(x(1)),  0];

          x_ss = NR_method(F, J, init_est, tol);
          init_est = x_ss;
          x(:,ii)=x_ss;

          t_init = t_init + step;
          ii = ii + 1; 
end

hold on
plot(t_plot,x(1,:),'r-', 'LineWidth', 1.4)
plot(t_plot,x(2,:),'b-', 'LineWidth', 1.4)

fi_r_vec = pi/6+omega*t_plot;

theta_diff = (a*cos(fi_r_vec)*fi_r_diff)./(0.2*cos(x(1,:)));
d_diff = -a*sin(fi_r_vec)*fi_r_diff-b*sin(x(1,:)).*theta_diff;

plot(t_plot,theta_diff,'c-', 'LineWidth', 1.4)
plot(t_plot,d_diff,'g-', 'LineWidth', 1.4)
title('Position and Velocity')
xlabel('t')
ylabel('\theta [rad], d [m],\theta velocity [rad/s], d velocity [m/s]')
legend('\theta position','d position','\theta velocity', 'd velocity')
hold off
