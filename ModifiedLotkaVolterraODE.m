function ModifiedLotkaVolterraODE
clear x;

% time span in weeks
tspan = [0 100];

% carrying capacity of prey
X_0 = 8;

% constants
k_x  = 0.3; % intrinsic rate of prey population increase
k_xy  = 1.325; % predation rate coefficient
k_y = 1; % mortality rate of predator
K_g  = 1;

% initial values
X   = 5; % prey population
Y   = 6; % predator population

% integrate with the ode solver
[t, x] = ode45(@(t, x) f(t, x, X_0, k_x, k_xy, k_y, K_g), tspan, [X, Y]);

% plot the results
plot(t, x,'-o')
xlabel('t'); 
legend('X','Y')

function funcs = f(t, x, X_0, k_x, k_xy, k_y, K_g)

% Define the function on the right side of the ode
% x(1) - X x(2) - Y
funcs    = zeros(2, 1);
funcs(1) = k_x * (X_0 - x(1)) * x(1) - k_xy * x(1) * x(2)/(K_g + x(1));
funcs(2) = k_xy * x(1) * x(2)/(K_g + x(1)) - k_y * x(2);
