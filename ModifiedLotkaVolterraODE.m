function ModifiedLotkaVolterraODE
clear x;

% time span in weeks
tspan = [0 100];

% carrying capacity of prey
G_0 = 8; 

% constants
k_0  = 0.6; % intrinsic rate of prey population increase
k_r  = 1.325; % predation rate coefficient
k_rd = 1; % mortality rate of predator
K_g  = 1;

% initial values
G   = 5; % prey population
R   = 6; % predator population

% integrate with the ode solver
[t, x] = ode45(@(t, x) f(t, x, G_0, k_0, k_r, k_rd, K_g), tspan, [G, R]);

% plot the results
plot(t, x,'-o')
xlabel('t'); 
legend('A','B')

function funcs = f(t, x, G_0, k_0, k_r, k_rd, K_g)

% Define the function on the right side of the ode
% x(1) - G x(2) - R
funcs    = zeros(2, 1);
funcs(1) = k_0 * (G_0 - x(1)) * x(1) - k_r * x(1) * x(2)/(K_g + x(1));
funcs(2) = k_r * x(1) * x(2)/(K_g + x(1)) - k_rd * x(2);