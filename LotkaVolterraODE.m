
function LotkaVolterraODE
clear x;

% time span in weeks
tspan = [0 100];

% starting values
G = 20;
R = 20;

%F = @(t,x) [(-1 * (2 * log(2)) / N) * x(2) * x(1); ((2 * log(2)) / N) * x(2) * x(1) - (log(2)) * x(2); (log(2)) * x(2)];

% integrate with the ode solver
[t, x] = ode45(@(t, x) odef(t, x), tspan, [G R]);

% plot the results
plot(t, x,'-o')
xlabel('t'); 

function funcs = odef(t, x)

% constants
k_0 = 1;
k_r = 1;
k_rd = 1;

% Define the function on the right side of the ode
% NOTE: x(1) - G x(2) - R
funcs = zeros(2, 1);
funcs(1) = k_0 * x(1) - k_r * x(1) * x(2);
funcs(2) = k_r * x(1) * x(2) - k_rd * x(2);

