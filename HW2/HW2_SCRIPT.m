global sigma beta rho
sigma = 10;
beta = 8/3;
rho = 28;
T = 100;
phi0 = [-7; 7; 25];

function dphidt = LorentzFunc(t, phi)
    global sigma beta rho
    dphidt = [sigma*(phi(2)-phi(1));        % dx/dt = sigma * (y-x)
              phi(1)*(rho-phi(3))-phi(2);   % dy/dt = x * (rho-z) - y
              phi(1)*phi(2)-beta*phi(3)];   % dz/dt = x * y - beta * z
end

[t, phi] = ode45(@LorentzFunc, [0, T], phi0);