% Script to solve homework 3 of computational physics
% Author: Haoxuan Zhu

% Params init
v_0 = 0;
theta_0 = pi / 4;
T = 2000;
h = 0.2;

% Array init
t = 0:0.2:2000;
v_midi = zeros(size(t));
theta_i = zeros(size(t));
v_midi(1) = -h / 2 * sin(theta_0);
theta_i(1) = theta_0;

% Recursive solve
[~, max_i] = size(t);
for i = 2:max_i
    theta_i(i) = theta_i(i-1) + h * v_midi(i-1);
    v_midi(i) = v_midi(i-1) + h * (-sin(theta_i(i)));
end

% Reconstruct v_i array
v_next_midi = circshift(v_midi, -1);
v_next_midi = v_next_midi(1:end-1);
v_form_midi = v_midi(1:end-1);
v_i = (v_form_midi + v_next_midi) ./ 2;
v_i = [v_0, v_i];

% Calculate energy
% Consider the lowest position of the pendulum as the zero potential line
E = v_i .^ 2 - 2 * cos(theta_i);
figure;
plot(t, E);
xlabel("Time");
ylabel("Energy");