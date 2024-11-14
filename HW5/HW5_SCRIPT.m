% 公共参数
N = 109;
Lmax = 11;
h = Lmax / (N + 1);
i = 1:N;
x = (i - (N + 1) / 2) .* h;
a = (N+1) / 2 - 1 / h;
b = (N+1) / 2 + 1 / h;

% 电子参数
V0_e = 48;                % 电子势阱深度
me = 0.067;              % 电子有效质量 (相对于自由电子质量)

% 构建电子哈密顿量
H_e = eye(N, N) .* (2 + 2 * h^2 * V0_e);
for j = 1:N-1
    H_e(j, j+1) = -1;
    H_e(j+1, j) = -1;
    if j > a && j < b
        H_e(j, j) = 2;
    end
end
H_e = H_e / (2 * me * h^2);
[V_e, D_e] = eig(H_e);
eigenvalues_e = sort(diag(D_e));
GroundEnergy_e = eigenvalues_e(1);

% 空穴参数
V0_h = 48;                % 空穴势阱深度
mh = 0.45;               % 空穴有效质量 (相对于自由电子质量)

% 构建空穴哈密顿量
H_h = eye(N, N) .* (2 + 2 * h^2 * V0_h);
for j = 1:N-1
    H_h(j, j+1) = -1;
    H_h(j+1, j) = -1;
    if j > a && j < b
        H_h(j, j) = 2;
    end
end
H_h = H_h / (2 * mh * h^2);
[V_h, D_h] = eig(H_h);
eigenvalues_h = sort(diag(D_h));
GroundEnergy_h = eigenvalues_h(1);

% 物理常数
h_planck = 6.62607015e-34;    % J·s
c_light = 3e8;                % m/s
eV_to_J = 1.602176634e-19;    % J/eV

% 假设归一化单位对应于eV
% 若不是，需要根据实际情况转换
% 这里简单假设归一化能量已经是eV
Delta_E = GroundEnergy_e - GroundEnergy_h;  % eV
lambda = (h_planck * c_light) / (Delta_E * eV_to_J);  % m
lambda_nm = lambda * 1e9;                      % nm

% 输出结果
fprintf('电子基态能量: %.4f eV\n', GroundEnergy_e);
fprintf('空穴基态能量: %.4f eV\n', GroundEnergy_h);
fprintf('能量差: %.4f eV\n', Delta_E);
fprintf('发射波长: %.2f nm\n', lambda_nm);

% 绘制基态波函数
figure;
subplot(2,1,1);
plot(x, V_e(:,1));
xlabel('位置');
ylabel('电子基态波函数');
title('电子基态波函数');

subplot(2,1,2);
plot(x, V_h(:,1));
xlabel('位置');
ylabel('空穴基态波函数');
title('空穴基态波函数');

% 绘制能带
num_levels = 20;
figure;
hold on;
plot(1:num_levels, eigenvalues_e(1:num_levels), 'b-o', 'DisplayName', '电子能级');
plot(1:num_levels, eigenvalues_h(1:num_levels), 'r-s', 'DisplayName', '空穴能级');
xlabel('能级编号');
ylabel('能量 (eV)');
title('电子与空穴的能带结构');
legend;
hold off;
