% 设定程序参数
M = linspace(0, 2*pi, 50);
epsilon = 0.7;
p = 1.1;
delta = 1e-3; % 相对精度

% 根据M求E
E = zeros(size(M)); % 创建E数组
index = 1; % 初始化索引
func_M_E = @(m, e) e - epsilon * sin(e) - m;
for m = M
    func = @(e) func_M_E(e, m);
    % 不动点法开始
    e0 = 0; % E初始值
    % 不动点法迭代
    while true
        e = func(e0);
        if abs((e - e0) / e) < delta % 相对精度条件
            root = e0;
            break;
        elseif e == 0 % 考虑除数为零情况
            root = e0;
            break;
        end
        e0 = e;
    end
    E(index) = root; % 存储数值解
    index = index + 1;
end

Theta = 2 * atan(sqrt((1+epsilon)/(1-epsilon)) * tan(E / 2));
R = p ./ (1+epsilon*cos(Theta));

x = R .* cos(Theta);
y = R .* sin(Theta);
plot(x, y, "r.");