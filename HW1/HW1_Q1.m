% 参数设置
I_c = 1e-11; % 理想电流
R = 10; % 内阻
T = 300; % 温度
k_0 = physconst("Boltzmann"); % 玻尔兹曼常数
q = 1.602176634e-19; % 电子带电量
Vmin = -1;
Vmax = 1;
V = Vmin:0.01:Vmax;
I = zeros(size(V)); % 创建电流数组
delta = 0.001; % 定义相对精度
max_iter = 1e3; % 定义最大迭代数
i_init = 0; % 定义电流初始值

i_fzero = zeros(size(V)); % 创建新电流数组
index = 1;
for v = Vmin:0.01:Vmax
   
    fun = @(x) I_c * (exp(q*(v - x*R)/(k_0*T)) - 1.0) - x;
    x0 = 0.0;
    root = fzero(fun, x0); % 调用fzero求解方程
    i_fzero(index) = root;
    index = index + 1;
end