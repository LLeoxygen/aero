%% 固定翼飞行器飞行包线 Flight Envelope
% 横坐标: TAS (真空速, m/s)
% 纵坐标: H (高度, m)
% 绘制内容包括：
%   - 失速速度边界 (V_S)
%   - 最大速度边界 (V_max)
%   - 实用升限 (Service Ceiling)
%   - 推力可用与推力需用平衡线
% 作者: 星野拾光 @ bilibili
% Author: LLeoxygen @ github.com

clear; clc; close all;

%% ==================== Aerocraft Parameters ====================

% 几何参数
S = 2.97; % 机翼参考面积, m^2
b = 6.4; % 翼展, m
AR = 13.8; % 展弦比
lambda = 2; % 后掠角, deg

% 质量与重量
m = 100; % 飞行器质量, kg
g = 9.81; % 重力加速度, m/s^2
W = m * g; % 重量, N

% 升力特性
CL_max = 1.4936; % 最大升力系数 (失速)
CL_min = 0.4851; % 最小升力系数
CLCD_max = 17.0527; % 最大升阻比

% 阻力特性
e = 0.85; % 奥斯瓦尔德效率因子
k = 1 / (pi * AR * e); % 诱导阻力因子
CD0 = CL_max / CLCD_max - CL_max ^ 2 * k; % 零升阻力系数

% 发动机推力特性 (螺旋桨或涡轮喷气式)
P0 = 9370; % 海平面最大可用功率 (W) 或 海平面最大推力 (N)
% 注：对于螺旋桨飞机用功率，对于喷气式飞机用推力
engine_type = 'prop'; % 发动机类型: 'prop'(螺旋桨) 或 'jet'(喷气)
eta_prop = 0.85; % 螺旋桨效率 (仅螺旋桨飞机使用)

% 速度限制
V_NE = 160; % 结构极限速度 - 永不超越速度, m/s (约140节)
use_V_NE = false; % 是否启用结构极限速度边界 (true/false)

% 高度范围
H_max = 5000; % 最大飞行高度, m

%% ==================== Flight Envelope ====================

% 高度网格
H_points = linspace(0, H_max, 100);
nH = length(H_points);

% 初始化边界数组
V_S = zeros(1, nH); % 失速速度
V_max_power = zeros(1, nH); % 推力限制的最大速度
V_min_power = zeros(1, nH); % 推力限制的最小速度
V_Y = zeros(1, nH); % 最大爬升率速度
ROC_max = zeros(1, nH); % 最大爬升率
Ceiling_H = []; % 实用升限

% 各高度的速度边界
for i = 1:nH
    H = H_points(i);
    [rho, ~, ~, ~] = atmos_ISA(H);

    % 1. 最小平飞速度
    V_S(i) = sqrt((2 * W) / (rho * S * CL_max));

    % 2. 最大爬升率速度 V_Y
    V_Y(i) = sqrt(2 / rho * W / S * sqrt(k / (3 * CD0)));

    % 3. 推力限制的速度边界
    if strcmp(engine_type, 'prop')
        % 螺旋桨飞机：功率分析
        [~, P_avail] = available_thrust_power(H, P0, engine_type, eta_prop);

        % 求解 P_req(V) = P_avail
        if use_V_NE
            V_test_max = V_NE * 1.2;
        else
            V_test_max = 500; % 禁用结构限速时，使用足够大的搜索上限
        end

        V_test = linspace(V_S(i) * 1.01, V_test_max, 500);
        P_diff = zeros(size(V_test));

        for j = 1:length(V_test)
            [~, P_req] = required_thrust_power(V_test(j), H, W, S, CD0, k, eta_prop);
            P_diff(j) = P_avail - P_req;
        end

        % 找到功率平衡的交点
        idx_pos = find(P_diff >= 0, 1, 'first');
        idx_neg = find(P_diff >= 0, 1, 'last');

        if ~isempty(idx_pos) && ~isempty(idx_neg) && idx_pos ~= idx_neg

            if idx_pos > 1
                V_min_power(i) = V_test(idx_pos - 1) - P_diff(idx_pos - 1) * (V_test(idx_pos) - V_test(idx_pos - 1)) / (P_diff(idx_pos) - P_diff(idx_pos - 1));
            else
                V_min_power(i) = V_test(idx_pos);
            end

            if idx_neg < length(V_test)
                V_max_power(i) = V_test(idx_neg) - P_diff(idx_neg) * (V_test(idx_neg + 1) - V_test(idx_neg)) / (P_diff(idx_neg + 1) - P_diff(idx_neg));
            else
                V_max_power(i) = V_test(idx_neg);
            end

        else
            V_min_power(i) = NaN;
            V_max_power(i) = NaN;
        end

        % 计算该高度下的最大爬升率 (ROC = (P_avail - P_req(V_Y))/W)
        [~, P_req_Y] = required_thrust_power(V_Y(i), H, W, S, CD0, k, eta_prop);
        ROC_max(i) = (P_avail - P_req_Y) / W;

    else
        % 喷气式飞机：推力分析
        [T_avail, ~] = available_thrust_power(H, P0, engine_type, eta_prop);

        % 求解 T_req(V) = T_avail
        if use_V_NE
            V_test_max = V_NE * 1.2;
        else
            V_test_max = 500; % 禁用结构限速时，使用足够大的搜索上限
        end

        V_test = linspace(V_S(i) * 1.01, V_test_max, 500);
        T_diff = zeros(size(V_test));

        for j = 1:length(V_test)
            [T_req, ~] = required_thrust_power(V_test(j), H, W, S, CD0, k, eta_prop);
            T_diff(j) = T_avail - T_req;
        end

        idx_pos = find(T_diff >= 0, 1, 'first');
        idx_neg = find(T_diff >= 0, 1, 'last');

        if ~isempty(idx_pos) && ~isempty(idx_neg) && idx_pos ~= idx_neg

            if idx_pos > 1
                V_min_power(i) = V_test(idx_pos - 1) - T_diff(idx_pos - 1) * (V_test(idx_pos) - V_test(idx_pos - 1)) / (T_diff(idx_pos) - T_diff(idx_pos - 1));
            else
                V_min_power(i) = V_test(idx_pos);
            end

            if idx_neg < length(V_test)
                V_max_power(i) = V_test(idx_neg) - T_diff(idx_neg) * (V_test(idx_neg + 1) - V_test(idx_neg)) / (T_diff(idx_neg + 1) - T_diff(idx_neg));
            else
                V_max_power(i) = V_test(idx_neg);
            end

        else
            V_min_power(i) = NaN;
            V_max_power(i) = NaN;
        end

        % 计算该高度下的最大爬升率 (ROC = (T_avail - T_req(V_Y))*V_Y/W)
        [T_req_Y, ~] = required_thrust_power(V_Y(i), H, W, S, CD0, k, eta_prop);
        ROC_max(i) = (T_avail - T_req_Y) * V_Y(i) / W;

    end

end

% 找到 ROC=0.5 m/s 时的实用升限
idx_ROC = find(ROC_max >= 0.5, 1, 'last');

if ~isempty(idx_ROC)
    Ceiling_H_ROC = H_points(idx_ROC);
end

%% ==================== Display ====================

figure('Position', [100 100 900 700], 'Color', 'white');
hold on;
grid on;

% 各边界线
plot(V_S, H_points, 'y-', 'LineWidth', 2, 'DisplayName', '失速边界 (V_{stall})');
plot(V_max_power, H_points, 'r-', 'LineWidth', 2, 'DisplayName', '推力限制边界');

% 结构限制速度
if use_V_NE
    V_NE_line = V_NE * ones(size(H_points));
    plot(V_NE_line, H_points, 'm-', 'LineWidth', 2, 'DisplayName', '结构极限速度 (V_{NE})');
end

% 填充飞行包线区域
if use_V_NE
    V_upper = min(V_max_power, V_NE);
else
    V_upper = V_max_power;
end

V_lower = max(V_S, V_min_power);
V_lower(isnan(V_lower)) = V_S(isnan(V_lower));

% 创建闭合包线进行填充
valid_range = ~isnan(V_upper) & ~isnan(V_lower);

if any(valid_range)
    V_fill = [V_lower(valid_range), fliplr(V_upper(valid_range)), V_lower(find(valid_range, 1))];
    H_fill = [H_points(valid_range), fliplr(H_points(valid_range)), H_points(find(valid_range, 1))];
    fill(V_fill, H_fill, [0.6 1 0.6], 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'DisplayName', '飞行包线区域');
end

% 实用升限@ROC=0.5 m/s
if ~isempty(Ceiling_H_ROC)
    idx_R = idx_ROC;
    plot([V_lower(idx_R), V_upper(idx_R)], [Ceiling_H_ROC, Ceiling_H_ROC], 'b-', 'LineWidth', 1.5, 'DisplayName', sprintf('实用升限 ROC=0.5m/s (%.0f m)', Ceiling_H_ROC));
end

xlabel('真空速 TAS (m/s)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('高度 H (m)', 'FontSize', 12, 'FontWeight', 'bold');

if strcmp(engine_type, 'prop')
    engine_str = '螺旋桨';
else
    engine_str = '喷气式';
end

title(sprintf('飞行包线@星野拾光@bilibili'), 'FontSize', 14, 'FontWeight', 'bold');

legend('Location', 'best', 'FontSize', 10);

if use_V_NE
    xlim([0, V_NE * 1.1]);
else
    xlim([0, max(V_max_power(~isnan(V_max_power))) * 1.1]);
end

ylim([0, H_max]);

% 参数标注
annotation('textbox', [0.15, 0.65, 0.25, 0.25], 'String', { ...
                                               '飞行器参数:', ...
                                               sprintf('质量: %.0f kg', m), ...
                                               sprintf('翼面积: %.1f m²', S), ...
                                               sprintf('展弦比: %.1f', AR), ...
                                               sprintf('C_{L,max}: %.2f', CL_max), ...
                                               sprintf('C_{D0}: %.3f', CD0), ...
                                               sprintf('功率: %.0f kW', P0 / 1000)}, ...
    'FitBoxToText', 'on', 'BackgroundColor', [1 1 0.9], 'EdgeColor', [0.5 0.5 0.5], ...
    'FontSize', 9);

hold off;

%% ==================== Write to console ====================
fprintf('\n========== 飞行包线关键参数 ==========\n');
fprintf('海平面失速速度: %.1f m/s (%.1f km/h)\n', V_S(1), V_S(1) * 3.6);
fprintf('海平面最大速度: %.1f m/s (%.1f km/h)\n', V_max_power(1), V_max_power(1) * 3.6);
fprintf('平飞实用升限: %.0f m\n', Ceiling_H);

if ~isempty(Ceiling_H_ROC)
    fprintf('实用升限 (ROC=0.5m/s): %.0f m\n', Ceiling_H_ROC);
end

if use_V_NE
    fprintf('结构极限速度: %.1f m/s (%.1f km/h)\n', V_NE, V_NE * 3.6);
end

fprintf('=====================================\n\n');

%% ==================== functions ====================

function [rho, T, p, c] = atmos_ISA(H)
    % 国际标准大气 (ISA) 模型
    % 输入: H - 高度 (m)
    % 输出: rho - 空气密度 (kg/m^3)
    %       T - 温度 (K)
    %       p - 压力 (Pa)
    %       c - 声速 (m/s)

    % 海平面参数
    rho0 = 1.225; % 海平面密度 (kg/m^3)
    T0 = 288.15; % 海平面温度 (K)
    p0 = 101325; % 海平面压力 (Pa)
    R = 287.05; % 气体常数 (J/kg·K)
    gamma = 1.4; % 比热比
    L = 0.0065; % 温度递减率 (K/m)
    g = 9.81; % 重力加速度

    if H <= 11000 % 对流层
        T = T0 - L * H;
        p = p0 * (T / T0) ^ (g / (R * L));
    else % 平流层（简化）
        T = 216.65; % 恒温
        p = p0 * (1 - L * 11000 / T0) ^ (g / (R * L)) * exp(-g * (H - 11000) / (R * T));
    end

    rho = p / (R * T);
    c = sqrt(gamma * R * T);
end

function [T_avail, P_avail] = available_thrust_power(H, P0, engine_type, eta_prop)
    % 计算可用推力或功率
    % 输入: H - 高度 (m), P0 - 海平面功率/推力
    %       engine_type - 发动机类型: 'prop' 或 'jet'
    %       eta_prop - 螺旋桨效率
    % 输出: T_avail - 可用推力 (N), P_avail - 可用功率 (W)

    [rho, ~, ~, ~] = atmos_ISA(H);
    rho0 = 1.225;

    sigma = rho / rho0; % 密度比

    if strcmp(engine_type, 'prop')
        % 螺旋桨飞机：功率随高度近似恒定（直至临界高度）
        P_avail = P0 * (0.8 + 0.2 * sigma); % 功率略有下降
        T_avail = []; % 螺旋桨推力与速度相关
    else
        % 喷气式飞机：推力随密度下降
        T_avail = P0 * sigma; % 近似推力与密度成正比
        P_avail = [];
    end

end

function [T_req, P_req] = required_thrust_power(V, H, W, S, CD0, k, eta_prop)
    % 计算给定速度下维持平飞所需的推力/功率
    % 输入: V - 真空速 (m/s)
    %       H - 高度 (m)
    %       W - 重量 (N)
    %       S - 机翼面积 (m^2)
    %       CD0 - 零升阻力系数
    %       k - 诱导阻力系数
    %       eta_prop - 螺旋桨效率
    % 输出: T_req - 需用推力 (N), P_req - 需用功率 (W)

    [rho, ~, ~, ~] = atmos_ISA(H);

    CL = (2 * W) / (rho * V ^ 2 * S); % 升力系数
    CD = CD0 + k * CL ^ 2; % 阻力系数

    T_req = 0.5 * rho * V ^ 2 * S * CD; % 需用推力
    P_req = T_req * V; % 需用功率
end
