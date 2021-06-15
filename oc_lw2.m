%% Практикум по ОУ

% Лабораторная работа №2
% Нелинейная задача ОУ - задача распределения капиталовложений
% Вариант 13

clear;
clc;

%% входные данные

dots_per_conj_grid = 15; % количество точек в сетке сопряжённых переменных

t_0 = 0; % начальный момент времени
T = 7; % конечный момент времени

% коэффициенты в производственной функции F
global lambda;
lambda = 0.8; % in (0, 1)
if lambda <= 0 || lambda >= 1
    error('Input parameters: lambda is in (0, 1)');
end

K_0 = 30; % начальный объём производственных фондов
global L; 
L = 50; % объём доступных трудовых ресурсов

% коэффиценты системы
global gamma; 
gamma = 0.1; % > 0, коэффициент естественной убыли загрязнений
if gamma <= 0
    error('Input parameters: gamma > 0');
end

global delta;
delta = 2; % коэффициент уменьшения загрязнений за счёт затрат продукции

global mu; 
mu = 0.2; % > 0, коэффициент амортизации основного капитала
if mu <= 0
    error('Input parameters: mu > 0');
end

global epsilon;
epsilon = 4; % > 0, доля объёма загрязнений относительно объёма производства
if epsilon <= 0
    error('Input parameters: epsilon > 0');
end

P_0 = 350; % начальный объём загрязнений, отходов производства

% коэффициент в совокупной пользе
global r; 
r = 10; % > 0, коэффициент дисконтирования
if r <= 0
    error('Input parameters: r > 0');
end

% коэффициенты в функции полезности U
global A;
A = 3; % > 0
if A <= 0
    error('Input parameters: A > 0');
end

global B; 
B = 1; % > 0
if B <= 0
    error('Input parameters: B > 0');
end

C = 1; % > 0
if C <= 0
    error('Input parameters: C > 0');
end

global alpha;
alpha = 5;

global beta;
beta = 0.4;
%% параметризация

theta = linspace(- pi / 2, 0, dots_per_conj_grid);
phi = linspace(-pi, 0, dots_per_conj_grid);

% psi_1, psi_2
psi_1_0_grid = cos(theta) .* cos(phi);
psi_2_0_grid = cos(theta) .* sin(phi);

% psi_0
psi_0_grid = sin(theta);

% производственная функция
global F;
F = @(K) (K .^ lambda) * (L ^ (1 - lambda));

% производная F по K
global F_K; 
F_K = @(K) lambda * (K .^ (lambda - 1)) * (L ^ (1 - lambda));

% функция полезности
U = @(c, P) C - A * exp(-alpha * c) + B * exp(-beta * P);

% параметры для ode45
options = odeset('RelTol', 1e-6, 'AbsTol', 1e-6);

sol_start = zeros(6, 1);
sol_start(1, 1) = K_0;   % K
sol_start(2, 1) = P_0;   % P

% управление
u = zeros(2, 1);

M_2_3_4_temp = @(t) exp(-r * t);

%%

for psi_idx = 1 : size(psi_1_0_grid, 2)
    sol_start(3, 1) = psi_1_0_grid(1, psi_idx); % psi_1
    sol_start(4, 1) = psi_2_0_grid(1, psi_idx); % psi_2
    
    psi_0 = psi_0_grid(1, psi_idx);
    
    % вычисление u_0
    [u] = calc_control(t_0, psi_0, sol_start(1, 1), sol_start(3, 1), sol_start(4, 1));
    sol_start(5, 1) = u(1); % доля продукции, выделяемая на потребление
    sol_start(6, 1) = u(2); % доля продукции, выделяемая на борьбу с загрязнениями
    
    [t, sol] = ode45(...
            @(t, y) system(t, y, psi_0), ...
            [t_0, T], sol_start, options);
   
    c = sol(:, 5)' .* F(sol(:, 1))'; % объём потребления продукции
    U_val = U(c, sol(:, 2)');  % функция полезности
    
    t_temp = t';
    total_benefit_fun = U_val .* M_2_3_4_temp(t_temp);
    total_benefit_fun_temp = cumsum(...
        (total_benefit_fun(1 : size(total_benefit_fun, 2) - 1) + ...
        total_benefit_fun(2 : size(total_benefit_fun, 2))) .* ...
        (t_temp(2 : size(t_temp, 2)) - t_temp(1 : size(t_temp, 2) - 1)), 2);
    
    total_benefit_int = 0.5 * total_benefit_fun_temp(1, size(total_benefit_fun_temp, 2));
    
    if psi_idx > 1
        if total_benefit_int > total_benefit_int_opt
            total_benefit_int_opt = total_benefit_int;
            sol_opt = sol';
            t_opt = t_temp;
        end
    else
        total_benefit_int_opt = total_benefit_int;
        sol_opt = sol';
        t_opt = t_temp;
    end
end


%% Результаты

% оптимальное управление
figure('Position', [100 100 1000 1000]);

subplot(2, 1, 1);
plot(t_opt, sol_opt(5, :)); % u_1(t)
title('u_1(t) - доля продукции, выделяемая на потребление');
xlabel('t');
ylabel('u_1');

subplot(2, 1, 2);
plot(t_opt, sol_opt(6, :)); % u_2(t)
title('u_2(t) - доля продукции, выделяемая на борьбу с загрязнениеми');
xlabel('t');
ylabel('u_2');

% сопряжённые переменные
figure('Position', [100 100 1000 1000]);

subplot(2, 1, 1);
plot(t_opt, sol_opt(3, :)); % psi_1(t)
title('\psi_1(t)');
xlabel('t');
ylabel('\psi_1');

subplot(2, 1, 2);
plot(t_opt, sol_opt(4, :)); % psi_2(t)
title('\psi_2(t)');
xlabel('t');
ylabel('\psi_2');

% оптимальная траектория
figure('Position', [100 100 1000 1000]);
plot(sol_opt(1, :), sol_opt(2, :)); % P(K)
title('P(K)');
xlabel('K - объём производственных фондов');
ylabel('P - объём загрязнений');


function [u] = calc_control(time, psi_0, K, psi_1, psi_2)
    global alpha;
    global A;
    global delta;
    global r;
    global F;
    
    F_val = F(K);
    
    % вспомогательные переменные для вычисления M_i
    
    M_array = zeros(1, 4);

    M_2_4_temp = @(psi_0_val) alpha * A * psi_0_val;
    M_2_temp = @(psi_1_val, temp) -psi_1_val / temp; % temp = M_2_4_temp_val
    M_4_temp = @(psi_2_val, temp) (delta * psi_2_val) / temp; % temp = M_2_4_temp_val

    M_2_3_4_temp = @(t) exp(-r * t);
    M_2_4_5_temp = @(F_val, temp) exp(-alpha * F_val) * temp; % F = F(K, L), temp = M_2_3_4_temp_val

    M_3_temp1 = @(psi_0_val, temp) A * psi_0_val * temp; % temp = M_2_3_4_temp_val
    M_3_temp2 = @(psi_1_val, psi_2_val) psi_1_val + delta * psi_2_val;

    M_2_4_temp_val = M_2_4_temp(psi_0);
    
    % u^* in AB
    
    M_2_temp_val = M_2_temp(psi_1, M_2_4_temp_val);
    M_2_3_4_temp_val = M_2_3_4_temp(time);
    M_2_4_5_temp_val = M_2_4_5_temp(F_val, M_2_3_4_temp_val);

    
    if (psi_1 > 0) && ...
        (M_2_4_5_temp_val < M_2_temp_val) && ...
        (M_2_temp_val < M_2_3_4_temp_val)
        M_2 = (-psi_1 / alpha) * (1.0 - r * time - log(M_2_temp_val));
    else
        M_2 = NaN;
    end
    M_array(1, 1) = M_2;
    
    % u^* : u_1^* = 0
    
    M_3_temp1_val = M_3_temp1(psi_0, M_2_3_4_temp_val);
    M_3_temp2_val = M_3_temp2(psi_1, psi_2);
    
    if M_3_temp2_val >= 0
        M_3 = M_3_temp1_val;
    else
        M_3 = M_3_temp1_val - F_val * M_3_temp2_val;
    end
    M_array(1, 2) = M_3;
        
    % u^* : u_1^* + u_2^* = 1
    
    M_4_temp_val = M_4_temp(psi_2, M_2_4_temp_val);
    
    if (M_2_4_5_temp_val < M_4_temp_val) && ... 
        (M_4_temp_val < M_2_3_4_temp_val)
        M_4 = -psi_1 * F_val + ...
            (delta * psi_2 / alpha) * (1 - alpha * F_val - r * time - log(M_4_temp_val));
    else
        M_4 = NaN;
    end
    M_array(1, 3) = M_4;
    
    % u^* = (1, 0)
    
    M_5 = A * psi_0 * M_2_4_5_temp_val - F_val * psi_1;
    M_array(1, 4) = M_5;
        
    % выбор наибольшего из значений M_2, M_3, M_4, M_5
    [M, M_idx] = max(M_array, [], 'omitnan');
    
    switch M_idx
        case 1
            u(1, 1) = (-1 / (alpha * F_val)) * (r * time + log(M_2_temp_val));
            u(2, 1) = 0;
        case 2
            u(1, 1) = 0;
            
            if M_3_temp2_val > 0
                u(2, 1) = 0;
            else % M_3_temp2_val < 0
                u(2, 1) = 1;
            end
            
        case 3
            u(1, 1) = (-1 / (alpha * F_val)) * (r * time + log(M_4_temp_val));
            u(2, 1) = 1 - u(1, 1);
            
        otherwise
            u(1, 1) = 1;
            u(2, 1) = 0;
    end 
end


function[sol] = system(t, y, psi_0)
    % y = [y_1, y_2, y_3, y_4, y_5, y_6]':
    % y_1 = K
    % y_2 = P
    % y_3 = psi_1
    % y_4 = psi_2
    % y_5 = u_1
    % y_6 = u_2
    
    global r;
    global A;
    global F;
    global F_K;
    global epsilon;
    global delta;
    global mu;
    global gamma;
    global alpha;
    global beta;
    global B;
    
    M_2_3_4_temp = @(t_val) exp(-r * t_val);
    M_3_temp1 = @(temp) A * psi_0 * temp; % temp = M_2_3_4_temp_val

    u = calc_control(t, psi_0, y(1, 1), y(3, 1), y(4, 1));
    sol(5, 1) = u(1, 1);
    sol(6, 1) = u(2, 1);
    
    F_val = F(y(1, 1));
    F_K_val = F_K(y(1, 1));
    
    M_2_3_4_temp_val = M_2_3_4_temp(t);
    M_3_temp1_val = M_3_temp1(M_2_3_4_temp_val);
    
    temp1 = 1 - u(1, 1) - u(2, 1);
    temp2 = epsilon - delta * u(2, 1);

    sol(1, 1) = temp1 * F_val - mu * y(1, 1);
    
    sol(2, 1) = temp2 * F_val - gamma * y(2, 1);
    
    sol(3, 1) = mu * y(3, 1) - ...
            F_K_val * (y(3, 1) * temp1 + y(4, 1) * temp2 - ...
            alpha * u(1, 1) * M_3_temp1_val * exp(-alpha * u(1, 1) * F_val));
        
    sol(4, 1) = gamma * y(4, 1) - beta * B * psi_0 * exp(-beta * y(2, 1)) * M_2_3_4_temp_val;
end