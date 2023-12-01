%% SEIR ESTOCÁSTICO: Algoritmo de Gillespie

clc; clear; close all;

% Parámetros del modelo SEIR
beta = 1.8;    % Tasa de transmisión
gamma = 0.1;   % Tasa de recuperación
sigma = 0.5;     % Inverso del período de incubación

% Condiciones iniciales
S1 = 4899;     % Susceptibles iniciales
E1 = 100;      % Expuestos iniciales
I1 = 1;        % Infectados iniciales
R1 = 0;        % Recuperados iniciales

% Número de iteraciones
N = 10;

% Tiempo total de simulación
Ts = 100;

% Inicializar vectores para almacenar resultados
tiempo = zeros(N, Ts);
PS = zeros(N, Ts);
PE = zeros(N, Ts);
PI = zeros(N, Ts);
PR = zeros(N, Ts);

% Inicializar el tiempo y las poblaciones
tiempo(:, 1) = 0;
PS(:, 1) = S1;
PE(:, 1) = E1;
PI(:, 1) = I1;
PR(:, 1) = R1;

% Ciclo de simulación utilizando el algoritmo de Gillespie
for i = 1:N
    for t = 2:Ts
        % Calcular tasas de eventos
        tasa_infeccion = beta * PS(i, t-1) * PI(i, t-1) / (PS(i, t-1) + PE(i, t-1) + PI(i, t-1) + PR(i, t-1));
        tasa_incubacion = sigma * PE(i, t-1);
        tasa_recuperacion = gamma * PI(i, t-1);
        tasa_total = tasa_infeccion + tasa_incubacion + tasa_recuperacion;
    
        % Calcular tiempo hasta el próximo evento
        tiempo(i, t) = tiempo(i, t-1) - log(rand)/tasa_total;
    
        % Determinar tipo de evento (infección, incubación o recuperación)
        r = rand;
        if r < tasa_infeccion/tasa_total
            % Evento de infección
            PS(i, t) = PS(i, t-1) - 1;
            PE(i, t) = PE(i, t-1) + 1;
        elseif r < (tasa_infeccion + tasa_incubacion)/tasa_total
            % Evento de incubación
            PE(i, t) = PE(i, t-1) - 1;
            PI(i, t) = PI(i, t-1) + 1;
        else
            % Evento de recuperación
            PI(i, t) = PI(i, t-1) - 1;
            PR(i, t) = PR(i, t-1) + 1;
        end
    end
end

% Graficar resultados de 10 corridas
figure;
hold on;
for i = 1:10
    plot(tiempo(i, :), PS(i, :), 'b-', tiempo(i, :), PE(i, :), 'g-', tiempo(i, :), PI(i, :), 'r-', tiempo(i, :), PR(i, :), 'm-');
end
hold off;
title('Modelo SEIR - Simulación Estocástica - 10 Corridas');
xlabel('Tiempo');
ylabel('Población');
legend('Susceptibles', 'Expuestos', 'Infectados', 'Recuperados');
grid on;
