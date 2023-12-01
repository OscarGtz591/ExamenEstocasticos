%% PARTE 2: SIR Model
clc; clear; close all
B = 1.8; % Beta
G = 0.1; % Gamma
sigma = 0.19;
N = 5000; % Población
SIR0 = [4899 100 1 0]; % Condiciones iniciales
tspan = [0 100];
h = 0.001;
f1 = @(t,S,E,I,R) -B*S.*(I/N);
f2 = @(t,S,E,I,R) (B*S.*I)/N - sigma*E;
f3 = @(t,S,E,I,R) sigma*E - G*I;
f4 = @(t,S,E,I,R) G*I;
SIR = {f1, f2, f3,f4};
[t,S,E,I,R] = RK4System(SIR, SIR0, tspan, h);
figure(1)
hold on
plot(t, S, 'b', 'LineWidth', 1);
plot(t, E, 'r', 'LineWidth', 1); 
plot(t, I, 'g', 'LineWidth', 1);
plot(t, R, 'y', 'LineWidth', 1);
hold off

title('Modelo SIR: $\beta = 0.5, \gamma = 0.5$','Interpreter','latex')
xlabel('Tiempo','Interpreter','latex')
ylabel('Poblaci\''on','Interpreter','latex')
legend('S: Suceptibles', 'E: Expuestos', 'I: Infectados', 'R: Recuperados', 'Interpreter', 'latex', 'Location', 'bestoutside')
grid on