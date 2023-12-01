function [t,e1,e2,e3,ne] = RK4System(ef, e0, tspan, h)
% RK4System resuelve un sistema de ecuaciones diferenciales de primer orden
% utilizando el método de Runge-Kutta de cuarto orden.
%   - ef: arreglo de celdas correspondiente al sistema de ecuaciones diferenciales.
% Las funciones en el vector deben tomar como entrada t, e1, e2 y e3 (en ese orden).
%   - e0: vector con las condiciones iniciales [e01 e02 e03].
%   - tspan: vector [t0 tf] que define el intervalo temporal.
%   - h: tamaño del paso.
% Output:
%   - t: vector con los valores de t en los que se calculó la solución.
%   - e1: vector con los valores de la solución de e1.
%   - e2: vector con los valores de la solución de e2.
%   - e3: vector con los valores de la solución de e3.
%   - ne: vector correspondiente a la norma de los vectores.

% Extraer los valores iniciales
t0 = tspan(1);
tf = tspan(2);

% Implementar el algoritmo de Runge-Kutta de cuarto orden
t = t0:h:tf;
e1 = zeros(size(t));
e2 = zeros(size(t));
e3 = zeros(size(t));
ne = zeros(size(t));
e1(1) = e0(1);
e2(1) = e0(2);
e3(1) = e0(3);
ne(1) = sqrt(e1(1)^2 + e2(1)^2 + e3(1)^2);

for i = 1:length(t)-1
    k1e1 = h*ef{1}(t(i),e1(i),e2(i),e3(i));
    k1e2 = h*ef{2}(t(i),e1(i),e2(i),e3(i));
    k1e3 = h*ef{3}(t(i),e1(i),e2(i),e3(i));
    k2e1 = h*ef{1}(t(i)+h/2,e1(i)+k1e1/2,e2(i)+k1e2/2,e3(i)+k1e3/2);
    k2e2 = h*ef{2}(t(i)+h/2,e1(i)+k1e1/2,e2(i)+k1e2/2,e3(i)+k1e3/2);
    k2e3 = h*ef{3}(t(i)+h/2,e1(i)+k1e1/2,e2(i)+k1e2/2,e3(i)+k1e3/2);
    k3e1 = h*ef{1}(t(i)+h/2,e1(i)+k2e1/2,e2(i)+k2e2/2,e3(i)+k2e3/2);
    k3e2 = h*ef{2}(t(i)+h/2,e1(i)+k2e1/2,e2(i)+k2e2/2,e3(i)+k2e3/2);
    k3e3 = h*ef{3}(t(i)+h/2,e1(i)+k2e1/2,e2(i)+k2e2/2,e3(i)+k2e3/2);
    k4e1 = h*ef{1}(t(i)+h,e1(i)+k3e1,e2(i)+k3e2,e3(i)+k3e3);
    k4e2 = h*ef{2}(t(i)+h,e1(i)+k3e1,e2(i)+k3e2,e3(i)+k3e3);
    k4e3 = h*ef{3}(t(i)+h,e1(i)+k3e1,e2(i)+k3e2,e3(i)+k3e3);
    e1(i+1) = e1(i) + (k1e1 + 2*k2e1 + 2*k3e1 + k4e1)/6;
    e2(i+1) = e2(i) + (k1e2 + 2*k2e2 + 2*k3e2 + k4e2)/6;
    e3(i+1) = e3(i) + (k1e3 + 2*k2e3 + 2*k3e3 + k4e3)/6;
    ne(i+1) = sqrt(e1(i)^2 + e2(i)^2 + e3(i)^2);
end
end