close all; clc; clear;

% Mesh spaziale
L = 1;
N = 80;
h = L / N;
x = linspace(0, L-h, N);

% Mesh temporale
Tf = L;
Nt = 160;
Dt = Tf / Nt;
t = linspace(0, Tf, Nt);

% BCs periodiche

% ICs
u0 = exp(-100/L^2 * (x-L/2).^2);
v0 = zeros(1, N);
ICs = [u0 v0];

%% Primo quesito
% Matrice D1
circ(1) = 0; circ(2) = 1; circ(N) = -1;
D1 = gallery('circul', circ) ./ (2*h);

% Matrice M del sistema di ODEs ottenuto
M = [-D1 D1; D1 zeros(N)];

% Runge-Kutta e diagramma della soluzione
figure('Name', 'Quesito 1: u', 'NumberTitle', 'off', ...
    'Units', 'normalized', 'Position', [.1 .3 .4 .5]);
axu = gca; set(gca, 'TickLabelInterpreter', 'latex');
uan0 = animatedline('Color', 'r', 'LineStyle', '--', 'LineWidth', 1);
uan = animatedline('Color', 'k', 'LineStyle', '-', 'LineWidth', 1);
axis([0 L -1 1]);
grid on;
xlabel('x'); ylabel('u');
legend('Iniziale', 'Attuale');
figure('Name', 'Quesito 1: v', 'NumberTitle', 'off', ...
    'Units', 'normalized', 'Position', [.55 .3 .4 .5]);
axv = gca;
van0 = animatedline('Color', 'r', 'LineStyle', '--', 'LineWidth', 1);
van = animatedline('Color', 'k', 'LineStyle', '-', 'LineWidth', 1);
addpoints(uan0, [x L], [u0 u0(1)]);
addpoints(van0, [x L], [v0 v0(1)]);
axis([0 L -1 1]);
grid on;
xlabel('x'); ylabel('u');
legend('Iniziale', 'Attuale');
y = ICs(:);
for it = 1:Nt
    tn = it * Dt;
    f1 = myRHS(tn, y, M);
    f2 = myRHS(tn + Dt/3, y + Dt/3*f1, M);
    f3 = myRHS(tn + 2*Dt/3, y + Dt*(-1/3*f1+f2), M);
    f4 = myRHS(tn + Dt, y + Dt*(f1-f2+f3), M);
    y = y + Dt*(1/8*f1+3/8*f2+3/8*f3+1/8*f4);
    u = y(1:N);
    v = y(N+1:end);
    pause(Dt);
    title(axu, ['t = ' num2str(round(it*Dt, 2))]);
    title(axv, ['t = ' num2str(round(it*Dt, 2))]);
    clearpoints(uan); clearpoints(van);
    addpoints(uan, [x L], [u; u(1)]); addpoints(van, [x L], [v; v(1)]);
    drawnow
end

%% Secondo quesito
% Risoluzione con ode45
options = odeset('AbsTol', 1e-12, 'RelTol', 1e-12);
[t, y] = ode45(@(t, y) myRHS(t, y, M), t, ICs, options);
u = y(:, 1:N);
v = y(:, N+1:end);

% Diagramma della soluzione
figure('Name', 'Quesito 2: u', 'NumberTitle', 'off', ...
    'Units', 'normalized', 'Position', [.1 .3 .4 .5]);
axu = gca;
uan0 = animatedline('Color', 'r', 'LineStyle', '--', 'LineWidth', 1);
uan = animatedline('Color', 'k', 'LineStyle', '-', 'LineWidth', 1);
axis([0 L -1 1]);
grid on;
xlabel('x'); ylabel('u');
legend('Iniziale', 'Attuale');
figure('Name', 'Quesito 2: v', 'NumberTitle', 'off', ...
    'Units', 'normalized', 'Position', [.55 .3 .4 .5]);
axv = gca;
van0 = animatedline('Color', 'r', 'LineStyle', '--', 'LineWidth', 1);
van = animatedline('Color', 'k', 'LineStyle', '-', 'LineWidth', 1);
addpoints(uan0, [x L], [u0 u0(1)]);
addpoints(van0, [x L], [v0 v0(1)]);
axis([0 L -1 1]);
grid on;
xlabel('x'); ylabel('v');
legend('Iniziale', 'Attuale');
for it = 1:Nt
    pause(Dt);
    title(axu, ['t = ' num2str(round(it*Dt, 2))]);
    title(axv, ['t = ' num2str(round(it*Dt, 2))]);
    clearpoints(uan); clearpoints(van);
    addpoints(uan, [x L], [u(it, :) u(it, 1)]);
    addpoints(van, [x L], [v(it, :) v(it, 1)]);
    drawnow
end

%% Terzo quesito
% Matrici
I = eye(N);
D1 = gallery('circul', circ);
circ(1) = -2; circ(2) = 1; circ(N) = 1;
D2 = gallery('circul', circ);

% Matrice di transizione
T = [I-Dt/(2*h)*D1+.5*D2 Dt/(2*h)*D1; Dt/(2*h)*D1 I+.5*D2];
nT = norm(T); disp(['nT = ' num2str(nT)]);

% Diagramma della soluzione
figure('Name', 'Quesito 3: u', 'NumberTitle', 'off', ...
    'Units', 'normalized', 'Position', [.1 .3 .4 .5]);
axu = gca;
uan0 = animatedline('Color', 'r', 'LineStyle', '--', 'LineWidth', 1);
uan = animatedline('Color', 'k', 'LineStyle', '-', 'LineWidth', 1);
axis([0 L -1 1]);
grid on;
xlabel('x'); ylabel('u');
legend('Iniziale', 'Attuale');
figure('Name', 'Quesito 3: v', 'NumberTitle', 'off', ...
    'Units', 'normalized', 'Position', [.55 .3 .4 .5]);
axv = gca;
van0 = animatedline('Color', 'r', 'LineStyle', '--', 'LineWidth', 1);
van = animatedline('Color', 'k', 'LineStyle', '-', 'LineWidth', 1);
addpoints(uan0, [x L], [u0 u0(1)]);
addpoints(van0, [x L], [v0 v0(1)]);
axis([0 L -1 1]);
grid on;
xlabel('x'); ylabel('v');
legend('Iniziale', 'Attuale');
y = ICs(:);
for it = 1:Nt
    y = T * y;
    u = y(1:N);
    v = y(N+1:end);
    pause(Dt);
    title(axu, ['t = ' num2str(round(it*Dt, 2))]);
    title(axv, ['t = ' num2str(round(it*Dt, 2))]);
    clearpoints(uan); clearpoints(van);
    addpoints(uan, [x L], [u; u(1)]); addpoints(van, [x L], [v; v(1)]);
    drawnow
end

%% Function del sistema di ODEs
function dydt = myRHS(~, y, M)
    y = y(:);
    dydt = M * y;
end