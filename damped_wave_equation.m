close all; clc; clear;

%% Parametri in ingresso
L = 1;
a = 1;
alfa = .5;

Tg = L / a;
n = 1;
Tf = n * Tg;

%% Definizione mesh spaziale, funzione b(x) e mesh temporale
N = 100;
x = linspace(0, L, N);
xi = x(2:end-1); % punti interni
b = alfa * xi.^2;
Nt = 200;
t = linspace(0, Tf, Nt);
dt = Tf / Nt;

%% Primo quesito
% BCs
phiB = 0;
% ICs
phi0 = exp(-100 * (x./L - .5).^2); phi0 = phi0(:);
Dphi0 = zeros(N, 1);

% Matrice B
B = diag(b);

% Matrice D2
D2 = zeros(N-2);
for i = 3:N-4
    xs = xi(i-2:i+2);
    xc = xi(i);
    w = PesiDer(xs, xc, 2);
    D2(i, i-2:i+2) = w;
end
xs = x(1:6); xc = x(2); w = PesiDer(xs, xc, 2);
D2(1, 1:5) = w(2:end);
xs = x(1:5); xc = x(3); w = PesiDer(xs, xc, 2);
D2(2, 1:4) = w(2:end);
xs = x(end-4:end); xc = x(end-2); w = PesiDer(xs, xc, 2);
D2(end-1, end-3:end) = w(1:end-1);
xs = x(end-5:end); xc = x(end-1); w = PesiDer(xs, xc, 2);
D2(end, end-4:end) = w(1:end-1);

% Matrice del sistema di equazioni differenziali
I = eye(N-2);
O = zeros(N-2);
M = [O I; a^2*D2 -B];

% Risoluzione con ode45
ICs = [phi0(2:end-1); Dphi0(2:end-1)];
options = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);
[t, y] = ode45(@(t, y) myRHS(t, y, M), t, ICs, options);
Y = y(:, 1:N-2);

% Diagramma dell'evoluzione temporale usando animatedline
figure('Name', 'Quesito 1', 'NumberTitle', 'off');
an0 = animatedline('Color', 'r', 'LineStyle', '-', 'LineWidth', 1);
an = animatedline('Color', 'k', 'LineStyle', '-', 'LineWidth', 1);
legend('Iniziale', 'Attuale');
grid on;
xlabel('x'); ylabel('\phi');
axis([0 L min(Y, [], 'all') max(phi0)]);
addpoints(an0, x, phi0);
for it = 1:Nt
    pause(dt);
    clearpoints(an);
    addpoints(an, x, [phiB Y(it, :) phiB]);
    title(['t = ' num2str(round(it*dt, 2))]);
    drawnow
end

%% Secondo quesito
% BCs
phiB0 = 0; % Dirichlet
DphiBL = 0; % Neumann
% ICs
phi0 = exp(-100 * (x./L - .5).^2); phi0 = phi0(:);
Dphi0 = zeros(N, 1);

% Pesi derivata prima per BC alla Neumann
xs = x(end-4:end); xc = x(end);
wN = PesiDer(xs, xc, 1);

% Matrice D2
xs = x(end-4:end); xc = x(end-2); w = PesiDer(xs, xc, 2);
D2(end-1, end-3:end) = w(1:end-1);
D2(end-1, end-3:end) = D2(end-1, end-3:end) - w(end)/wN(end)*wN(1:end-1);
xs = x(end-5:end); xc = x(end-1); w = PesiDer(xs, xc, 2);
D2(end, end-4:end) = w(1:end-1);
D2(end, end-3:end) = D2(end, end-3:end) - w(end)/wN(end)*wN(1:end-1);

% Matrice del sistema di equazioni differenziali
M = [O I; a^2*D2 -B];

% Risoluzione con ode45
[t, y] = ode45(@(t, y) myRHS(t, y, M), t, ICs, options);
Y = y(:, 1:N-2);

% Diagramma dell'evoluzione temporale usando animatedline
figure('Name', 'Quesito 2', 'NumberTitle', 'off');
an0 = animatedline('Color', 'r', 'LineStyle', '-', 'LineWidth', 1);
an = animatedline('Color', 'k', 'LineStyle', '-', 'LineWidth', 1);
legend('Iniziale', 'Attuale');
grid on;
xlabel('x'); ylabel('\phi');
axis([0 L min(Y, [], 'all') max(phi0)]);
addpoints(an0, x, phi0);
for it = 1:Nt
    pause(dt);
    clearpoints(an);
    phiBL = -1/wN(end)*wN(1:end-1)*Y(it, end-3:end).';
    addpoints(an, x, [phiB0 Y(it, :) phiBL]);
    title(['t = ' num2str(round(it*dt, 2))]);
    drawnow
end

%% Funzione che definisce il sistema di ODEs
function dydt = myRHS(~, y, M)
    y = y(:);
    dydt = M * y;
end