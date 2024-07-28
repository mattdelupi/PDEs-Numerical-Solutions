close all; clc; clear;

% Parametri in ingresso
a = 1;
k = .1;
L = 1;

% Mesh spaziale
N = 150;
x = linspace(0, L, N);
xi = x(2:end-1); % punti interni
Dx = x(2) - x(1);

% Mesh temporale
Nt = 300;
Tf = 2 * L / a;
t = linspace(0, Tf, Nt);
Dt = t(2) - t(1);

% BCs
DphiB0 = cos(2*pi*a/L * t); DphiB0 = [DphiB0 DphiB0(1)]; % Neumann in x = 0
phiBL = pi * ones(1, Nt); % Dirichlet in x = L

% IC
phi0 = 2*pi * ones(1, N); phi0 = phi0(:);

%% Primo quesito
% Matrici dello schema
I = eye(N-1);
D1 = gallery('tridiag', N-1, -1, 0, 1);
D2 = gallery('tridiag', N-1, 1, -2, 1);
C = a * Dt / Dx;
beta = k * Dt / Dx^2;
A = I - beta*D2;
B = I - .5*C*D1 + .5*C^2*D2;

% Implementazione delle BCs
B(1, :) = zeros(1, N-1);
xs = x(1:3); xc = x(1); w = PesiDer(xs, xc, 1);
A(1, 1:3) = w;

% Matrice di Transizione
T = A \ B;

% Termine noto
q = zeros(N-1, 1);
q(end) = q(end) - .5*C*pi + beta*pi + .5*C^2*pi;

% Diagramma con evoluzione temporale
figure('Name', 'Quesito 1: Discretizzazione globale', 'NumberTitle', 'off');
an0 = animatedline('Color', 'r', 'LineStyle', '-', 'LineWidth', 1);
an = animatedline('Color', 'k', 'LineStyle', '-', 'LineWidth', 1);
addpoints(an0, x, phi0);
xlabel('x'); ylabel('\phi');
legend('Iniziale', 'Attuale');
axis([0 L 0 7.5]);
grid on;
phi = phi0(1:end-1);
for it = 1:Nt
    q(1) = DphiB0(it+1);
    phi = T*phi + A\q;
    pause(it/(it-.99)*Dt);
    title(['t = ' num2str(round(t(it), 2))]);
    clearpoints(an);
    addpoints(an, x, [phi; phiBL(it)]);
end

%% Secondo quesito
% Matrice D1
D1 = 1/(2*Dx) * gallery('tridiag', N-2, -1, 0, 1);
xs = x(1:3); xc = x(1); w = PesiDer(xs, xc, 1);
D1(1, 1:2) = D1(1, 1:2) + w(2:end)./(2*Dx*w(1));

% Matrice D2
D2 = 1/Dx^2 * gallery('tridiag', N-2, 1, -2, 1);
D2(1, 1:2) = D2(1, 1:2) - w(2:end)./(w(1)*Dx^2);

% Termine noto
q = zeros(N-2, 1);
q(end) = q(end) + (-a/(2*Dx) + k/Dx^2)*pi;
c1 = a / (2*Dx*w(1));
c2 = k / (Dx^2*w(1));
c = c1 + c2; % serve nella function

% Soluzione con ode45
ICs = phi0(2:end-1);
options = odeset('AbsTol', 1e-12, 'RelTol', 1e-12);
[t, y] = ode45(@(t, y) myRHS(t, y, L, a, k, q, c, D2, D1), t, ICs, options);

% Diagramma con evoluzione temporale
figure('Name', 'Quesito 2: Semidiscretizzazione', 'NumberTitle', 'off');
an0 = animatedline('Color', 'r', 'LineStyle', '-', 'LineWidth', 1);
an = animatedline('Color', 'k', 'LineStyle', '-', 'LineWidth', 1);
addpoints(an0, x, phi0);
xlabel('x'); ylabel('\phi');
legend('Iniziale', 'Attuale');
axis([0 L 0 1.2*max(y, [], 'all')]);
grid on;
for it = 1:Nt
    pause(it/(it-.99)*Dt);
    title(['t = ' num2str(round(t(it), 2))]);
    clearpoints(an);
    phiB0 = DphiB0(it)/w(1) - 1/w(1)*w(2:end)*y(it, 1:2).';
    addpoints(an, x, [phiB0 y(it, :) phiBL(it)]);
end

%% Function per il sistema di ODEs
function dydt = myRHS(t, y, L, a, k, q, c, D2, D1)
    y = y(:);
    q(1) = c * cos(2*pi*a/L * t); % Implementazione della BC alla Neumann
    
    % Matrice del sistema di ODEs
    M = -a*D1 + k*D2;

    dydt = M*y + q;
end