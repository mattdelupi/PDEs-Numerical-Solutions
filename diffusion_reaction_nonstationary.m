close all; clc; clear;

% Parametri in ingresso
L = 1;
b = 1;
k = @(x) 1 + sin(pi/L*x);

% Mesh spaziale
N = 50;
csi = linspace(0, 1, N);
x = L/2 * (1 - cos(pi*csi));
xi = x(2:end);

% Mesh temporale
tf = 2;
Nt = 120;
t = linspace(0, tf, Nt);
dt = t(2) - t(1);

% IC
phi0 = zeros(N, 1);

% BCs
% Dirichlet omogenea in x=0 e Neumann omogenea in x=L

% Matrice I
I = eye(N-1);

% Matrice Dk
Dk = zeros(N-1);
for i = 2:N-2
    xs = xi(i-1:i+1);
    xc = xi(i); wi = PesiDer(xs, xc, 1);
    xc = xi(i-1); wim1 = PesiDer(xs, xc, 1);
    xc = xi(i+1); wip1 = PesiDer(xs, xc, 1);
    wiki = wi .* k(xs);
    Dk(i, i-1:i+1) = wiki(1)*wim1 + wiki(2)*wi + wiki(3)*wip1;
end
xs = x(1:3);
xc = x(2); wi = PesiDer(xs, xc, 1);
xc = x(1); wim1 = PesiDer(xs, xc, 1);
xc = x(3); wip1 = PesiDer(xs, xc, 1);
wiki = wi .* k(xs);
Dk(1, 1:2) = wiki(1)*wim1(2:end) + wiki(2)*wi(2:end) + wiki(3)*wip1(2:end);

% Matrice A
A = (1+.5*b*dt)*I - .5*dt*Dk;
xs = x(end-2:end); xc = x(end); wN = PesiDer(xs, xc, 1);
A(end, end-2:end) = wN;

% Matrice B
B = (1-.5*b*dt)*I + .5*dt*Dk;
B(end, end) = 0;

% Termine noto
q = dt * ones(N-1, 1);
q(end) = 0;
q = A \ q;

% Matrice di transizione
T = A \ B;

% Soluzione
PHI = zeros(Nt, N-1);
PHI(1, :) = phi0(2:end);
for it = 2:Nt
    PHI(it, :) = T*PHI(it-1, :).' + q;
end
PHI = [zeros(Nt, 1) PHI];

% Diagramma
figure;
an0 = animatedline('Color', 'r', 'LineStyle', '--', 'LineWidth', 1);
an = animatedline('Color', 'k', 'LineStyle', '-', 'LineWidth', 1);
xlabel('x'); ylabel('\phi'); grid on;
axis([0 L min(PHI(:)) max(PHI(:))]);
legend('Iniziale', 'Attuale', 'Location', 'northwest');
addpoints(an0, x, phi0);
for it = 1:Nt
    pause(dt);
    title(['t = ' num2str(round(t(it), 2))]);
    clearpoints(an);
    addpoints(an, x, PHI(it, :));
    drawnow
end