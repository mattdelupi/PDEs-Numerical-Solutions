close all; clc; clear;

%% Parametri in ingresso
L = 1;
a = 1;
N = 70;
h = L / N;
x = 0:h:L-h;
Nt = 600;
Tg = L / a;
n = 1; % numero di giri
Tf = n * Tg;
dt = Tf / Nt;

%% Primo quesito
cou = a * dt / h;
if abs(cou) > 1
    warning(['Attenzione! Numero di Courant instabile. C = ' num2str(cou)]);
end

phi0 = exp(-10 * (x-.5).^2); phi0 = phi0(:);
T = Tmatrix(L, N, a, dt);
nT = norm(T);
if nT > 1
    warning(['Attenzione! Norma instabile. nT = ' num2str(nT)]);
end

phi = phi0;
figure('Name', 'Quesito 1', 'NumberTitle', 'off');
an0 = animatedline('Color', 'r', 'LineStyle', '-', 'LineWidth', 1);
an = animatedline('Color', 'k', 'LineStyle', '-', 'LineWidth', 1);
legend('Iniziale', 'Attuale');
grid on;
xlabel('x'); ylabel('\phi');
axis([0 L min(phi0) max(phi0)]);
addpoints(an0, [x L], [phi0; phi0(1)]);
for it = 1:Nt
    phi = T * phi;
    clearpoints(an);
    addpoints(an, [x L], [phi; phi(1)]);
    title(['t = ' num2str(round(it*dt, 2))]);
    drawnow
end

%% Secondo quesito
C = linspace(0, 3, 300);
Dt = h / a * C;
norma = zeros(300, 1);
for i = 1:300
    T = Tmatrix(L, N, a, Dt(i));
    norma(i) = norm(T);
end
figure('Name', 'Quesito 2', 'NumberTitle', 'off');
plot(C, norma, 'k-', 'LineWidth', 1); hold on;
plot(C, ones(300, 1), 'r-');
axis([0 3 0 max(norma)]);
grid on;
legend('Norma', '1');
xlabel('Numero di Courant'); ylabel('Norma di T');

%% Terzo quesito
Cou = cou; % fissato uguale a quello del primo quesito
Nx = 10:1:110;
h = L ./ Nx;
Dt = Cou / a * h;
err = zeros(length(Nx), 1);
for i = 1:length(Nx)
    x = 0:h(i):L-h(i);
    phi0 = exp(-10 * (x-.5).^2); phi0 = phi0(:);
    T = Tmatrix(L, Nx(i), a, Dt(i));
    [phif, phif_an] = soluz(phi0, T, Dt(i), L, a);
    err(i) = norm(phif - phif_an, "inf");
end
figure('Name', 'Quesito 3', 'NumberTitle', 'off');
plot(Nx, err, 'k.-');
grid on;
xlabel('Punti di discretizzazione'); ylabel('Norma inf Errore');
title(['Valutazione errore dopo un giro. C = ' num2str(Cou)]);

%% Funzioni utili
function T = Tmatrix(L, N, a, dt)
    h = L / N;

    I = eye(N);
    
    v(1) = 4;
    v(2) = 1;
    v(N) = 1;
    A = gallery('circul', v);
    v(1) = 0;
    v(2) = 1;
    v(N) = -1;
    B = 3/h * gallery('circul', v);
    D = A \ B;

    alfa = .1;
    beta = 6 / (5*h^2);
    v(1) = 1;
    v(2) = alfa;
    v(N) = alfa;
    A = gallery('circul', v);
    v(1) = -2;
    v(2) = 1;
    v(N) = 1;
    B = beta * gallery('circul', v);
    D2 = A \ B;

    T = I - a*dt*D + .5*a^2*dt^2*D2;
end

function [phif, phif_an] = soluz(phi0, T, dt, L, a)
    Tg = L / a;
    Nt = Tg / dt;
    
    phi = phi0;
    for it = 1:Nt
        phi = T * phi;
    end
    phif = [phi; phi(1)];
    phif_an = [phi0; phi0(1)];
end