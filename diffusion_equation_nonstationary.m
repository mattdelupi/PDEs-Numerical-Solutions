close all; clc; clear;

%% Primo quesito
eqSolver(1, 1, 25, 1.7, .25, .65, 1);

%% Secondo quesito
N = 100;
alpha = linspace(0, .5, N);
beta = linspace(0, .5, N);
[A, B] = meshgrid(alpha, beta);
norma = zeros(N);
for i = 1:N
    for j = 1:N
        norma(i, j) = eqSolver(1, 1, 25, 1.7, A(i, j), B(i, j), 0);
    end
end
figure('Units', 'normalized', 'Position', [.1 .2 .4 .7]);
surf(A, B, norma);
title('Analisi di stabilità');
xlabel('\alpha'); ylabel('\beta'); zlabel('norma');
grid on;
shading interp;
hold on;
contour3(A, B, norma, linspace(min(norma(:)), max(norma(:)), 15), 'k-');

beta = linspace(0, 1, N);
[A, B] = meshgrid(alpha, beta);
norma = zeros(N);
for i = 1:N
    for j = 1:N
        norma(i, j) = eqSolver(1, 1, 25, 1.7, A(i, j), B(i, j), 0);
    end
end
figure('Units', 'normalized', 'Position', [.55 .2 .4 .7]);
contour(A, B, norma, [1 1], 'r-');
title('Curva per norma = 1');
xlabel('\alpha'); ylabel('\beta');
grid on;
text(.1, .2, "Stabile all'interno");

%% Function che risolve e calcola la norma
function norma = eqSolver(k, L, N, final_t, alpha, beta, ifplot)
    % Mesh spaziale
    x = linspace(0, L, N);
    dx = x(2) - x(1);

    % Mesh temporale
    dt = beta * dx^2 / k;
    t = 0:dt:final_t;
    Nt = length(t);

    % BCs
    % Neumann omogenea in x = 0
    phiBL = 1;

    % IC
    phi0 = 3*cos(3*pi/(2*L) * x) + (x./L).^2;

    % Matrice D2s
    D2s = gallery('tridiag', N-2, 1, -2, 1); D2s = full(D2s);
    xs = x(2:4); xc = x(2); w = dx^2 * PesiDer(xs, xc, 2);
    D2s(1, 1:3) = w;
    xs = x(end-3:end-1); xc = x(end-1); w = dx^2 * PesiDer(xs, xc, 2);
    D2s(end, end-2:end) = w;

    % Matrice I
    I = eye(N-2);
    I(1, 1) = I(1, 1) + beta - 2*alpha*beta^2;
    I(2, 2) = I(2, 2) + alpha*beta^2;

    % Matrice D2
    D2 = gallery('tridiag', N-2, 1, -2, 1); D2 = full(D2);
    
    % Matrice di transizione
    T = I + beta*D2 + alpha*beta^2*D2*D2;
    norma = norm(T);

    % Termine noto
    q = zeros(N-2, 1);
    q(end) = q(end) + beta*phiBL - 2*alpha*beta^2*phiBL;
    q(end-1) = q(end-1) + alpha*beta^2*phiBL;

    if ifplot == 1
        % Soluzione
        PHI = zeros(Nt, N-2);
        PHI(1, :) = phi0(2:end-1);
        for it = 2:Nt
            PHI(it, :) = T*PHI(it-1, :).' + q;
        end
        PHI = [PHI(:, 1) PHI ones(Nt, 1)];
        PHI(1, :) = phi0;
    
        % Diagramma
        figure;
        an0 = animatedline('Color', 'r', 'LineStyle', '--', 'LineWidth', 1);
        an = animatedline('Color', 'k', 'LineStyle', '-', 'LineWidth', 1);
        xlabel('x'); ylabel('\phi');
        legend('Iniziale', 'Attuale');
        grid on;
        axis([0 L min(PHI(:)) max(PHI(:))]);
        addpoints(an0, x, phi0);
        for it = 1:Nt
            pause(dt);
            clearpoints(an);
            addpoints(an, x, PHI(it, :));
            title(['t = ' num2str(round(t(it), 2))]);
        end
    end
end