close all; clc; clear;

%% Primo quesito
eqSolver(50, 4e-3, 1, .5, .1);

%% Secondo quesito
eqSolver(50, 4e-3, 5, 0, 0);

%% Terzo quesito
eqSolver(30, 4e-3, 1, .5, .1);
eqSolver(50, 4e-3, 1, .5, .1);
eqSolver(100, 1e-3, 1, .5, .1);

%% Function per risolvere l'equazione
function eqSolver(N, dt, final_t, gamma, epsilon)
    % Parametri da fissare
    a = 1.5;
    k = .02;

    % Mesh spaziale
    psi = linspace(0, 1, N);
    x = psi.^(pi/4) ./ max(psi.^(pi/4));
    x = x(1:end-1);
    N = N-1;

    % Mesh temporale
    t = 0:dt:final_t;
    Nt = length(t);

    % IC
    phi0 = zeros(1, N);
    for i = 1:N
        if abs(x(i) - .15) <= .05
            phi0(i) = -2;
        elseif abs(x(i) - .4) <= .1
            phi0(i) = 2.2;
        elseif abs(x(i) - .8) <= .1
            phi0(i) = .6;
        end
    end
    
    % Matrice D1
    D1 = zeros(N);
    for i = 3:N
        xs = x(i-2:i); xc = x(i); w = PesiDer(xs, xc, 1);
        D1(i, i-2:i) = w;
    end
    xs = [x(end-1:end)-1 x(1)]; xc = x(1); w = PesiDer(xs, xc, 1);
    D1(1, [end-1:end 1]) = w;
    xs = [x(end)-1 x(1:2)]; xc = x(2); w = PesiDer(xs, xc, 1);
    D1(2, [end 1:2]) = w;

    % Matrice D2
    D2 = zeros(N);
    for i = 3:N-1
        xs = x(i-2:i+1); xc = x(i); w = PesiDer(xs, xc, 2);
        D2(i, i-2:i+1) = w;
    end
    xs = [x(end-1:end)-1 x(1:2)]; xc = x(1); w = PesiDer(xs, xc, 2);
    D2(1, [end-1:end 1:2]) = w;
    xs = [x(end)-1 x(1:3)]; xc = x(2); w = PesiDer(xs, xc, 2);
    D2(2, [end 1:3]) = w;
    xs = [x(end-2:end)-1 x(1)]; xc = x(end)-1; w = PesiDer(xs, xc, 2);
    D2(end, [end-2:end 1]) = w;

    % Matrice P
    p = gamma + epsilon*cos(2*pi*x);
    P = diag(p);

    % Array di Butcher
    b = [2/9 1/3 4/9];
    A = [0 0 0; 1/2 0 0; 0 3/4 0];

    % Soluzione
    PHI = zeros(Nt, N);
    PHI(1, :) = phi0;
    for it = 2:Nt
        PHI(it, :) = RKsolv(t(it), PHI(it-1, :), a, D1, k, D2, P, dt, b, A);
    end
    PHI = [PHI PHI(:, 1)];
    x = [x 1];

    % Diagramma con evoluzione temporale
    figure;
    an0 = animatedline('Color', 'r', 'LineStyle', '--', 'LineWidth', 1);
    an = animatedline('Color', 'k', 'LineStyle', '-', 'LineWidth', 1);
    legend('Iniziale', 'Attuale');
    xlabel('x'); ylabel('\phi');
    grid on;
    axis([0 1 min(PHI(:)) max(PHI(:))]);
    addpoints(an0, x, PHI(1, :));
    for it = 1:Nt
        pause(dt);
        clearpoints(an);
        addpoints(an, x, PHI(it, :));
        title(['t = ' num2str(round(t(it), 2))]);
    end
end

%% Function per il Runge Kutta
function y_new = RKsolv(t, y_old, a, D1, k, D2, P, dt, b, A)
    y_old = y_old(:);
    N = length(y_old);
    s = length(b);
    c = sum(A, 2);

    F = zeros(N, s);
    for j = 1:s
        F(:, j) = myRHS(t+dt*c(j), y_old+dt*F*A(j, :).', a, D1, k, D2, P);
    end

    y_new = y_old + dt*F*b(:);
end

%% Function che definisce il RHS del sistema di ODEs
function dydt = myRHS(~, y, a, D1, k, D2, P)
    y = y(:);
    dydt = (-a*D1 + k*D2)*y + P*(y.^2);
end