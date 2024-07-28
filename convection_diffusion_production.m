close all; clc; clear;

%% Primo quesito
x = linspace(0, 1, 25);
eqSolver(x, 1, .1, 3, 1, "beta", .3);
pause;

%% Secondo quesito
x = linspace(0, 1, 65);
eqSolver(x, 10, .1, .5, 1, "cou", .5);
pause;
th = linspace(0, pi, 35);
x = 1/2 * (1 - cos(th));
eqSolver(x, 10, .1, .5, 0, "cou", 5);

%% Function che risolve l'equazione
function eqSolver(x, a, k, final_t, ifRK, which, param)
    % Mesh spaziale
    L = x(end);
    N = length(x);
    dx = min(diff(x));
    xi = x(2:end-1);

    % Mesh temporale
    if which == "beta"
        dt = param * dx^2 / k;
    elseif which == "cou"
        dt = param * dx / a;
    end
    Nt = ceil(final_t / dt);
    t = linspace(0, final_t, Nt);
    dt = t(2) - t(1);

    % BCs
    phiB0 = 1;
    phiBL = pi;

    % IC
    phi0 = 3 * exp(-250 * (x-L/2).^4);

    % Matrice D1
    D1 = zeros(N-2);
    for i = 3:N-3
        xs = xi(i-2:i+1); xc = xi(i); w = PesiDer(xs, xc, 1);
        D1(i, i-2:i+1) = w;
    end
    W = zeros(2, 4);
    xs = x(1:4); xc = x(3); w = PesiDer(xs, xc, 1);
    D1(2, 1:3) = w(2:end); W(1, 2) = w(1);
    xs = x(1:4); xc = x(2); w = PesiDer(xs, xc, 1);
    D1(1, 1:3) = w(2:end); W(1, 1) = w(1);
    xs = x(end-3:end); xc = x(end-1); w = PesiDer(xs, xc, 1);
    D1(end, end-2:end) = w(1:end-1); W(1, 4) = w(end);

    % Matrice D2
    D2 = zeros(N-2);
    for i = 3:N-4
        xs = xi(i-2:i+2); xc = xi(i); w = PesiDer(xs, xc, 2);
        D2(i, i-2:i+2) = w;
    end
    xs = x(1:5); xc = x(3); w = PesiDer(xs, xc, 2);
    D2(2, 1:4) = w(2:end); W(2, 2) = w(1);
    xs = x(1:6); xc = x(2); w = PesiDer(xs, xc, 2);
    D2(1, 1:5) = w(2:end); W(2, 1) = w(1);
    xs = x(end-4:end); xc = x(end-2); w = PesiDer(xs, xc, 2);
    D2(end-1, end-3:end) = w(1:end-1); W(2, 3) = w(end);
    xs = x(end-5:end); xc = x(end-1); w = PesiDer(xs, xc, 2);
    D2(end, end-4:end) = w(1:end-1); W(2, 4) = w(end);

    % Matrice del sistema di ODEs
    M = -a*D1 + k*D2;

    % Termine noto
    q = .5 + cos(2*pi/L * xi);
    q = q(:);
    q([1:2 end-1:end]) = q([1:2 end-1:end]) + diag([phiB0 phiB0 phiBL phiBL]) * W.' * [-a; k];

    % Soluzione
    if ifRK == 1
        b = [2/9 1/3 4/9];
        A = [0 0 0; 1/2 0 0; 0 3/4 0];
        PHI = zeros(Nt, N-2);
        PHI(1, :) = phi0(2:end-1);
        for it = 2:Nt
            PHI(it, :) = RKsolve(t(it), PHI(it-1, :), M, q, dt, b, A);
        end
    else
        ICs = phi0(2:end-1);
        options = odeset('AbsTol', 1e-12, 'RelTol', 1e-12);
        [t, PHI] = ode45(@(t, y) myRHS(t, y, M, q), t, ICs, options);
    end
    PHI = [phiB0*ones(Nt, 1) PHI phiBL*ones(Nt, 1)];
    PHI(1, :) = phi0;

    % Diagramma con evoluzione temporale
    figure;
    an0 = animatedline('Color', 'r', 'LineStyle', '--', 'LineWidth', 1);
    an = animatedline('Color', 'k', 'LineStyle', '-', 'LineWidth', 1);
    legend('Iniziale', 'Attuale');
    axis([0 L min(PHI(:)) max(PHI(:))]);
    grid on;
    xlabel('x'); ylabel('\phi');
    addpoints(an0, x, phi0);
    for it = 1:Nt
        pause(dt);
        clearpoints(an);
        addpoints(an, x, PHI(it, :));
        title(['t = ' num2str(round(t(it), 2))]);
    end
end

%% Function del RHS del sistema di ODEs
function dydt = myRHS(~, y, M, q)
    y = y(:);
    q = q(:);
    dydt = M*y + q;
end

%% Function per applicare il RK method
function y_new = RKsolve(t, y_old, M, q, dt, b, A)
    y_old = y_old(:);
    N = length(y_old);
    s = length(b);
    c = sum(A, 2);

    F = zeros(N, s);
    for j = 1:s
        F(:, j) = myRHS(t+dt*c(j), y_old+dt*F*A(j, :).', M, q);
    end

    y_new = y_old + dt*F*b(:);
end