close all; clc; clear;

%% Primo quesito
% Mesh spaziale
N = 15;
L = 1;
x = linspace(0, L, N); y = flip(x);
h = x(2) - x(1);
[X, Y] = meshgrid(x, y);

% BCs
zetWNE = 0; % lati Ovest, Nord, Est
zetS = 1; % lato Sud

% Matrici per la discretizzazione
G = numgrid('S', N);
figure('Name', 'Quesiti 1&2: Numgrid', 'NumberTitle', 'off'); spy(G); title('Discretizzazione');
L = -delsq(G) ./ h^2; L = full(L);

% Termine noto (omogeneo)
q = zeros(N-2); % forma matriciale

% Implementazione delle BCs nel termine noto
q(end, :) = q(end, :) - zetS/h^2;
q = q(:);

% Soluzione nei punti interni
zet = L \ q; % forma vettoriale
ZET = reshape(zet, N-2, N-2); % forma matriciale

% Orlare ZET per inserire le BCs
o = zeros(N-2, 1);
ZET = [o ZET o];
o = [0; o; 0].';
ZET = [o; ZET; zetS*ones(1, N)];

% Diagramma della soluzione
figure('Name', 'Quesito 1: Surface plot', 'NumberTitle', 'off');
surf(X, Y, ZET);
shading interp;
axis square;
xlabel('x'); ylabel('y'); zlabel('\zeta');
title('Soluzione');
hold on;
levels = linspace(min(zet), max(zet), 10);
contour3(X, Y, ZET, levels, 'k-');
figure('Name', 'Quesito 1: 2D Contour plot', 'NumberTitle', 'off');
contour(X, Y, ZET, levels);
axis square;
xlabel('x'); ylabel('y');
title('Curve di isolivello');

%% Secondo quesito
% BCs Dirichlet omogenee su tutti i bordi

% termine noto (soluzione del quesito precedente)
q = ZET(2:end-1, 2:end-1);
q = q(:); % non vanno implementate le BCs poiché tutte Dirichlet omogenee

% soluzione nei punti interni
psi = L \ q;
PSI = reshape(psi, N-2, N-2);

% Orlare PSI per inserire le BCs
o = zeros(N-2, 1);
PSI = [o PSI o];
o = [0; o; 0].';
PSI = [o; PSI; o];

% Diagramma della soluzione
figure('Name', 'Quesito 2: Surface plot', 'NumberTitle', 'off');
surf(X, Y, PSI);
shading interp;
axis square;
xlabel('x'); ylabel('y'); zlabel('\psi');
title('Soluzione');
hold on;
levels = linspace(min(psi), max(psi), 10);
contour3(X, Y, PSI, levels, 'k-');
figure('Name', 'Quesito 2: 2D Contour plot', 'NumberTitle', 'off');
contour(X, Y, PSI, levels);
axis square;
xlabel('x'); ylabel('y');
title('Curve di isolivello');

%% Terzo quesito
Hdomain(25, 2, 1.25, 2, 1);

A = [.15 .5 .75 1 1.15 1.5 2 2.5]; % A = Lx / Ly  |  DEVE AVERE LENGTH PARI;
Ly = ones(1, length(A)); Lx = Ly .* A;
zetB = 1;

figure('Name', 'Quesito 3: Studio topologico', 'NumberTitle', 'off');
for i = 1:length(A)
    [X, Y, PSI] = Hdomain(25, Lx(i), Ly(i), zetB, 0);
    levels = linspace(min(PSI, [], 'all'), max(PSI, [], 'all'), 8);
    subplot(2, length(A)/2, i);
    contour(X, Y, PSI, levels);
    axis square;
    xlabel('x'); ylabel('y');
    title(['Lx / Ly = ' num2str(A(i))]);
end

function [X, Y, PSI] = Hdomain(N, Lx, Ly, zetB, plt)
    % N = numero di punti di discretizzazione sul lato corto
    % Lx = lunghezza del lato lungo x
    % Ly = lunghezza del lato lungo y
    % zetB = condizione al contorno alla Dirichlet non omogenea
    % plt = variabile booleana (1 se si vuole il surf plot, 0 altrimenti)

    if Lx <= Ly
        Nx = N;
        h = Lx / Nx;
        Ny = floor(Ly / h);
    else
        Ny = N;
        h = Ly / Ny;
        Nx = floor(Lx / h);
    end

    x = linspace(0, Lx, Nx);
    y = linspace(0, Ly, Ny); y = flip(y);
    [X, Y] = meshgrid(x, y);

    G = zeros(Ny, Nx);
    count = 1;
    for j = 2:Nx-1
        for i = 2:Ny-1
            G(i, j) = count;
            count = count + 1;
        end
    end

    L = -delsq(G) ./ h^2; L = full(L);

    q = zeros(Ny-2, Nx-2);
    if Lx <= Ly
        q(end, :) = -zetB / h^2;
        q(1, :) = q(end, :);
    else
        q(:, 1) = -zetB / h^2;
        q(:, end) = q(:, 1);
    end
    q = q(:);
    
    zet = L \ q;
    ZET = reshape(zet, Ny-2, Nx-2);

    if Lx <= Ly
        o = zeros(Ny-2, 1);
        ZET = [o ZET o];
        b = zetB * ones(1, Nx);
        ZET = [b; ZET; b];
    else
        o = zeros(1, Nx-2);
        ZET = [o; ZET; o];
        b = zetB * ones(Ny, 1);
        ZET = [b ZET b];
    end

    if plt == 1
        figure('Name', 'Quesito 3: Zeta', 'NumberTitle', 'off');
        surf(X, Y, ZET);
        shading interp;
        axis square;
        xlabel('x'); ylabel('y'); zlabel('\zeta');
        title('Soluzione');
        hold on;
        levels = linspace(min(zet), max(zet), 10);
        contour3(X, Y, ZET, levels, 'k-');
    end

    q = ZET(2:end-1, 2:end-1);
    q = q(:);

    psi = L \ q;
    PSI = reshape(psi, Ny-2, Nx-2);
    o = zeros(Ny-2, 1);
    PSI = [o PSI o];
    o = zeros(1, Nx);
    PSI = [o; PSI; o];

    if plt == 1
        figure('Name', 'Quesito 3: Psi', 'NumberTitle', 'off');
        surf(X, Y, PSI);
        shading interp;
        axis square;
        xlabel('x'); ylabel('y'); zlabel('\psi');
        title('Soluzione');
        hold on;
        levels = linspace(min(psi), max(psi), 10);
        contour3(X, Y, PSI, levels, 'k-');
    end
end