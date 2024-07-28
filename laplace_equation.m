close all; clc; clear

%% Quesito 1
% Mesh spaziale
L = 1;
N = 100;
x = linspace(0, L, N);
y = x;
h = x(2) - x(1);

% BCs
phiE = cos(2*pi/L * y); phiE = flip(phiE);
phiN = x ./ L;

% Numgrid
G = numgrid('S', N);
nz = max(G(:));

% Operatore delsq
L = -delsq(G) ./ h^2;

% Termine noto
q = zeros(nz, 1);

% BCs alla Neumann in L
k = G(2, 2:end-1);
for i = 1:N-2
    L(k(i), k(i)) = L(k(i), k(i)) + 1/h^2; % Lato Sud
    L(i, i) = L(i, i) + 1/h^2; % Lato Ovest
end

% BCs alla Dirichlet nel termine noto
kN = G(end-1, 2:end-1);
kE = G(2:end-1, end-1);
for i = 1:N-2
    q(kN(i)) = q(kN(i)) - phiN(i+1)/h^2;
    q(kE(i)) = q(kE(i)) - phiE(i+1)/h^2;
end

% Soluzione
phi = L \ q;
PHI = reshape(phi, N-2, N-2);
phiN = phiN(2:end-1);
PHI = [PHI(1, :); PHI; phiN(:).'];
PHI = [PHI(:, 1) PHI phiE(:)];

% Diagramma
[X, Y] = meshgrid(x, y);
figure
surf(X, Y, PHI);
xlabel('x'); ylabel('y'); zlabel('\phi')
shading interp
hold on
contour3(X, Y, PHI, linspace(min(PHI(:)), max(PHI(:)), 25), 'k-')