close all; clc; clear

% Mesh spaziale
lungh = 1;
Ny = 35;
y = linspace(0, lungh, Ny);
h = y(2) - y(1);
x = 0:h:2*lungh;
Nx = length(x);
[X, Y] = meshgrid(x, y);

% BCs
phiE = 1;
phiW = 0;
DphiN = 1;
DphiS = 0;

% Numgrid
m = Ny - 2;
n = Nx - 2;
nz = m * n;
G = zeros(m, n);
G(G==0) = 1:nz;
G = [zeros(1, n); G; zeros(1, n)];
G = [zeros(Ny, 1) G zeros(Ny, 1)];

% Termine noto
q = ones(nz, 1);

% Operatore delsq
L = -delsq(G) ./ h^2;

% BCs alla Neumann in L
kS = G(2, 2:end-1);
kN = G(end-1, 2:end-1);
for j = 1:n
    L(kS(j), kS(j)) = L(kS(j), kS(j)) + 1/h^2;
    L(kN(j), kN(j)) = L(kN(j), kN(j)) + 1/h^2;
end

% BCs alla Dirichlet in q
kE = G(2:end-1, end-1);
kW = G(2:end-1, 2);
for i = 1:m
    q(kE(i)) = q(kE(i)) - phiE/h^2;
    q(kW(i)) = q(kW(i)) - phiW/h^2;
end

% BCs alla Neumann in q
for j = 1:n
    q(kN(j)) = q(kN(j)) - DphiN/h;
    q(kS(i)) = q(kS(i)) + DphiS/h;
end

% Soluzione
phi = L \ q;
PHI = reshape(phi, m, n);
phiN = PHI(end, :) + DphiN*h;
phiS = PHI(1, :) - DphiS*h;
PHI = [phiS; PHI; phiN];
phiE = phiE * ones(Ny, 1);
phiW = phiW * ones(Ny, 1);
PHI = [phiW PHI phiE];

% Diagramma
figure
surf(X, Y, PHI)
shading interp
xlabel('x'); ylabel('y'); zlabel('\phi')
grid on
hold on
contour3(X, Y, PHI, linspace(min(PHI(:)), max(PHI(:)), 25), 'k-')