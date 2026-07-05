function [deltaV, v1_trans] = Lambert(r1, v1_init, r2, dt, mu, lw)
% r1: Posizione iniziale [1x3]
% v1_init: Velocità attuale [1x3]
% r2: Posizione finale [1x3]
% dt: Tempo di volo desiderato
% mu: Parametro gravitazionale standard (es. Terra: 3.986e5)
% lw: Long Way flag (0 per arco corto < 180°, 1 per arco lungo > 180°)

r1_mag = norm(r1);
r2_mag = norm(r2);

% 1. Calcolo dell'angolo di trasferimento (cos theta)
cos_theta = dot(r1, r2) / (r1_mag * r2_mag);

if lw == 0
    theta = acos(cos_theta);
else
    theta = 2*pi - acos(cos_theta);
end

% 2. Costante geometrica A
A = sin(theta) * sqrt(r1_mag * r2_mag / (1 - cos_theta));

% 3. Risoluzione iterativa per la variabile universale 'z'
% Definiamo la funzione del tempo di volo in base a z
z = 0; % Valore iniziale
tol = 1e-8;
max_iter = 1000;

% Limiti per bisezione (z < (2*pi)^2 per orbite ellittiche)
z_up = 4 * pi^2;
z_low = -100; % Per orbite iperboliche

for i = 1:max_iter
    [C, S] = stumC_S(z);

    % Calcolo di y(z)
    y = r1_mag + r2_mag + A * (z * S - 1) / sqrt(C);

    % Controllo per evitare numeri immaginari con A e y
    if A > 0 && y < 0
        % Regolazione del limite inferiore se y diventa negativo
        z_low = z; 
    end

    chi = sqrt(y / C);
    dt_check = (chi^3 * S + A * sqrt(y)) / sqrt(mu);

    % Metodo della bisezione per trovare la radice
    if abs(dt_check - dt) < tol
        break;
    end

    if dt_check <= dt
        z_low = z;
    else
        z_up = z;
    end
    z = (z_up + z_low) / 2;
end

% 4. Calcolo dei coefficienti di Lagrange f, g, g_dot
f = 1 - y / r1_mag;
g = A * sqrt(y / mu);
g_dot = 1 - y / r2_mag;

% 5. Velocità di trasferimento al punto 1
v1_trans = (r2 - f * r1) / g;

% 6. Calcolo del Delta-V richiesto
deltaV = v1_trans - v1_init;

%fprintf('Iterazioni completate: %d\n', i);
%fprintf('Delta-V totale: %.4f km/s\n', norm(deltaV));
end

% Funzione ausiliaria: Funzioni di Stumpff C e S
function [C, S] = stumC_S(z)
if z > 0
    S = (sqrt(z) - sin(sqrt(z))) / (sqrt(z)^3);
    C = (1 - cos(sqrt(z))) / z;
elseif z < 0
    S = (sinh(sqrt(-z)) - sqrt(-z)) / (sqrt(-z)^3);
    C = (cosh(sqrt(-z)) - 1) / (-z);
else
    S = 1/6;
    C = 1/2;
end
end
