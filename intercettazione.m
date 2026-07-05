% --- MISSIONE ---
% 1. fuga dalla soi terrestre
% 2. inclinazione dell'orbita in afelio
% 3. hohmann
% 4. phasing

clear; clc;

% --- DATI INIZIALI ---
sat_position = [-7368.038574853538, -7231.584293256432, -148.523707822187];
sat_velocity = [4.126512186315761, -3.956371322777358, -0.490613661500991];

earth_position = [-1.211877736660399E+08, -8.932426658555530E+07, 6.074880469534546E+03];
earth_velocity = [1.719435858581588E+01, -2.407934077970530E+01, 2.444511952605311E-03];

asteroide_position = [5.7220e+7, 2.1137e+8, 1.4634e+8];
asteoide_velocity = [-19.9150495, 3.1435630, 4.3690766];

%lo stato deve essere un vettore COLONNA (18x1)
sat_state0 = [sat_position(:); sat_velocity(:)]; 
earth_state0 = [earth_position(:); earth_velocity(:)];
asteroide_state0 = [asteroide_position(:); asteoide_velocity(:)];

state0=[earth_state0; sat_state0; asteroide_state0];

consts.mu_terra = 398600;
consts.mu_sole = 1.32712440018e11;
au=149597870.7; % 1 Unità Astronomica in km
options = odeset('RelTol', 1e-8, 'AbsTol', 1e-10);



launch_date_init=300*24*60*60;%3.6314e+07 - 14*24*60*60;
launch_date_max=500*24*60*60; %3.6314e+07 + 14*24*60*60;
min_tof=50*24*60*60;
max_tof=200*24*60*60;
samples=100;
date_lancio= linspace(launch_date_init, launch_date_max, samples);
tof = linspace(min_tof, max_tof, samples);
[T_lancio, TOF]=meshgrid(date_lancio, tof);
c3_map = zeros(size(T_lancio));


t_init=0;
t_propagation = launch_date_max+max_tof;



state_sol=ode45(@(t, y) solver(t, y, consts), [t_init, t_propagation], state0, options);

for i=1:size(T_lancio, 1)
    for j=1:size(T_lancio, 2)
        t_partenza=T_lancio(i, j);
        tof=TOF(i, j);
        t_arrivo=t_partenza+tof;
        earth_position = deval(state_sol, t_partenza, 1:3);
        earth_velocity = deval(state_sol, t_partenza, 4:6);
        sat_position = deval(state_sol, t_partenza, 7:9) + earth_position;
        sat_velocity_rel = deval(state_sol, t_partenza, 10:12);
        sat_velocity=sat_velocity_rel+earth_velocity;
        ast_position = deval(state_sol, t_arrivo, 13:15);
        [deltaV, v_transfer] = Lambert(earth_position, earth_velocity, ast_position, tof, consts.mu_sole, 0);

        v_inf=norm(v_transfer - earth_velocity);
        c3_map(i, j) = v_inf^2;
        V_TRANS_MAP(i, j, :)=v_transfer;

    end
end

[minC3, linear_index] = min(c3_map(:));
fprintf('min C3 calcolato: %.2f Km^2/s^2\n', minC3);
[i_min, j_min] = ind2sub(size(c3_map), linear_index);

fprintf("data di lancio trovata: giorno %d\n", T_lancio(i_min, j_min)/(60*60*24));
fprintf("durata di viaggio: %d giorni\n", TOF(i_min, j_min)/(60*60*24));

figure('Color', 'k', 'Name', 'Porkchop Plot');
livelli = [0, 55, 60, 65, 70, 75, 80, 85, 100, 120, 140, 150];
[C, h] = contourf(T_lancio/(24*3600), TOF/(24*3600), c3_map, livelli);
clabel(C, h);
colorbar;
xlabel('Giorni dalla partenza');
ylabel('Durata del viaggio (Giorni)');
title('Porkchop Plot: Energia di Lancio C3 [km^2/s^2]');
colormap('jet');


% --- INTEGRAZIONE ---
% vettore stato::  
%   1:3 posizione terra
%   4:6 velocità terra
%   7:9 posizione satellite (rel. terra)
%   10:12 velocita satellite (rel. terra)
%   13:15 posizione asteroide (ass, sole)
%   16:18 velocita asteroide (ass. sole)


v_transfer=[V_TRANS_MAP(i_min, j_min, 1), V_TRANS_MAP(i_min, j_min, 2), V_TRANS_MAP(i_min, j_min, 3)];
manovra_1.t=T_lancio(i_min, j_min);
manovra_1.get_deltaV = @(state) utils.get_manovra_delta_velocity(state, v_transfer);

manovre=[manovra_1];

ast_position=deval(state_sol, T_lancio(i_min, j_min)+TOF(i_min, j_min), 13:15);
correction_manouvres=10;


%{
frazione = TOF(i_min, j_min)/(correction_manouvres+1);
for correction=1:correction_manouvres
    t_correction = T_lancio(i_min, j_min)+frazione*correction;
    tof_rimasto=TOF(i_min, j_min)-frazione*correction;
    correction_manouvre.t=t_correction;
    correction_manouvre.get_deltaV = @(state) utils.get_correction_deltaV(state, ast_position, tof_rimasto);
    manovre=[manovre, correction_manouvre];
end
%}

make_callback = @(ast_pos, tof) @(state) utils.get_correction_deltaV(state, ast_pos, tof);
for i=1:correction_manouvres
    frazione_tof=TOF(i_min, j_min)*(1-0.5^i );
    t_correction= T_lancio(i_min, j_min)+frazione_tof;
    tof_rimasto=TOF(i_min, j_min)-frazione_tof;
    correction_manouvre.t=t_correction;
    correction_manouvre.get_deltaV = make_callback(ast_position, tof_rimasto);
    manovre=[manovre, correction_manouvre];
end



T_tot=[];
state_sol=[];
state0=[earth_state0; sat_state0; asteroide_state0];
t_propagation =T_lancio(i_min, j_min)+TOF(i_min, j_min); %T_lancio(i_min, j_min)+ TOF(i_min, j_min)+22000;
deltaV_tot=0;
for manovra=manovre
    t_i=[t_init, manovra.t];
    [T_i, state_i] = ode45(@(t, y) solver(t, y, consts), t_i, state0, options);
    
    if isempty(T_tot)
        T_tot = T_i;
        state_sol = state_i;
    else
        T_tot = [T_tot; T_i(2:end)];
        state_sol = [state_sol; state_i(2:end, :)];
    end

    deltaV_vec=manovra.get_deltaV(state_i(end, :)); % vettore
    fprintf("manovra eseguita con deltaV: %d\n", norm(deltaV_vec))
    deltaV_tot=deltaV_tot+norm(deltaV_vec);
    v_sat=state_i(end, 10:12);
    v_sat_new=v_sat+deltaV_vec;
    state_i(end, 10:12)=v_sat_new;
    t_init=t_i(end);
    state0=state_i(end, :)'; 
end
[Tf, state_sol_f] = ode45(@(t, y) solver(t, y, consts), [t_init, t_propagation], state0, options);
if isempty(T_tot)
    T_tot = Tf;
    state_sol = state_sol_f;
else
    T_tot = [T_tot; Tf(2:end)];
    state_sol = [state_sol; state_sol_f(2:end, :)];
end

fprintf("deltaV totale: %d Km/s\n", deltaV_tot);

pos_terra = state_sol(:, 1:3);
vel_terra= state_sol(:, 4:6);
pos_sat_rel = state_sol(:, 7:9);
vel_sat_rel = state_sol(:, 10:12);
pos_sat_assoluta = pos_terra + pos_sat_rel;
vel_sat_assoluta = vel_terra + vel_sat_rel;
pos_asteroide = state_sol(:, 13:15);
vel_asteroide = state_sol(:, 16:18);

r_sat_rel=vecnorm(pos_sat_rel, 2, 2);
v_sat_rel=vecnorm(state_sol(:, 10:12), 2, 2);
r_sat_ass=vecnorm(pos_sat_assoluta, 2, 2);
v_sat_ass=vecnorm(vel_sat_assoluta, 2, 2);

distanze = vecnorm(pos_sat_assoluta - pos_asteroide, 2, 2);
[min_dist, idx_min] = min(distanze);
fprintf('Distanza minima allincontro: %.2f km (%.4f AU)\n', min_dist, min_dist/au);

% --- VISUALIZZAZIONE GRAFICO VELOCITA E TEMPO (sistema terrestre) ---
%figure('Color', 'k', 'Name', 'r(t) e v(t) assoluta');
%subplot(2, 1, 1);
%plot(T_tot, r_sat_ass);
%title("r(t)"); xlabel("T (s)"); ylabel("r (km)");
%subplot(2, 1, 2);
%plot(T_tot, v_sat_ass)
%title("v(t)"); xlabel("T (s)"); ylabel("v (km/s)");

%figure('Color', 'k', 'Name', 'r(t) e v(t) relativa');
%subplot(2, 1, 1);
%plot(T_tot, r_sat_rel);
%title("r(t)"); xlabel("T (s)"); ylabel("r (km)");
%subplot(2, 1, 2);
%plot(T_tot, v_sat_rel)
%title("v(t)"); xlabel("T (s)"); ylabel("v (km/s)");

% --- VISUALIZZAZIONE (sistema solare) ---
figure('Color', 'k', 'Name', 'sistema solare');
hold on; axis equal; grid off;
set(gca, 'Color', 'k', 'XColor', 'w', 'YColor', 'w', 'ZColor', 'w');
view(3);
raggio_sole = 0.2;
[xe, ye, ze] = sphere(50);
surface(xe*raggio_sole, ye*raggio_sole, ze*raggio_sole)

plot3(pos_terra(1, 1)/au, pos_terra(1, 2)/au, pos_terra(1, 3)/au, 'ro');
plot3(pos_terra(:, 1)/au, pos_terra(:, 2)/au, pos_terra(:, 3)/au, 'yellow', 'LineWidth', 1.5);
plot3(pos_sat_assoluta(:, 1)/au, pos_sat_assoluta(:, 2)/au, pos_sat_assoluta(:, 3)/au, "cyan", 'LineWidth', 1.5)
plot3(pos_asteroide(1, 1)/au, pos_asteroide(1, 2)/au, pos_asteroide(1, 3)/au, "ro")
plot3(pos_asteroide(:, 1)/au, pos_asteroide(:, 2)/au, pos_asteroide(:, 3)/au, "green", 'LineWidth', 1.5)

% --- VISUALIZZAZIONE (sistema terrestre) ---
figure('Color', 'k', 'Name', 'sistema terrestre');
hold on; axis equal; grid off;
set(gca, 'Color', 'k', 'XColor', 'w', 'YColor', 'w', 'ZColor', 'w');
view(3);

raggio_terra = 6371;
[xe, ye, ze] = sphere(50);

topo = load('topo.mat');
surface(xe*raggio_terra, ye*raggio_terra, ze*raggio_terra, ...
    'FaceColor', 'texturemap', 'CData', topo.topo, 'EdgeColor', 'none');
colormap(demcmap(topo.topo));

%plot3(-pos_terra(:, 1), -pos_terra(:, 2), -pos_terra(:, 3), "yellow");
plot3(pos_sat_rel(:, 1), pos_sat_rel(:, 2), pos_sat_rel(:, 3), 'cyan', 'LineWidth', 1.5);
plot3(pos_sat_rel(1,1), pos_sat_rel(1,2), pos_sat_rel(1,3), 'ro', 'MarkerFaceColor', 'r');
function dY = solver(~, Y, consts) % Y vettore stati 1x18
    % --- STATO TERRA (Rispetto al Sole) ---
    r_terra = Y(1:3);
    v_terra = Y(4:6);
    dist_terra = norm(r_terra);
    acc_terra = -(consts.mu_sole / dist_terra^3) * r_terra;

    % --- STATO SATELLITE (Rispetto alla Terra) ---
    r_sat_rel = Y(7:9);
    v_sat_rel = Y(10:12);
    dist_sat = norm(r_sat_rel);

    r_sat_ass=r_terra+r_sat_rel;
    dist_sat_ass=norm(r_sat_ass);

    acc_sat_terra = -(consts.mu_terra / dist_sat^3) * r_sat_rel;
    acc_sat_sole_perturb = -consts.mu_sole * ( (r_sat_ass / dist_sat_ass^3) - (r_terra / dist_terra^3) );
    acc_sat_rel = acc_sat_terra + acc_sat_sole_perturb;

    % --- STATO ASTEROIDE ---
    r_asteroide=Y(13:15);
    v_asteroide = Y(16:18);
    acc_asteroide = -(consts.mu_sole / norm(r_asteroide)^3) * r_asteroide;

    % --- OUTPUT ---
    dY = [v_terra; acc_terra; v_sat_rel; acc_sat_rel; v_asteroide; acc_asteroide];
end

function [] = porkchop_plot()
    state_sol=ode45(@(t, y) solver(t, y, consts), [t_init, t_propagation], state0, options);
    for i=1:size(T_lancio, 1)
        for j=1:size(T_lancio, 2)
            t_partenza=T_lancio(i, j);
            tof=TOF(i, j);
            t_arrivo=t_partenza+tof;
            earth_position = deval(state_sol, t_partenza, 1:3);
            earth_velocity = deval(state_sol, t_partenza, 4:6);
            ast_position = deval(state_sol, t_arrivo, 13:15);
            [~, v_transfer] = Lambert(earth_position, earth_velocity, ast_position, tof, consts.mu_sole, 0);
    
            v_inf=norm(v_transfer - earth_velocity);
            c3_map(i, j) = v_inf^2;
            V1_REQ_MAP(i, j, :)=v_transfer;
    
        end
    end
    
    minC3 = min(c3_map(:));
    fprintf('Range C3 calcolato: [%.2f, %.2f]\n', minC3, maxC3);
    
    figure('Name', 'Porkchop Plot');
    %livelli = linspace(minC3, minC3 + 100, 20);
    livelli = [0, 55, 60, 65, 70, 75, 80, 85, 100, 120, 140, 150]
    [C, h] = contourf(T_lancio/(24*3600), TOF/(24*3600), c3_map, livelli);
    clabel(C, h);
    colorbar;
    xlabel('Giorni dal riferimento (Partenza)');
    ylabel('Durata del viaggio (Giorni)');
    title('Porkchop Plot: Energia di Lancio C3 [km^2/s^2]');
    colormap('jet');

end
