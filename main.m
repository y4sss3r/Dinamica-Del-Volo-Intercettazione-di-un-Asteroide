figure('Color', 'k', 'Name', 'Simulatore Orbitale');
hold on; grid off; axis equal;
view(3);
lunghezza_vettori=10000;
versore_z=[0, 0, 1];
versore_y=[0, 1, 0];
versore_x=[1, 0, 0];

raggio_terra = 6371;
[xe, ye, ze] = sphere(50);
topo = load('topo.mat');
terra = surface(xe*raggio_terra, ye*raggio_terra, ze*raggio_terra, 'FaceColor', 'texturemap', 'CData', topo.topo, 'EdgeColor', 'none', HandleVisibility='off');
colormap(demcmap(topo.topo));

dt=0.1
propagationTime=8497;
mezzora=30*60/dt;
[X, Y, Z, T, MAG_V, MAG_R] = get_orbit(dt, ceil(propagationTime/dt));


init_position = [-7368.038574853538, -7231.584293256432, -148.523707822187];
init_velocity = [4.126512186315761, -3.956371322777358, -0.490613661500991];

versore_init_velocity=init_velocity/norm(init_velocity)*lunghezza_vettori;
momento_specifico=cross(init_position, init_velocity);
versore_h=momento_specifico/norm(momento_specifico)*lunghezza_vettori;

posizioni_satellite=[];

orbit=plot3(X, Y, Z, 'w-', 'LineWidth', 1.5);
for point=0:4
    posizione=[X(point*mezzora+1), Y(point*mezzora+1), Z(point*mezzora+1)];
    scatter=plot3(posizione(1), posizione(2), posizione(3), "ro", LineWidth=2);
    posizioni_satellite=[posizioni_satellite; posizione];
end


quiver3(0, 0, raggio_terra, 0, 0, 1*lunghezza_vettori, 0, 'LineWidth', 2, 'Color', 'r', 'MaxHeadSize', 0.5, HandleVisibility='off');
quiver3(0, raggio_terra, 0, 0, 1*lunghezza_vettori, 0, 0, 'LineWidth', 2, 'Color', 'g', 'MaxHeadSize', 0.5, HandleVisibility='off');
quiver3(raggio_terra, 0, 0, 1*lunghezza_vettori, 0, 0, 0, 'LineWidth', 2, 'Color', 'b', 'MaxHeadSize', 0.5, HandleVisibility='off');
h_vec=quiver3(0, 0, raggio_terra, versore_h(1), versore_h(2), versore_h(3), 0, 'LineWidth', 2, 'Color', [1 153/255 1], 'MaxHeadSize', 0.5);
vi_vec=quiver3(init_position(1), init_position(2), init_position(3), versore_init_velocity(1), versore_init_velocity(2), versore_init_velocity(3), 0, 'LineWidth', 2, 'Color', [1 1 0.2], 'MaxHeadSize', 0.5);
xlabel("X (Km)");
ylabel('Y (km)');
zlabel('Z (km)');
title('Orbital Trajectory');
legend([orbit, scatter, h_vec, vi_vec], ["orbita", "posizione ogni mezz'ora", "momento specifico", "velocità iniziale"])



hold off;

