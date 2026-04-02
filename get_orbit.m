function [X, Y, Z, T, MAG_V, MAG_R] = get_orbit(dt, itermax)
    mu = 398600;
    position = [-7368.038574853538, -7231.584293256432, -148.523707822187];
    velocity = [4.126512186315761, -3.956371322777358, -0.490613661500991];
    X=[-7368.038574853538];
    Y=[-7231.584293256432];
    Z=[-148.523707822187];
    T=[0];
    MAG_V=[norm(velocity)];
    MAG_R=[norm(position)];
    for iter=1:itermax
        fprintf("status: %.2f\n", iter/itermax*100);
        acc = - (mu / norm(position)^3) * position;
        velocity = velocity + acc * dt;
        MAG_V=[MAG_V, norm(velocity)];
        position = position + velocity * dt;
        MAG_R=[MAG_R, norm(position)];
        X=[X, position(1)];
        Y=[Y, position(2)];
        Z=[Z, position(3)];
        T=[T, iter*dt];
    end
end
