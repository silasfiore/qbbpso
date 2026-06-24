function n = randuvec()

    while true

        u = 2*rand() -1;
        v = 2*rand() -1;

        if (u*u + v*v) < 1
            break
        end
    end

    w = u * u + v * v;

    x = 2 * u * sqrt(1 - w);
    y = 2 * v * sqrt(1 - w);
    z = 1 - 2 * w;
    n = [x;y;z];
end