function sim_der = inf_cal(A, B, L, c)
    n1 = size(A, 1);
    n2 = size(B, 1);
    sim_der = zeros(1, n2 - 1);
    
    q = {ones(n1,1)/n1,ones(n2,1)/n2};
    p = {ones(n1,1)/n1,ones(n2,1)/n2};
    p1 = p{1};
    p2 = p{2};
    q1 = q{1};
    q2 = q{2};
    qx = kron(q1, q2);
    px = kron(p1, p2);
    Wx = kron(A, B);
    inv_part = inv(eye(n1 * n2) - c * L * Wx);
    first_part = qx' * inv_part * c * L;
    second_part = inv_part * L * px;
    
    for i = 1 : (n2 - 1)
        tmp = zeros(n2);
        tmp(n2, i) = 1;
        tmp(i, n2) = 1;
        core = kron(A, tmp);
        sim_der(1, i) = first_part * core * second_part;
    end
    %disp("sim_der size " + size(sim_der));
end