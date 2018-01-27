%graph kernel influcial matrix score computation
function [inf_old_team, inf_new_team, inf_time] = inf_matrix_cal(A, B, L, c)
    
    start = cputime;
    
    n1 = size(A, 1);
    n2 = size(B,1);
    inf_old_team = zeros(n1);
    inf_new_team = zeros(n2);
    
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
    
    %get the matrix for every edge's influence score in old team
    for i = 1:n1
        for j = 1:n1
            if i ~= j
                tmp = zeros(n1);
                tmp(i, j) = 1;
                tmp(j, i) = 1;
                core = kron(tmp, B);           
                inf_old_team(i, j) = first_part * core * second_part;
                inf_old_team(j, i) = inf_old_team(i, j);
            end
        end
    end
    
    %get the matrix for every edge's influence score in new team
    for i = 1:n2
        for j = 1:n2
            if i ~= j
                tmp = zeros(n2);
                tmp(j, i) = 1;
                tmp(i, j) = 1;
                core = kron(A, tmp);
                inf_new_team(i, j) = first_part * core * second_part;
                inf_new_team(j, i) = inf_new_team(i, j);
            end
        end
    end
    
    inf_time = cputime - start;
end