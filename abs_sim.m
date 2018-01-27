function [abs_old_team, abs_new_team, abs_time] = abs_sim(A,B,L,c)
    
    start = cputime;
    
    n1 = size(A,1);
    n2 = size(B,1);   
    abs_old_team = zeros(n1);
    abs_new_team = zeros(n2);
    %old similarity before a memeber in B is replaced
    sim_old = label_gs(A,B,L,c);
    for i = 1:n1
        for j = 1:n1
            tmp = A;
            if i ~= j
                tmp(i, j) = 0;
                tmp(j, i) = 0;
                sim_new = label_gs(tmp,B,L,c);
                abs_old_team(i, j) = sim_new - sim_old;
                abs_old_team(j, i) = sim_new - sim_old;
            end
        end
    end
    
    for i = 1:n2
        for j = 1:n2
            tmp = B;
            if i ~= j
                tmp(i, j) = 0;
                tmp(j, i) = 0;
                sim_new = label_gs(A,tmp,L,c);
                abs_new_team(i, j) = sim_new - sim_old;
                abs_new_team(j, i) = sim_new - sim_old;
            end
        end
    end
    
    abs_time = cputime - start;
end
