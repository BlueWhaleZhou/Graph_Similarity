function [score_matrix, sim_derivation]= Score_Cal(A, B)



l1 = size(A, 1);
l2 = size(B, 1);

N = l1 + l2;
score_matrix = zeros(N, N);
sim_derivation = zeros(N, N);
%[A_norm_1, A_norm_2] = normalize(A);
%[B_norm_1, B_norm_2] = normalize(B);
%sim_0_t = GS_RW_Plain(A, B);
sim_0_t_norm = norm_kernel(A, B);%cosine similarity/graph kernel
%temp1 = GS_RW_Plain(A_norm_1, B_norm_1);
%disp("after normalization1: " + temp1(1));
%temp2 = GS_RW_Plain(A_norm_2, B_norm_2);
%disp("after normalization2: " + temp2(1));

%sim_0 = sim_0_t(1);%result of direct inverse random walk
%disp("no normalization " + sim_0);
 
%derivation part
lamda1 = eigs(A, 1, 'lm');
lamda2 = eigs(B, 1, 'lm');
lamda_new = 1 / (abs(lamda1 * lamda2) + 1);
init_dist = {ones(l1,1)/l1,ones(l2,1)/l2};
stop_dist = {ones(l1,1)/l1,ones(l2,1)/l2};
init1 = init_dist{1};
init2 = init_dist{2};
stop1 = stop_dist{1};
stop2 = stop_dist{2};
%disp(init1);
%disp(stop1);
init_x = kron(init1, init2);%kronecker product
stop_x = kron(stop1, stop2);
k_product = kron(A, B);
inv_part = inv(eye(l1 * l2) - lamda_new * k_product);

%if mod(l1, 2) == 1
%    l1_new = floor(l1 / 2) + 1;
%else
%    l1_new = l1 / 2;
%end

%if mod(l2, 2) == 1
%    l2_new = floor(l2 / 2) + 1;
%else
%    l2_new = l2 / 2;
%end
t0 = cputime;

for i2 = 1:l1
    for j = 1:l1
        A_tmp = A;
        if A(i2, j) == 0
            if i2 ~= j
                A_tmp(i2, j) = 1;%0 -> 0.01
                A_tmp(j, i2) = 1;
                score_normalized = norm_kernel(A_tmp, B);
                score_matrix(i2, j) = score_normalized - sim_0_t_norm;
                score_matrix(j, i2) = score_normalized - sim_0_t_norm;
        %        %A_tmp_norm = normalize(A_tmp);
        %        score_t = GS_RW_Plain(A_tmp, B);
        %        score = score_t(1);
        %        score_matrix(i2, j) = score - sim_0;
        %        score_matrix(j, i2) = score - sim_0;
            end
        end
        %else
        if A(i2, j) == 1
            if i2 ~= j
                A_tmp(i2, j) = 0;% 1 -> 0.99
                A_tmp(j, i2) = 0;
                score_normalized = norm_kernel(A_tmp, B);
                score_matrix(i2, j) = score_normalized - sim_0_t_norm;
                score_matrix(j, i2) = score_normalized - sim_0_t_norm;
                %A_tmp_norm = normalize(A_tmp);
                %score_t = GS_RW_Plain(A_tmp, B);
                %score = score_t(1);
                %score_matrix(i2, j) = score - sim_0;
                %score_matrix(j, i2) = score - sim_0;
            end
        end
    end
end

for i1 = 1:l2
    for j = 1:l2
        B_tmp = B;
        if B(i1, j) == 0
            if i1 ~= j
                B_tmp(i1, j) = 1; % 0 -> 0.01
                B_tmp(j, i1) = 1;
                score_normalize = norm_kernel(A, B_tmp);
                score_matrix(i1 + l1, j + l1) = score_normalize - sim_0_t_norm;
                score_matrix(j + l1, i1 + l1) = score_normalize - sim_0_t_norm;
        %        %B_tmp_norm = normalize(B_tmp);
        %        score_t = GS_RW_Plain(A, B_tmp);
        %        score = score_t(1);
        %        score_matrix(i1 + l1, j + l1) = score - sim_0;
        %        score_matrix(j + l1, i1 + l1) = score - sim_0;
            end
        end
        %else
        if B(i1, j) == 1
            if i1 ~= j
                B_tmp(i1, j) = 0; % 1 -> 0.99
                B_tmp(j, i1) = 0;
                score_normalize = norm_kernel(A, B_tmp);
                score_matrix(i1 + l1, j + l1) = score_normalize - sim_0_t_norm;
                score_matrix(j + l1, i1 + l1) = score_normalize - sim_0_t_norm;
                %B_tmp_norm = normalize(B_tmp);
                %score_t = GS_RW_Plain(A, B_tmp);
                %score = score_t(1);
                %score_matrix(i1 + l1, j + l1) = score - sim_0;
                %score_matrix(j + l1, i1 + l1) = score - sim_0;
            end
        end
    end
end
disp("score processing time is " + (cputime -t0));

t1 = cputime;
for row_no = 1:l1
    for col_no = 1:l1
        if row_no ~= col_no
            A_tmp = zeros(l1, l1);
            A_tmp(row_no, col_no) = 1;
            A_tmp(col_no, row_no) = 1;
        %disp(A_tmp);
            derivation_tmp = stop_x' * inv_part * (lamda_new * kron(A_tmp, B)) * inv_part * init_x;
            sim_derivation(row_no, col_no) = derivation_tmp;
        end
    end
end
for row_no = 1:l2
    for col_no = 1:l2
        if row_no ~= col_no
            B_tmp = zeros(l2, l2);
            B_tmp(row_no, col_no) = 1;
            B_tmp(col_no, row_no) = 1;
        %disp(B_tmp);
            derivation_tmp = stop_x' * inv_part * (lamda_new * kron(A, B_tmp)) * inv_part * init_x;
            sim_derivation(row_no + l1, col_no + l1) = derivation_tmp;
        end
    end
end

disp("deirvative processing time is " + (cputime - t1));

end
function norm_value = norm_kernel(A, B)
    A_B = GS_RW_Plain(A, B);
    A_B_1 = A_B(1);
    A_A = GS_RW_Plain(A, A);
    A_A_1 = A_A(1);
    B_B = GS_RW_Plain(B, B);
    B_B_1 = B_B(1);
    norm_value = A_B_1/(sqrt(A_A_1) * sqrt(B_B_1));
end
function sim = GS_RW_Plain(A,B,c,flag,p,q)
% % % % using random walk kernel to compute the similarity between A & B
n1 = size(A,1);
n2 = size(B,1);

if nargin<4
    flag = [1 1 1 1];%directly inverse;linear system;sylvester;eigens
end
if nargin<6
    q = {ones(n1,1)/n1,ones(n2,1)/n2};
end
if nargin<5
    p = {ones(n1,1)/n1,ones(n2,1)/n2};
end
if nargin<3
    [u1,lam1] = eigs(A,1,'lm');
    [u2,lam2] = eigs(B,1,'lm');
    c = 1/(abs(lam1 * lam2)+1);
end
p1 = p{1};
p2 = p{2};
q1 = q{1};
q2 = q{2};
   
if flag(1)==1%direct inverse
    X = kron(A,B);
    qx = kron(q1,q2);
    px = kron(p1,p2);
    sim(1,1) = qx' * inv(eye(n1 * n2) - c * X) * px;
end
if flag(2)==1%linear system
    X = kron(A,B);
    qx = kron(q1,q2);
    px = kron(p1,p2);
    x = rand(n1*n2,1);
    for i=1:200
        x = px + c * (X * x);
    end
    sim(1,2) = qx' * x;
end

if flag(4)==1%eigen-decomposition, for symmetric matrix only
    r1 = rank(A);
    r2 = rank(B);
    [U1,Lam1] = eigs(A,r1,'lm');
    [U2,Lam2] = eigs(B,r2,'lm');
    
    x1 = q1' * U1;
    x2 = q2' * U2;
    y1 = U1' * p1;
    y2 = U2' * p2;
    
    l = 1./(1-c*kron(diag(Lam1),diag(Lam2)));
    x = kron(x1,x2);
    y = kron(y1,y2);
    sim(1,4) = sum(x'.*l.*y);
%     L1 = diag(kron(Lam1,Lam2));
%     
%     H2 = diag(1./(1-c*L1));
%     sim0 = kron(x1,x2) * H2 * kron(y1,y2);
%     disp(sim0-sim(1,4))
end    
end

function [norm_matrix_1, norm_matrix_2] = normalize(A)

length = size(A, 1);
D = zeros(length, length);
for i = 1:length
    deg = 0;
    for j = 1:length
        deg = deg + A(i, j);
    end
    %disp(deg);
    D(i, i) = deg;
end%degree matrix 1

norm_matrix_1 = A * D^(-1);
norm_matrix_2 = D^(-1) * A;
end