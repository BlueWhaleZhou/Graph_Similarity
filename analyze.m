load DBLP;

%smola 
%currentTeam = [916232, 250177, 219532, 545532, 756929]; %paper index
%i0= 250177;

%hinton_1: 
%currentTeam = [208818, 554802, 455386, 218580, 822000]; %paper index
%i0= 218580;

%jiawei han_1
currentTeam = [750967, 124372, 167178, 620387, 3603, 582475];
i0 = 124372;

%jiawei han_2
%currentTeam = [37832, 785879, 721587, 124372, 260647];
%i0 = 124372;

%C. Faloutsos_1
%currentTeam = [603161, 175073, 206733, 185952, 314018];
%i0 = 185952;

%Philip S. Yu_1
%currentTeam = [414491, 722335, 912598, 124372, 247326];
%i0 = 247326;

%Philip S. Yu_2
%currentTeam = [399690, 216073, 431073, 247326, 551591];
%i0 = 247326;

%decay factor
c = 0.00000001;

dn = 43; % number of skills
L = cell(1,dn);
for i = 1:dn 
    L{i} = diag(count_label(:,i));
end

fileID=fopen('authorDict.txt');
authorDict=textscan(fileID,'%s','delimiter','\n');
authorDict=authorDict{1};

fprintf('We need to replace %s ...\n', authorDict{i0});

score = label_direct_recommend(aa,L,currentTeam,i0,true);
top6 = topsix(score);
disp ('Using TEAMREP-BASIC after pruning, the top six candidates are:');

A1 = aa(currentTeam, currentTeam);
A1 = (triu(A1, 1) + tril(A1, -1));

existingTeam = setdiff(currentTeam, i0, 'stable');
disp("remained team members...")
disp(existingTeam);
top_cand_len = length(top6);

abs_old_cell = cell(1, top_cand_len);
abs_new_cell = cell(1, top_cand_len);
inf_old_cell = cell(1, top_cand_len);
inf_new_cell = cell(1, top_cand_len);

abs_time = zeros(1, top_cand_len);
inf_time = zeros(1, top_cand_len);

for i = 1:top_cand_len
    newTeam = [existingTeam, top6(i)];
    disp(newTeam);
    LL = zeros(length(newTeam)^2);
    for j = 1:dn
        LL = LL + kron(L{j}(currentTeam, currentTeam), L{j}(newTeam, newTeam));
    end
    A2 = aa(newTeam, newTeam);
    A2 = (triu(A2, 1) + tril(A2, -1));
    
    %abs_cell filling
    [tmp_abs_old, tmp_abs_new, tmp_abs_time] = abs_sim(A1, A2, LL, c);
    abs_old_cell{i} = tmp_abs_old;
    abs_new_cell{i} = tmp_abs_new;
    abs_time(i) = tmp_abs_time;
    
    %inf_cell filling
    [tmp_inf_old, tmp_inf_new, tmp_inf_time] = inf_matrix_cal(A1, A2, LL, c);
    inf_old_cell{i} = tmp_inf_old;
    inf_new_cell{i} = tmp_inf_new;
    inf_time(i) = tmp_inf_time;
    
end

fprintf('%s \n', authorDict{top6});