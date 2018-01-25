function [score, edge_weight_matrix] = label_direct_recommend(aa,L,currentTeam,i0,prune)
%The proposed TEAMREP-BASIC algorithm
%input:
%aa: the whole social network, e.g., author-author network, size nxn
%L: label matrice cell, e.g., if there are dn skills, then L is a cell of size dn, 
%L{i} is a nxn diagonal matrix, L{i}(j,j) shows the strength of j-th person
%having i-th skill
%
%currentTeam: current members of a team, e.g., authors for a paper
%i0: the member to be replaced
%prune: prune or not?
%
%output:
%score: each row is a score and its candidate id, note it's not sorted

if nargin < 4
    prune = false;
end

n=size(aa,1);
remainTeam = setdiff(currentTeam,i0); %i0 is the index of author to be replaced
currentTeam = [remainTeam, i0];
%disp(currentTeam);
%disp("check point 1");

A1 = aa(currentTeam,currentTeam);
A1 = (triu(A1,1) + tril(A1,-1));  % remove diagonal elements

cand = setdiff((1:n),currentTeam); % this is the set of candidates after removing current team members

% prune those unpromising candidates
if prune == true
    cand = cand(sum(aa(cand,remainTeam),2)>0);
    %summ = sum(aa(cand,remainTeam),2);
    %index = sum(aa(cand,remainTeam),2)>0; 
    %save('aa_sum_filtered.mat', 'summ', 'cand', 'index');
end

% decay factor in RWR, set it so that it converges
c=0.00000001;

% number of skills
dn=length(L);

score = zeros(length(cand), 2);
edge_weight_matrix = zeros(length(cand), length(currentTeam));

for i=1:length(cand)
    newTeam = [remainTeam,cand(i)];
    LL = zeros(length(newTeam)^2);
    for j=1:dn
        LL = LL + kron(L{j}(currentTeam,currentTeam),L{j}(newTeam,newTeam));
    end
    A2 = aa(newTeam,newTeam);
    A2 = (triu(A2,1) + tril(A2,-1));
    
    score(i,1) = label_gs(A1,A2,LL,c);
    score(i,2) = cand(i);
    edge_weight_matrix(i, length(A2)) = cand(i);
    result = inf_cal(A1,A2,LL,c);
    edge_weight_matrix(i, 1:(length(A2)-1)) = result; 
end

end

function sim = label_gs(A,B,L,c)
% graph kernel computation for attributes graphs
% input:
%   A: first graph
%   B: second graph
%   L: label matrix
%   c: decay factor
% output:
%   sim: similarity between A and B

    n1 = size(A,1);
    n2 = size(B,1);
    q = {ones(n1,1)/n1,ones(n2,1)/n2};
    p = {ones(n1,1)/n1,ones(n2,1)/n2};
    p1 = p{1};
    p2 = p{2};
    q1 = q{1};
    q2 = q{2};
    
    X = kron(A,B);
    qx = kron(q1,q2);
    px = kron(p1,p2);
    sim = qx' * inv(eye(n1 * n2) - c * L * X) * L * px;
end

%graph kernel influcial score computation
function sim_der = inf_cal(A, B, L, c)
    n1 = size(A, 1);
    n2 = size(B,1);
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
        tmp = zeros(n2, n2);
        tmp(n2, i) = 1;
        tmp(i, n2) = 1;
        core = kron(A, tmp);
        sim_der(1, i) = first_part * core * second_part;
    end
    %disp("sim_der size " + size(sim_der));
end