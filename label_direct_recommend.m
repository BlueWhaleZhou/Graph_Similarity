function [score, edge_weight_matrix, result_max] = label_direct_recommend(aa,L,currentTeam,i0,prune)
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

%get author name with author id
fileID=fopen('authorDict.txt');
authorDict=textscan(fileID,'%s','delimiter','\n');
authorDict=authorDict{1};

%expand current team for a incomplete graph
l1 = length(currentTeam);
disp(currentTeam);

%old team adjacency matrix
old_team_new = zeros([l1, l1]);
old_team = aa(currentTeam, currentTeam);

old_team =(triu(old_team, 1) + tril(old_team, -1));

for i = 1:l1
    for j = 1:l1
        old_team_new(i, j) = old_team(i, j);
    end
end

s = [];
t = [];
weights = [];
for i = 1:l1
    for j = i+1:l1
        s = [s, i];
        t = [t, j];
        weights = [weights, old_team_new(i, j)];
    end
end

%disp(s);
%disp(t);
%disp(weights);
%disp(old_team_new);
%plot old team
name_old = strings([1, l1]);

%node name
for i = 1:length(currentTeam)
    name_old(i) = string(authorDict{currentTeam(i)});
end
name_old = cellstr(name_old);
G = graph(s, t, weights, name_old);
LWidths = 5 * G.Edges.Weight / max(G.Edges.Weight);
%plot(G, 'EdgeLabel', G.Edges.Weight, 'LineWidth', LWidths);

n=size(aa,1);
remainTeam = setdiff(currentTeam,i0,'stable'); %i0 is the index of author to be replaced
currentTeam = [remainTeam, i0];

A1 = aa(currentTeam,currentTeam);
A1 = (triu(A1,1) + tril(A1,-1));  % remove diagonal elements


cand = setdiff((1:n),currentTeam,'stable'); % this is the set of candidates after removing current team members

% prune those unpromising candidates
if prune == true
    cand = cand(sum(aa(cand,remainTeam),2)>0);
end

% decay factor in RWR, set it so that it converges
c=0.00000001;

% number of skills
dn=length(L);

score = zeros(length(cand), 2);
result_max = zeros(length(cand), 2);
edge_weight_matrix = zeros(length(cand), length(currentTeam));

for i=1:length(cand)
    newTeam = [remainTeam, cand(i)];
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
    result_max(i, 1) = max(result);
    result_max(i, 2) = cand(i);
    edge_weight_matrix(i, 1:(length(A2)-1)) = result; 
end

end
