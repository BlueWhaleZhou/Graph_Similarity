currentTeam = [750967, 124372, 167178, 620387, 3603, 582475];
newTeam = [750967, 167178, 620387, 3603, 582475, 247326];
g_1 = inf_old_cell{1} * 10000;
disp(g_1);
g_2 = inf_new_cell{1} * 10000;
disp(g_2);
len = length(g_1);

s_1 = [];
t_1 = [];
weights_1 = [];

s_2 = [];
t_2 = [];
weights_2 = [];

for i = 1:len
    for j = i+1:len
        s_1 = [s_1, i];
        t_1 = [t_1, j];
        weights_1 = [weights_1, g_1(i, j)];
     end
end

for i = 1:len
    for j = i+1:len
        s_2 = [s_2, i];
        t_2 = [t_2, j];
        weights_2 = [weights_2, g_2(i, j)];
     end
end

name_old = strings([1, len]);
name_new = strings([1, len]);

for i = 1:length(currentTeam)
    name_old(i) = string(authorDict{currentTeam(i)});
end

for i = 1:length(newTeam)
    name_new(i) = string(authorDict{newTeam(i)});
end

name_old = cellstr(name_old);
name_new = cellstr(name_new);
G_1 = graph(s_1, t_1, weights_1, name_old);
G_2 = graph(s_2, t_2, weights_2, name_new);
LWidths_1 = 10 * G_1.Edges.Weight / max(G_1.Edges.Weight);
LWidths_2 = 10 * G_2.Edges.Weight / max(G_2.Edges.Weight);

%plot(G_1, 'EdgeLabel', G_1.Edges.Weight, 'LineWidth', LWidths_1);
plot(G_2, 'EdgeLabel', G_2.Edges.Weight, 'LineWidth', LWidths_2);