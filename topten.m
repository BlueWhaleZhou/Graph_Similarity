function ids = topten(score)
    score = sortrows(score,-1);
    ids = score(1:10,2);
end