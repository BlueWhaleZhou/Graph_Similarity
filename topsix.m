function ids = topsix(score)
    score = sortrows(score,-1);
    ids = score(1:6,2);
end