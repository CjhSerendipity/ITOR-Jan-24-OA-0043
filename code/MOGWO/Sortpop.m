function pop = Sortpop(pop)
    %拥挤度降序排列
    [~, CDSO] = sort([pop.crowdingdistance], 'descend'); 
    pop = pop(CDSO);
    %非支配等级升序排列
    [~, RSO] = sort([pop.rank]);
    pop = pop(RSO);
end
