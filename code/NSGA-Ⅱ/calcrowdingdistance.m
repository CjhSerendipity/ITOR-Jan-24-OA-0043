function pop = calcrowdingdistance(pop,F)
    nf = numel(F);
    for i = 1 : nf
        objs = [(pop(F{i}).obj)];
        [nobj, n] = size(objs);
        d = zeros(n,nobj);
        for j = 1 : nobj
            [cj, so] = sort(objs(j,:));
            d(so(1), j) = inf;
            for k = 2 : n - 1
                d(so(k), j) = abs(cj(k+1) - cj(k-1))/abs(cj(1) - cj(end));
            end
            d(so(end), j) = inf;
        end
        for m = 1 : n
            pop(F{i}(m)).crowdingdistance = sum(d(m,:));
        end
    end

end