function newChromo = mutate(chromo,mu)
    if rand<=mu
        newChromo = initial_population(1);
    else
        newChromo = chromo;
    end
end
