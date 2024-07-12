function z = dominate(p, q)

    c1 = [p.Cost];
    c2 = [q.Cost];
    
    z = all(c1 <= c2) & any(c1 < c2);

end