function z = dominate(p, q)

    c1 = [p.obj];
    c2 = [q.obj];
    
    z = all(c1 <= c2) & any(c1 < c2);

end