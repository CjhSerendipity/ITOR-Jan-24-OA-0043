for a= 1:30
    for b = 1:7
        v2(a,b) = sum(V_ltBp(:,b,a))+sum(V_ltHp(:,b,a));
    end
end
for a= 1:7
    for b = 1:7
        y2(a,b) = sum(y_tpp(a,b,:));
    end
end
for a= 1:7
    for b = 1:7
        z2(a,b) = sum(z_tpp(a,b,:));
    end
end