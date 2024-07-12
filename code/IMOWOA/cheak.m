function [x_ktt] = cheak(x_kt)
global T X_1t X_2t
    for t = 1:1:T
        x_k = x_kt(t,:);
        s_x = sum(x_k(1:4));
        s_x2 = sum(x_k(5:6));
        for h = 1:4
            x_k(h) = fix((x_k(h)/s_x)*X_1t(t));
        end
        for h = 1:2
            if fix(s_x2) >= X_2t(t)
                x_k(4+h) = fix((x_k(4+h)/s_x2)*X_2t(t));
            else
                x_k(4+h) = fix(x_k(4+h));
            end  
        end
        s_x3 = sum(x_k(7:8));
        if fix(x_k(7)+x_k(8)) > x_k(2)
            for h = 1:2
                x_k(6+h) = fix((x_k(6+h)/s_x3)*x_k(2));
            end
        else
            x_k(7) = fix(x_k(7));
            x_k(8) = fix(x_k(8));
        end
        if fix(x_k(9)) > x_k(1)+x_k(8)+x_k(3)
            x_k(9) = x_k(1)+x_k(8)+x_k(3);
        else
            x_k(9) = fix(x_k(9));
        end
        if fix(x_k(10)) > x_k(1)+x_k(2)+x_k(3) + 3*x_k(5)
            x_k(10) = x_k(1)+x_k(2)+x_k(3) + 3*x_k(5);
        end
            x_k(10) = fix(x_k(10));
        if fix(x_k(11)) > x_k(3)
            x_k(11) = x_k(3);
        else
            x_k(11) = fix(x_k(11));
        end
        x_ktt(t,:) = x_k;
    end   
end

