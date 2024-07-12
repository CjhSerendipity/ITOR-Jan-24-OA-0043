function [hv] = HV(F1obj)
% Efficient method for 3D objective function values
%
% Implementation after:
%
%   'B. Naujoks, N. Beume, M. Emmerich. Multi-objective
%    optimisation using S-metric Selection: Application to
%    three-dimensional Solution Spaces. CEC 2005, 1282-1289,
%    2005.'
    F = F1obj';
    F(1,:) = F(1,:)./-1;
    F(2,:) = F(2,:)./-1;
    F(3,:) = F(3,:)./1;
%     f2 = F(2,:);
%     f3 = F(3,:);
%     F(3,:) = f2;
%     F(2,:) = f3;
    % ub = [max(F1obj(:,1));max(F1obj(:,2));max(F1obj(:,3))];
    ub = [0;0;5];
	[M, l] = size(F);

	a = sort(F(1,:), 'ascend');
	b = sort(F(2,:), 'ascend');
	a(end+1) = ub(1);
	b(end+1) = ub(2);

	best1_f3 = ub(3) * ones(l,l);
	best2_f3 = ub(3) * ones(l,l);

	for t = 1:l
		for i = 1:l
			for j = 1:l
				if (F(1,t) <= a(i) && F(2,t) <= b(j))
					if (F(3,t) < best1_f3(i,j))
						best2_f3(i,j) = best1_f3(i,j);
						best1_f3(i,j) = F(3,t);
					elseif (F(3,t) < best2_f3(i,j))
						best2_f3(i,j) = F(3,t);
					end
				end
			end
		end
	end

	hypervolume_contributions = zeros(1, l);
	for t = 1:l
		for i = 1:l
			for j = 1:l
				if (F(1,t) <= a(i) && F(2,t) <= b(j) && F(3,t) == best1_f3(i,j))
					hypervolume_contributions(t) = hypervolume_contributions(t) + (a(i+1) - a(i)) * (b(j+1) - b(j)) * (best2_f3(i,j) - best1_f3(i,j));
				end
			end
        end
    end
    hv = sum(hypervolume_contributions);
end


