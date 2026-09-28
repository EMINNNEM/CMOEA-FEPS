function DP = Diversity(DP, Offspring, zmin, Ns)
    %% Dominance relation calculation
    H  = [DP, Offspring];
    DP = [];
    c_NonDominated = Domination(H);
    NonDominated = NDSort(H,1);
    S1   = H(c_NonDominated);
    S2 = H(NonDominated);
    S = [S1,S2];
    S = unique(S);
    H = setdiff(H,S);
    Obj = S.objs;
    Con = S.cons;
    CV  = sum(max(0,Con),2);
    count = 0;
    %% Minimum angular distance calculation
    [W,Ns] = UniformPoint(Ns, size(Obj,2));
    Angle_W_to_W = acos(1 - pdist2(W,W,'cosine'));
    Angle_W_to_W(eye(Ns) == 1) = inf;
    h = mean(min(Angle_W_to_W));
    zmax = max(Obj,[],1);
    Obj          = (Obj - zmin) ./ (zmax - zmin);
    Angle_S_to_W = sin(acos(1 - pdist2(W,Obj,'cosine')));

    for i = 1 : Ns
        Angle = Angle_S_to_W(i,:);
        list  = Angle <= h;
        if sum(list) == 0
            [~,index]   = min(Angle_S_to_W(i,:));
            list(index) = true;
        end
        T        = inf(length(S), 1);
        T(list)  = CV(list);
        feasible = T <= 0;
        if sum(feasible) <= 0
            [~,index] = min(T);
        else
            T = inf(size(Angle));
            Fitness     = Angle_S_to_W(i,:);
            T(feasible) = Fitness(feasible);
            [~,index]   = min(T);
        end
        DP = [DP,S(index)];
        Angle_S_to_W(:,index) = inf;
        if (i - count)==size(Angle_S_to_W,2)
            c_NonDominated = Domination(H);
            S   = H(c_NonDominated);
            H = setdiff(H,S);
            Obj = S.objs;
            Con = S.cons;
            CV  = sum(max(0,Con),2);

            if sum(CV <= 0) > 0
                zmax = max(Obj(CV <= 0, :), [], 1);
            else
                [~,index] = min(CV);
                zmax = Obj(index,:);
            end
            Obj          = (Obj - zmin) ./ (zmax - zmin);
            Angle_S_to_W = sin(acos(1 - pdist2(W,Obj,'cosine')));
        end
        count = count+1;
    end
end
