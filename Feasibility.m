function MP = Feasibility(MP, Offspring, N)
    H = [MP, Offspring];
    NonDominated = Domination(H);
    MP = H(NonDominated);
    PopObj = MP.objs;
    [P, ~] = size(PopObj);
    if P <= N
        while length(MP) < N
            H = setdiff(H,MP);
            NonDominated = Domination(H);
            MP = [MP,H(NonDominated)];
        end
    end
    if length(MP) > N
        PopCon = MP.cons;
        Cons   = sum(max(0,PopCon),2);
        if sum(Cons <= 0) <= N
            [~, index] = sort(Cons);
            remained_solution = index(1: N);
            MP = MP(remained_solution);
        else
            MP    = MP(Cons <= 0);
            PopObj = MP.objs;
            [P, ~] = size(PopObj);
            zmin   = min(PopObj, [], 1);
            zmax   = max(PopObj, [], 1);
            PopObj = (PopObj - zmin) ./ (zmax - zmin);
            Del    = Truncation(PopObj, P - N);
            remained_solution = 1 : P;
            remained_solution(Del) = [];
            MP = MP(remained_solution);
        end
    end
end

function Del = Truncation(PopObj,K)
% Select part of the solutions by truncation
    Distance = pdist2(PopObj,PopObj);
    Distance(logical(eye(length(Distance)))) = inf;
    Del = false(1,size(PopObj,1));
    while sum(Del) < K
        Remain   = find(~Del);
        Temp     = sort(Distance(Remain,Remain),2);
        [~,Rank] = sortrows(Temp);
        Del(Remain(Rank(1))) = true;
    end
end
