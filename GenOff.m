function Offspring = GenOff(Problem,UP,DP,MP)
    Offspring1 = [];
    Pop1 = UP;
    Pop2 = DP;
    Pop3 = [UP, MP,DP];

    MatingPool_Pop1_1 = randperm(length(Pop1));
    MatingPool_Pop1_2 = randperm(length(Pop1));
    MatingPool_Pop2_1 = randperm(length(Pop2));
    MatingPool_Pop2_2 = randperm(length(Pop2));
    MatingPool_Pop3_1 = randperm(length(Pop3));
    MatingPool_Pop3_2 = randperm(length(Pop3));

    if rand() < 0.5
        Offspring1 = OperatorDE(Problem, Pop1, Pop1(MatingPool_Pop1_1), Pop1(MatingPool_Pop1_2),{1,0.5,1,1});
        Offspring2 = OperatorDE(Problem, Pop2, Pop2(MatingPool_Pop2_1), Pop2(MatingPool_Pop2_2),{1,0.5,1,1});
        Offspring3 = OperatorDE(Problem, MP, Pop3(MatingPool_Pop3_1), Pop3(MatingPool_Pop3_2),{1,0.5,1,1});
    else
        Offspring1 = OperatorGA(Problem, Pop1(MatingPool_Pop1_1), {1,20,1,1});
        Offspring2 = OperatorGA(Problem, Pop2(MatingPool_Pop2_1), {1,20,1,1});
        Offspring3 = OperatorGA(Problem, Pop3(MatingPool_Pop3_1), {1,20,1,1});
    end
    Offspring = [Offspring1, Offspring2, Offspring3];
end
