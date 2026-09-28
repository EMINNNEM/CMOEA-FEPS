classdef CMOEA < ALGORITHM
% <2025> <multi/many> <real/binary/permutation><constrained/none>
% [LOG-ONLY PATCH 2026-09-22] A.metric.traj 由 7 列扩展为 14 列，行为零改动：
%   1..7  = [FE, gen, size, FR, VR2, CR, IGDp]  （与旧格式逐位一致）
%   8..14 = [VR, c, tempRaw, tempClamped, sizePrev, logTerm, stallFlag]
%   tempRaw/logTerm/stallFlag 为夹取前目标规模、对数项、停滞分支标志；
%   预热期（gen<=last_gen）无更新，第 8..14 列为 NaN。
%   备份：/tmp/opencode/CMOEA.m.bak-20260922


    methods
        function main(Algorithm,Problem)
            %% Parameter setting
            Ns = floor(Problem.N/3);
            %% Generate random population
            Population = Problem.Initialization();
            %设置FR可行解比例阈值a1和change_rate变化度比率a2
            a1 = 0.95;
            a2 = 0.05;
            EP  = Population(1:ceil(Problem.N/2));
            DP  = [];
            MP = Population;
            Offspring = Population;
            zmin      = min(Population.objs,[],1);
            ideal_points = zeros(ceil(Problem.maxFE/Problem.N),Problem.M);
            nadir_points = zeros(ceil(Problem.maxFE/Problem.N),Problem.M);

            ideal_points2 = zeros(ceil(Problem.maxFE/Problem.N),Problem.M);
            nadir_points2 = zeros(ceil(Problem.maxFE/Problem.N),Problem.M);

            VR = 1;
            VR2 = 1;

            gen = 1;
            last_gen = 20;
            FR = sum(sum(max(MP.cons,0),2)==0)/length(MP); %FEA的可行解比例
            CR = 0;
            a = 0.1;
            b = 0.2;
            c = 1;
            size = Problem.N;
            vec_size = [];
            score =0;

            %% Optimization
            while Algorithm.NotTerminated2(MP,DP)

                ideal_points(gen,:) = min(MP.objs,[],1);
                nadir_points(gen,:) = max(MP.objs,[],1);

                ideal_points2(gen,:) = min(EP.objs,[],1);
                nadir_points2(gen,:) = max(EP.objs,[],1);

                zmin = min(zmin, min(Offspring.objs, [], 1) - 1e-6);

                FR = sum(sum(max(MP.cons,0),2)==0)/length(MP);  %更新FR
                % [LOG-ONLY 2026-09-22] 观测变量（不改变任何算法行为）
                sizePrev = size; logTerm = NaN; tempRaw = NaN; temp = NaN; stallFlag = NaN;
                if gen > last_gen
                    VR = calc_maxchange(ideal_points,nadir_points,gen,last_gen);
                    VR2 = calc_maxchange(ideal_points2,nadir_points2,gen,last_gen);

                    disp(ideal_points(gen,:));
                    disp(ideal_points2(gen,:));
                    %更新辅助种群大小
                    logTerm = log2((1+CR)*(1+VR2)/(1+FR/2));  % [LOG-ONLY]
                    tempRaw = ceil(size*logTerm*c);            % [LOG-ONLY] 夹取前目标规模
                    temp = tempRaw;

                    disp((1+CR)*(1+VR2)/(1+FR));
                    if temp < 5
                        temp = 5;
                    elseif temp > Problem.N
                        temp = Problem.N;
                    end
                    size = ceil(size*0.7+temp*0.3);
                    disp(size);
                    if VR2 < a2 && VR < a2
                        c = 1 - Problem.FE/Problem.maxFE;  % [LOG-ONLY] stallFlag
                        stallFlag = 1;
                    else
                        c = c+0.1;                         % [LOG-ONLY] stallFlag
                        stallFlag = 0;
                    end
                end


                MP  = Feasibility(MP, Offspring, Problem.N);
                DP   = Diversity(DP, Offspring, zmin, Ns);
                UPop_Off = Offspring(1:size);
                CR = sum(ismember(MP.objs,UPop_Off.objs,"rows")) / length(EP);
                EP   = Forward(EP, Offspring ,size);

                Offspring = GenOff(Problem,EP,DP,MP);
                gen = gen+1;


                score = IGDp(MP,Problem.optimum);
                % [LOG-ONLY 2026-09-22] traj 追加观测列 8..14：
                %   [VR, c, tempRaw, tempClamped, sizePrev, logTerm, stallFlag]
                %   前 7 列 [FE, gen, size, FR, VR2, CR, IGDp] 与旧格式逐位一致。
                vec_size = cat(1,vec_size,[Problem.FE,gen,size,FR,VR2,CR,score, ...
                    VR, c, tempRaw, temp, sizePrev, logTerm, stallFlag]);
                if Problem.FE >= Problem.maxFE
                    Algorithm.metric.traj = vec_size;
                end
                % if mod(gen,60)==0
                %     filename = sprintf('Result_gen_%04d.csv', gen);
                %     writematrix(UP.objs,filename);
                % end
                % if Problem.FE >= Problem.maxFE
                %     writematrix(vec_size,'C:\Code\res\size.csv');
                % end

            end
        end
    end
end

function VR = calc_maxchange(ideal_points,nadir_points,gen,last_gen)
    delta_value = 1e-6 * ones(1,size(ideal_points,2));
    rz = abs((ideal_points(gen,:) - ideal_points(gen - last_gen + 1,:)) ./ max(ideal_points(gen-last_gen+1,:),delta_value));
    nrz = abs((nadir_points(gen,:) - nadir_points(gen - last_gen + 1,:)) ./ max(nadir_points(gen-last_gen + 1,:),delta_value));
    VR = max([rz,nrz]);
end
