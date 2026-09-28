function EP = Forward(EP, Offspring, N)
    %% Dominance relation calculation
    S = [EP, Offspring];
    [FrontNo,~] = NDSort(S.objs,inf);
    EP = S(FrontNo==1);
    PopObj = EP.objs;
    [P, ~] = size(PopObj);
    W = UniformPoint(N,size(PopObj,2));
    if P <= N
        i=2;
        while length(EP) < N
            EP = [EP,S(FrontNo==i)];
            i = i+1;
        end
    end

    while length(EP) > N
        PopObj = EP.objs;
        x_wrost=SelectWrost(EP,W,PopObj);
        EP=setdiff(EP,x_wrost);
    end
end

function x_wrost = SelectWrost(EP,W,PopObj)
    zmax = max(PopObj, [], 1);
    zmin = min(PopObj,[],1);
    PopObj = (EP.objs-zmin) ./ (zmax - zmin);
    PopObj = NearZero(PopObj);
    W = NearZero(W);
    [~,Region] = max(1-pdist2(PopObj,W,'cosine'),[],2);
    [value,~]=sort(Region,'ascend');
    flag=max(value);
    counter=histc(value,1:flag);
    [~,most_crowded]=max(counter);
    idx_crowdest=EP(Region==most_crowded);

    dist=pdist2(idx_crowdest.objs,idx_crowdest.objs);
    dist(dist==0)=inf;
    [row,~]=find(min(min(dist))==dist);
    St=idx_crowdest(row);
    [~,Region_FEA] = max(1-pdist2(St.objs,W,'cosine'),[],2);
    Z = min(St.objs,[],1);
    g_tch=max(abs(St.objs-repmat(Z,length(St),1))./W(Region_FEA,:),[],2);
    [~,order]=max(g_tch);
    x_wrost=St(order);
end
function mat = NearZero(mat)
    threshold = 1e-11;
    norm = sqrt(sum(mat.^2,2));
    idx_mat_small = norm < threshold;
    mat(idx_mat_small) = 0;
end
