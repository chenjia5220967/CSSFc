function   [cls_idx,vec,cls_num]  =  My_kmeans(Y, cls_num, itn)
Y         =   Y';
[L b2]    =   size(Y);
P         =   randperm(L);
P2        =   P(1:cls_num);
vec       =   Y(P2(1:end), :);
m_num     =   2000;

for i = 1 : itn
    cnt       =  zeros(1, cls_num);    
    
    v_dis    =   zeros(L, cls_num);
    for  k = 1 : cls_num
        v_dis(:, k) = (Y(:,1) - vec(k,1)).^2;
        for c = 2:b2
            v_dis(:,k) =  v_dis(:,k) + (Y(:,c) - vec(k,c)).^2;
        end
    end      

    [val cls_idx]     =   min(v_dis, [], 2);
    
    [s_idx, seg]   =  Proc_cls_idx( cls_idx );
    for  k  =  1 : length(seg)-1
        idx    =   s_idx(seg(k)+1:seg(k+1));    
        cls    =   cls_idx(idx(1));    
        vec(cls,:)    =   mean(Y(idx, :));
        cnt(cls)      =   length(idx);
    end        
    
    if (i==itn-2)
        [val ind]  =  min( cnt );       % Remove these classes with little samples        
        while (val<m_num) && (cls_num>=40)
            vec(ind, :)    =  [];
            cls_num       =  cls_num - 1;
            cnt(ind)      =  [];
            [val  ind]    =  min(cnt);
        end        
    end
end