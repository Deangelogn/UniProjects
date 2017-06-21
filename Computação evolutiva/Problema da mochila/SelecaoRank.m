function [ p_out, fnp ] = SelecaoRank(fp, n, p);
%UNTITLED Summary of this function goes here
% fitness da população
    t = length(fp);
    [rp,I] = sort(fp,'ascend');
    j=t;
    i=1;
    lrp = logical(zeros(1,t));
    
    while (i<=n)
        c = rand(1);
        if(c<(j/t) && lrp(j)==0)
            Idx(i)=I(j);
            i=i+1;
            lrp(j)=1;
        end
        j=j-1;
        
        if(j==0)
            j = t;
        end
    end
    p_out = p(1:n,:);
    for k =1:n
        p_out(k,:) = p(Idx(k),:);
        fnp(k) = fp(Idx(k));
    end
    
end

