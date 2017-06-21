function [ out ] = SelecaoRoleta( fp,n )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
total = sum(fp);
w = fp/total;
w = cumsum(w);
lrp = logical(zeros(1,length(fp)));
    
    for i=1:n
        c = rand(1);
        Idx(i) = find(w>c,1,'first');
    end
    out = Idx;
end

