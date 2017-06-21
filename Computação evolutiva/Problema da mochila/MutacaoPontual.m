function [ cadeia_out ] = MutacaoPontual( k,cadeia_in )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
l = k(2)-k(1);
pm = randi([0,1],1,l);
for i=1:l
    if(pm(i)==1)
        if(cadeia_in(k(1)+i)==1)
            cadeia_in(k(1)+i)=0;
        else
            cadeia_in(k(1)+i)=1;
        end
    end
end
cadeia_out = cadeia_in;
end

