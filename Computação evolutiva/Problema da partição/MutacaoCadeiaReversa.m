function [ cadeia_out ] = MutacaoCadeiaReversa( k,cadeia_in )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
    cadeia_out = cadeia_in(1:k(1));
    for i=2:length(k)
        cadeia_out = [cadeia_out flip(cadeia_in(k(i-1)+1:k(i)))];
    end
        cadeia_out = [cadeia_out cadeia_in(k(i)+1:end)];
 end

