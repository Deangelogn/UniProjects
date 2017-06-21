function [ cadeia_out ] = MutacaoInvBit( k,cadeia_in )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
if(cadeia_in(k)==1)
    cadeia_in(k)=0;
else
    cadeia_in(k)=1;
end
    cadeia_out = cadeia_in;
end

