function [np, fnp] = SelecaoTorneio( q,p,fp,n )
%UNTITLED Summary of this function goes here
%fp = fitnessda popula��o
%n = n�mero de indiv�duos da nova gera��o
%q = n�mero de indiv�duos para o torneio
%p = popula��o

    np = logical(zeros(n,size(p,2))); % Nova gera��o de indiv�duos
    fnp = zeros(n,1);
    t = randi(size(p,1),n,q); % t = torneio com q indiv�duos
    [fnp, maxIdx]=max(fp(t),[],2);
   
    for j=1:n; 
       np(j,:) = p(t(j,maxIdx(j)),:);
       w(j)=(t(j,maxIdx(j)));  
    end
%      pause;
%      
end

