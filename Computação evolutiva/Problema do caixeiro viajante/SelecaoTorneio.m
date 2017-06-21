function [np, fnp] = SelecaoTorneio( q,p,fp,n )
%UNTITLED Summary of this function goes here
%fp = fitnessda população
%n = número de indivíduos da nova geração
%q = número de indivíduos para o torneio
%p = população

    np = zeros(n,size(p,2)); % Nova geração de indivíduos
    fnp = zeros(n,1);
    t = randi(size(p,1),n,q); % t = torneio com q indivíduos
    [fnp, minIdx]=min(fp(t),[],2);
   
    for j=1:n; 
       np(j,:) = p(t(j,minIdx(j)),:);
       w(j)=(t(j,minIdx(j)));  
    end
  
end

