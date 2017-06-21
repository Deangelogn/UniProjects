function [ np, fnp ] = TorneioPE( q,pop,fp, pop_size)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

t = randi(length(fp),length(fp),q);
idx = 1:length(fp);
idx = meshgrid(idx)';
idx = idx(1:length(idx),1:q);
find(idx==t);

while(length(find(idx==t))>0)
    x = length(find(idx==t));
    t(find(idx==t)) = randi(length(fp),x,1);
end

opp = fp(t);
vMatch = meshgrid(fp)';
vMatch = vMatch(1:length(fp),1:q);
numWin = sum(vMatch<opp,2);

[rp,I] = sort(numWin,'descend');

for i=1:pop_size
    fnp(i) = fp(I(i));
    np(i,:) = pop(I(i),:);
end

end

