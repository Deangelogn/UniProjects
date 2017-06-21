function [ out ] = crossoverOX( k1, k2, c)
%crossoverOX Summary of this function goes here
%   Detailed explanation goes here
out = zeros(size(c));

out(1,:)= [c(1,1:k1) c(2,k1+1:k2) c(1,k2+1:end)];
out(2,:)= [c(2,1:k1) c(1,k1+1:k2) c(2,k2+1:end)];

end

