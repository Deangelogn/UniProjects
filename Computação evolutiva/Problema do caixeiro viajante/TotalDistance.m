function [ dist ] = TotalDistance( v )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
v(end+1,:) = v(1,:);
dist = sum(sqrt(sum((v(1:end-1,:) -v(2:end,:)).^2,2)));
end

