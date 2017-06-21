function [ s1 s2 ] = CrossoverPosition( f1, f2, vp)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
strSize = length(f1);
s1 = zeros(1,strSize);
s2 = s1;

for i=1:2:length(vp)
    s1(vp(i): vp(i+1)) = f1(vp(i): vp(i+1));
    s2(vp(i): vp(i+1)) = f2(vp(i): vp(i+1));
end

emptyPos = find(s1==0); 

idx = 1;
while(length(emptyPos)>0)
    if( s1(emptyPos(1))==0 && isempty(find(s1==f2(idx))) )
        s1(emptyPos(1)) = f2(idx);
    else
        idx = mod(idx,strSize)+1;
    end
    
    if( s1(emptyPos(1))~=0 )
        emptyPos(1)=[];
    end
end

emptyPos = find(s2==0);
idx = 1;
while(length(emptyPos)>0)
    if( s2(emptyPos(1))==0 && isempty(find(s2==f1(idx))) )
        s2(emptyPos(1)) = f1(idx);
    else
        idx = mod(idx, strSize)+1;
    end
    
    if( s2(emptyPos(1))~=0 )
        emptyPos(1)=[];
    end
end


end

