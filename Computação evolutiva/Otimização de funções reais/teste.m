x =-1:0.01:1;
y =-1:0.01:1;

[X,Y] = meshgrid(x,y)

for i =1:201
    for j =1:201
        f(i,j) = x(i)^2 - y(j)^2;
    end
end

plot3(x,y,f);