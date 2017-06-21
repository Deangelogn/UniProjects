clear all; clc;

  [X,Y] = meshgrid(-1:.02:2, -1:.02:2);
        Z =sin(4*pi*X).*X-sin(4*pi*Y+pi).*Y+1;
        
        % sin(4*pi*xx).*xx-sin(4*pi*yy+pi).*yy+1
        figure
        surf(X,Y,Z)
        xlabel('x')
        ylabel('y')
        zlabel('f(x,y)')
        %legend('Fitness function')
        
        figure
       contour(X,Y,Z)
       xlabel('x')
        ylabel('y')
 