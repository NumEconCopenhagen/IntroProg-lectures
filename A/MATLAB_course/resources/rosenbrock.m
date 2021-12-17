function [ z ] = rosenbrock( x,y )
%ROSENBROCK 

a = 1; 
b = 1; 
z = (a - x)^2 + b*(y-x^2)^2; 

end

