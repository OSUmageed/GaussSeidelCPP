clear
close all

x = dlmread('GS_output.txt');
Lx = x(1);
Ly = x(2);
ds = x(3);
Temperature = x(4:end)';
x1 = 0:ds:Lx;
y1 = 0:ds:Ly;
[X,Y] = meshgrid(x1,y1);
Temper = reshape(Temperature,size(X));
surf(X,Y,Temper,'edgecolor','none')
colorbar
zlabel('Temperature')
