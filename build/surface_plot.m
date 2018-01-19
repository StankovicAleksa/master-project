[X,Y] = meshgrid(0:0.1:10,0:1.0/29:1);
Z = sin(X) + cos(Y);
filename = 'test.txt';
M = csvread(filename);
Z = M(:,2:31);
size(Z)
size(X)
size(Y)
surf(X,Y,Z')
pbaspect([10 1 1])