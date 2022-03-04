clear;
clf;

x1=[0  -4  -4   0  4     1     4     1        4   1   4     1     4   ];
y1=[2   2  -2  -2  -1.6 -1.2  -0.8   -0.4     0  0.4  0.8   1.2   1.6 ];
     


r=1.8;
theta2=2*pi:-pi/32: pi/32;
x2= -2 +r*cos(theta2);
y2= 0 +r*sin(theta2);

pgon = polyshape({x1,x2},{y1,y2});
figure (1)
plot(pgon);
axis equal off

T = triangulation(pgon);

figure (2)
triplot(T);
axis equal off

ps = T.Points';
ts = T.ConnectivityList';
fileID = fopen('Geometry.txt','w');
fprintf(fileID,'Points\r\n');
np=size(ps,2);
fprintf(fileID,'%i \r\n',np);
fprintf(fileID,'%i \t %e \t %e  \r\n', [1:np; ps]);

fprintf(fileID,'Triangles\r\n');
nt=size(ts,2);
fprintf(fileID,'%i \r\n',nt);
fprintf(fileID,'%i \t %i \t %i \t %i \t \r\n', [1:nt ; (ts-1)]);
fclose(fileID);



