clear;
clf;

x1=[-4   4   4  -4 ];
y1=[-2  -2   2   2];
     
r1=1.4;
theta1=2*pi:-pi/32: pi/32;
x2=-2 +r1*cos(theta1);
y2= 0 +r1*sin(theta1);

r2=1.8;
theta2=2*pi:-pi/32: pi/32;
x3= 2 +r2*cos(theta2);
y3= 0 +r2*sin(theta2);

pgon = polyshape({x1,x2,x3},{y1,y2,y3});
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



