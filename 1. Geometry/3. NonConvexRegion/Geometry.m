clear;
clf;


x1=[0.3 0  0   1   1  0.6 0.6 0.2 0.5 0.5 ];
y1=[1   1  0   0   1   1  0.2 0.2 0.5 0.8 ];


N2=15;
theta2=pi/2*(1:(N2-1))/N2;
x1=[x1,(0.3+0.2*cos(theta2))];
y1=[y1,(0.8 +0.2*sin(theta2))];


pgon = polyshape(x1,y1);
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



