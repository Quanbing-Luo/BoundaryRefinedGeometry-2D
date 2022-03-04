clear;
clf;

x1=[-3    1   1 ];
y1=[0   -1   1  ];

N2=50;r=0.65;
theta2=-2*pi*(0:(N2-1))/N2;
x2=r*cos(theta2);
y2=r*sin(theta2);

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



