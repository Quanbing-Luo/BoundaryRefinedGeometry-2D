clear;
clf;

x1=[-0.5 -0.5 0.5 0.5];
y1=[0.5 -0.5 -0.5 0.5];

r=0.45;
theta=pi:-pi/12:-pi/2;
x2=[r*cos(theta),-r];
y2=[r*sin(theta),-r];

pgon = polyshape({x1,x2},{y1,y2});


figure (1)
plot(pgon);
axis equal off

tr = triangulation(pgon);

% figure (2)
% triplot(tr);
% axis equal off

model = createpde;
tnodes = tr.Points';
telements = tr.ConnectivityList';

geometryFromMesh(model,tnodes,telements);

figure (2)
pdegplot(model)
axis equal off

figure (3)
pdemesh(model);
axis equal off

generateMesh(model);
figure (4)
pdemesh(model)
axis equal off


ps = tr.Points';
ts = tr.ConnectivityList';
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



