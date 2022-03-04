clear;
clf;
fid = fopen('InitialMesh.txt');
string1=fgets(fid);
ns = fscanf(fid,'%i',1);
M= fscanf(fid,'%*i %f %f %f \n',[3 ns]);
ps=M(1:2,:); rs=M(3,:);

string2=fgets(fid);
nts = fscanf(fid,'%i',1);
ts = fscanf(fid,'%*i %i %i %i \n',[3 nts])+1;

fclose(fid);

figure (1)
pdemesh(ps,ts);
%pdemesh(ps,ts,'NodeLabels','on');
%pdemesh(nodes,elements,'ElementLabels','on');
% xlabel('$x$','interpreter','latex');
% ylabel('$y$','interpreter','latex');
% xlim([-1 1]);
% ylim([-1 1]);
% xticks([-1 0 1]);yticks([-1 0 1]);

hold on  
for i=1:length(rs)
rectangle('Position',[ps(1,i)-rs(i),ps(2,i)-rs(i),2*rs(i),2*rs(i)], ...
    'Curvature', [1 1], 'EdgeColor','r');
end
%axis equal
axis off



