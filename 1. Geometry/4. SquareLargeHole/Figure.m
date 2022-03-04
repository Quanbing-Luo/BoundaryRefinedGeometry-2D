clear;
clf;
fid = fopen('Geometry.txt');
string1=fgets(fid);
nps = fscanf(fid,'%i',1);
ps= fscanf(fid,'%*i %f %f \n',[2 nps]);

string2=fgets(fid);
nts = fscanf(fid,'%i',1);
ts= fscanf(fid,'%*i %i %i %i\n',[3 nts])+1;

fclose(fid);


figure (1)
pdeplot(ps,ts);
%pdeplot(ps,ts,'NodeLabels','on');
axis equal off




