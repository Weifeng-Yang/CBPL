%% Read data set function
function [ngmar,dimension]=readfile(i)
if i==1
    E=load('.\Data\lp_ship12l.mat');
    E=E.Problem;
    ngmar=double(E.A);
elseif i==2
    E=load('.\Data\BASEHOCK.mat');
    ngmar=double(E.X);
end
ngmar=ngmar*(10/norm(ngmar(:)));
dimension=size(ngmar);
end
