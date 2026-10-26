%% Read data set function
function [ngmar,num]=readfile(i)
if i==1
    E=load('.\Data\breastmnist_224.mat');
    ngmar=cat(1,E.train_images,E.test_images,E.val_images);
    id=sort(randperm(size(ngmar,1),100));
    ngmar=double(ngmar(id,:,:));
    ngmar=reshape(ngmar,100,8,28,8,28);
    ngmar=mean(mean(ngmar,2),4);
    ngmar=reshape(ngmar,100,28,28);
    ngmar=permute(ngmar,[2,3,1]);
elseif i==2
    E=load('.\Data\PNW-ML.mat');
    ngmar=double(E.fea);
    ngmar=reshape(ngmar,100,size(ngmar,1)/100,3,size(ngmar,3));
    id=round(linspace(1,size(ngmar,2),20));
    ngmar=ngmar(:,id,:,:);
end
ngmar=ngmar*(10/norm(ngmar(:)));
num=length(size(ngmar));
end
