clearvars
clc
maxNumCompThreads(1);
cd(fileparts(mfilename('fullpath')));

%% Parameter settings
index=[1,2];
R=[20,50];
seed=[20261026,20261028];
percent=0.3;
lastpercent=1;
alpha=0.01;
epsilon=0.01;
delta=1e-4;
timeend=60;
maxiteropt=1e8;
names={'lp_ship12l','BASEHOCK'};

%% Algorithm iteration starts
for iter=1:length(index)
    rng(seed(index(iter)),'twister');
    [ngmar,dimension]=readfile(index(iter));
    Lapk=graphmat(ngmar);
    num=length(dimension);
    var=cell(1,num);
    aa=zeros(1,num);
    for j=1:num
        var{j}=rand(dimension(j),R(index(iter)));
        aa(j)=numel(var{j});
    end
    for j=1:num-1
        aa(j)=floor(percent*numel(var{j}));
        var{j}=PROX(zeros(size(var{j})),var{j},aa(j));
        var{j}=var{j}./max(sqrt(sum(var{j}.^2,1)),eps);
    end
    aa(num)=floor(lastpercent*numel(var{num}));
    if aa(num)<numel(var{num})
        active=1:min(size(var{num},2),aa(num));
        inactive=setdiff(1:size(var{num},2),active);
        for k=1:dimension(num)
            e=zeros(dimension(num),1); e(k)=1;
            if norm(project(e,Lapk))>0, break; end
        end
        var{num}=zeros(size(var{num}));
        var{num}(:,active)=repmat(e,1,length(active));
        for j=1:num-1, var{j}(:,inactive)=0; end
    end
    p=project(var{num},Lapk);
    q=sum(p.^2,1); s=sum(var{num}.^2,1);
    active=find(s>0);
    var{num}(:,active)=var{num}(:,active).*sqrt(s(active)./q(active).^2);

    [datas,vars]=ALGOchoose(var,ngmar,aa,maxiteropt,Lapk,alpha,epsilon,delta,timeend);
    Fdatas{iter}=datas;
    if ~isfolder('Data_run'), mkdir('Data_run'); end
    save(fullfile('Data_run',[names{index(iter)},'.mat']),'datas','vars','Fdatas');
    valueplot(datas,timeend,index(iter));
end
