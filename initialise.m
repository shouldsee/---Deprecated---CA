    for i=1:length(scalars)
        eval([scalars{i} 'o=[];'])
       
    end
    if ~exist('p0','var');
        p0=0.5;
    end
   
    if ~exist('pop')
        pop=n;
    end
    if ~exist('work','var')
        work=0;
    end

    LS_input_log=zeros(warmup,n^2);
    MIs=zeros(1,4);
    S_flow=ones(n,n);
%     S_flow(500:600)=20;
    S_changeND_old=zeros(n,n,9,lag+1);
%     LS_input_old=zeros(n,n,lag+1);
    S_pop_old=zeros(n,n,lag+1);
    S_input_old=zeros(n,n,lag+1);
    S_dinput_old=zeros(n,n,lag+1);
    S_linput_old=zeros(n,n,lag+1);
    cells_old=zeros(n+2,n+2,lag+1);
    LS_input_old=zeros(lag,n^2);
    LS_dinput_old=zeros(lag+1,n^2);
    DS_input_old=zeros(n,n,lag+1,elenum.S_input);
    S_pop_lagged=zeros(n,n);
    S_pop_lagged_old=zeros(n,n);
    S_linput=zeros(n,n);
    SD_Tlinput=zeros(n,n);
    SI_input_old=zeros(n,n,lag+1);
    SI_dinput_old=zeros(n,n,lag+1);
    H_Tinput=zeros(1,lag+1);
    axs=[];
    inds=[];
    
    if ~allrand && ~randrule
    S_popo=S_poplist{k};
    else
    S_popo=rand(n,n)<0.5;
    end
    
% S_popo=reshape(S_popo,[],n,n);
cells=zeros(n+2,n+2);
cellsM=zeros(n+2,n+2,uni);
D_transio=[];
D_inputo=[];
S_defectT=zeros(size(cells));
% S_defect=randi([0 1],n,n);
dall=zeros(warmup,uni);
S_defect=rand(n,n)<0.01;
S_defect=zeros(n,n);
% S_defect=ones(n,n);
S_defect(14:16,14:16)=1;
S_defect=ones(n,n);
S_pop0=rand(n+2,n+2)<p0;
siz=size(cells);
[xc yc]=ndgrid(2:pop+1,2:pop+1);
id=sub2ind(siz,xc,yc);
cells(id)=S_pop0(id);
mass=0.5;
mass=1;
hedge=0:0.01:30;
% hedge=0.01:0.01:4;
% hedge=0:0.01:1;
hedge=-30:0.1:30;
% S_defect([450,650,250])=1;
% S_defect(500)=1;
% S_defect(14:16,14:16)=1;

% S_defectT(x,y)=S_defect;
% S_defectT(1,:)=S_defectT(x(end),:);
% S_defectT(end,:)=S_defectT(x(1),:);
% S_defectT(:,1)=S_defectT(:,y(end));
% S_defectT(:,end)=S_defectT(:,y(1));
% S_pop0=squeeze(S_popo(1,:,:));

