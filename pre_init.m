if ~exist('n','var')
n=50;
end
if ~exist('stepmax','var')
stepmax=500;
end

kstep=6;
if ~exist('po','var')
po=0:0.05:1;
end
% warmup=length(po);
if ~exist('warmup','var')
    warmup=10;
end
if ~exist('store','var')
    store=1;
end

if ~exist('soupmax','var')
soupmax=6;
end

if ~exist('randrule','var')
    randrule=0
end
    S_poplist=repmat({rand(n,n)<0.5},length(ruletemp)*soupmax,1);
% warmup=10;
global x y

x=2:n+1;
y=2:n+1;
eso=[];
egs=[];
v1s=[];
v2s=[];
timp_D_inputo=[];
timp_H_popo=[];
H_ego=[];
inds=[];
E_dmat_V_inputos=[];
lists={'v1','v2'};
scalars={...'H_input','H_pop','D_input',...
    'rec',...
    'LS_pop',...
   ...'stds',...'HLE','H_flow','LS_flow','D_flow','MI_vic','LS_pop','MIs',...'FH_input'...,'LS_defectL','N_defect',...'H_input',...
    ...'H_input_old','S_inputp',...
    ...'D_pop_lagged','HD_pop_lagged'...
...    'H_dinput_input',
...    'H_dinput','MI_dinput',...
...%     'sH_input','Hx2','Hxy2',
...%    'HH_Tinput',...
...%     'HD_input',...
...%     'H_pop','D_Tinput','VSD_Tinput','V_input','va','vr','H_Tinput',...
...     'MI_pop',...
%      'LSI_input',...'LS_pop','LS_input','HS_vic2','HSS_vic2',...
...%     'KL_div'
     };
% cmap=colormap('colorcube');
if ~exist('dirname')
    dirname='gallery/temp/';
end
mkdir(dirname)
fprintf('saveing to %s \n',dirname);
trigger=1;
if ~exist('k0')
    k0=1;
end

if ~exist('dm')
    dm=1;
end

dmm=floor((dm+1)/2);
if ~exist('lag')
        lag=1;
end

if ~exist('randrule')
    randrule=0;
end
if ~exist('allrand','var')
    allrand=1;
end
