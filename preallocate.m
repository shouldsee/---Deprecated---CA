clear elenum
global FIR
clear FIR
if ~exist('n','var')
    n=30;
end
x = (1:n)+1;
y = (1:n)+1;
%% FIRs & elenums
elenum.S_pop=2;
a=-1:1;
a=combvec(a,a)';
sa=sum(abs(a)==1,2);
a=a((sa>0),:);
b=-2:2;
b=combvec(b,b)';
% sa=sum(abs(b)==1,2);
% a=a((sa>0),:);
neighbor.n5=b;
neighbor.S_pop=a;
neighbor.S_flow=neighbor.S_pop;
neighbor.S_change=[0 0;a];
d=pdist2(neighbor.S_change,-neighbor.S_change,'hamm');
[neighbor.project,b]=find(d==0);

% neighbor.S_pop=a;

FIR.S_input=[1 1 1;1 9 1;1 1 1];
FIR.outtot=[1 1 1; 1 0 1;1 1 1]/8;
FIR.outtot_3d=cat(3,ones(3,3),[ 1 1 1; 1 0 0; 0 0 0],zeros(3,3));
FIR.outtot_3da=cat(3,zeros(3,3,2),ones(3,3));
FIR.S_dinput=2.^[0 1 2;3 8 4;5 6 7];
FIR.S_dinput(FIR.S_dinput<1)=0;
FIR.dm=elenum.S_pop.^reshape(0:dm^2-1,dm,dm);
FIR.Sinput=FIR.S_input;


elenum.S_input=sum(FIR.S_input(:))+1;
elenum.S_pop2input=(elenum.S_pop-1)+(elenum.S_input-1)*elenum.S_pop+1;

elenum.S_dinput=sum(FIR.S_dinput(:)*1)+1;
elenum.S_tot2d=(elenum.S_input-1)*1+(elenum.S_dinput-1)*elenum.S_input+1;


% FIR.outtot_3d=ones(3,3,3);
% FIR.outtot_3d=padarray(ones(2,2,2),[1 1 1],0,'pre');
% FIR.outtot_3d(2,2,2)=0;


%     pop_elenumV=permute(1.^(lag:-1:0),[1,3,2]);
       pop_elenumV=permute(2.^(lag:-1:0),[1,3,2]);
 
    elenum.S_pop_lagged=sum(pop_elenumV)+1;
FIR.S_linput=reshape(elenum.S_pop_lagged.^([0 1 2 3]),2,2);
%     FIR.S_linput=elenum.S_pop_lagged.^[0 1];
    elenum.S_linput=sum(FIR.S_linput(:).*(elenum.S_pop_lagged-1))+1;     


SD_Tlinput=zeros(n,n);
FIR.vic=makeFIR(elenum.S_input,zeros(2,2,2));
%% dicts
% dict.S_linput=maketable(elenum.S_pop_lagged.^([0 1;2 3]),elenum.S_pop_lagged);
% elenum.S_linput=max(dict.S_linput);

%% change
% 1:elenum.S_input
sinputs=0:elenum.S_input-1;
% sinputs=0:elenum.S_input-1;
change1=rulecurr(sinputs+1)~=rulecurr(max(sinputs,1));
change0=rulecurr(sinputs+1)~=rulecurr(min(sinputs+2,18));
change2=rulecurr(sinputs+1)~=rulecurr(mod(sinputs+9,18)+1);
change=cat(1,change0,change1,change2)';

% change0=rulecurr(sinputs+1)~=rulecurr(max(sinputs,1));
% change1=rulecurr(sinputs+1)~=rulecurr(min(sinputs+2,18));
% change2=rulecurr(sinputs+1)~=rulecurr(mod(sinputs+9,18)+1);
% change=cat(1,change0,change1,change2)';
% change01=cat(4,change0,change1);

% 
        
%%
%% makedict
ndict=[];
compdict=[];
for s=0:2^9-1
    nb=dec2base(s,2,9);
    nb=str2num(nb');
    ndict(:,:,s+1)=reshape(nb,3,3);
    compdict(s+1)=sum(nb);
end


