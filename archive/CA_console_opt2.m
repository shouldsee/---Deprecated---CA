

%% start GUI

%=============================================
% set(0,'DefaultFigureWindowStyle','docked')
figure(1)
%build the GUI
%define the plot button
plotbutton=uicontrol('style','pushbutton',...
   'string','Run', ...
   'fontsize',12, ...
   'position',[100,400,50,20], ...
   'callback', 'run=1;');

%define the stop button
erasebutton=uicontrol('style','pushbutton',...
   'string','Stop', ...
   'fontsize',12, ...
   'position',[200,400,50,20], ...
   'callback','freeze=1;');

%define the Quit button
quitbutton=uicontrol('style','pushbutton',...
   'string','Quit', ...
   'fontsize',12, ...
   'position',[300,400,50,20], ...
   'callback','stop=1;close;');

number = uicontrol('style','text', ...
    'string','1', ...
   'fontsize',12, ...
   'position',[20,400,50,20]);

dC_display = uicontrol('style','text', ...
    'string','0', ...
   'fontsize',12, ...
   'position',[20,450,50,20]);

random = uicontrol('style','pushbutton',...
   'string','Randomise', ...
   'fontsize',12, ...
   'position',[400,400,50,20],...
    'callback','random=1;');
   
uicontrol('style','text', ...
    'string','Update toggle', ...
   'fontsize',12, ...
   'position',[200,475,50,20]);
update_on=uicontrol('style','pushbutton',...
   'string','on', ...
   'fontsize',12, ...
   'position',[200,450,50,20],...
    'callback','update=1;');
update_off=uicontrol('style','pushbutton',...
   'string','off', ...
   'fontsize',12, ...
   'position',[150,450,50,20],...
    'callback','update=0;');
interrupt_button=uicontrol('style','pushbutton',...
   'string','interrupt', ...
   'fontsize',12, ...
   'position',[250,450,50,20],...
    'callback','interrupt=1;stop=1;');


%=============================================
% if ~exist('k')
%    
% end
k=8;

FIR=repmat('000000000',[100 1]);
FIR(8,:)='111101111';
FIR(9,:)='111111111';


[cmap]=colormap(colorcube);
n=30;
%CA setup
FIR=reshape(str2num(FIR(k,:).').',[3,3]);
x = (1:n)+1;
y = (1:n)+1;
z = zeros(n+2,n+2,1);
cells = z;
part = z;
width=20;



imh = image(cat(3,cells,z,z));
set(imh, 'erasemode', 'none')



%% initiate some extra plots (including frequencey visualiser, live entropy plot, live input histo)


maxstep=242;
%Main event loop
run = 1; %no wait for a draw 
freeze = 0; %wait for a freeze
interrupt = 0;
% window=zeros(n+2,n+2,width);
% window_linear=zeros(10,n^2);
% window_part=zeros(n,n,width);
stepnumber=0;
soupnumber=0;
booknumber=0;

% rulenumber=0;
% mhstda=[];
% hstda=zeros(length(binf),1);
% temp1=zeros(3,1);
% temp2=zeros(3,1);
% mmtemp=zeros(3,1);
% snrtemp=zeros(3,1);
% 

%% decomment here to enable batch processing 

% storeM=1;
% rrule=1;
% store=1;
% import=0;
% fname='tst';
% v=2;
% k=8;
% haulsize=5;
% haulnumber=0;



% lists={'sharplist','mmlist','rulelist','snrlist','apclist','welist','lfplist'};
% lists={'rrlist','detlist','entrlist','l1list','mmlist','popmlist','popmmlist'};
% lists={'rr1list','rr2list','mmlist'};
% lists={'timp1list','timp2list','timp3list','timp4list'};
% lists={'list1','list2','list3','list4','list5'};
% scalars={'St1_Ent','St2_Ent','H_input'};
% scalars={'hno','popo','popio','pt2o','pt2io','pt4io','pt4o','cello'};
% scalars={'S_pop','H_pop','H_input','H_iinput'};
% scalars={'H_pop','H_input','S_pop'};
% scalars={'D_input','H_input'};

% outnumber=length(lists);
%% default values for manual processing
rng('shuffle','twister')

if ~exist('stepmax')
stepmax=142;
end
if ~exist('soupmax')
soupmax=3;
end
if ~exist('update')
update=0;
end


if ~exist('rulemax')
rulemax=100000;
end

if ~exist('storeM')   
storeM=0;
end
if ~exist('rulenumber')
    rulenumber=0;
end
if ~exist('rrule')
    rrule=0;
end
if ~exist('store')
    store=0;
end
if ~exist('ruleB')
    ruleB=str2num(dec2base(randi([0 2^(k+1)],1),2,k+1).')';
    ruleS=str2num(dec2base(randi([0 2^(k+1)],1),2,k+1).')';
end
    if ~exist('rrule');
    rrule=0;
    end
if ~exist('freq')
    freq=repmat(z,1,1,3); 
end
if ~exist('merge')
    merge=0;
end
if ~exist('import_rule')
    import_rule=0;
end
if ~exist('start')
    start=1;
end
if ~exist('haulsize')
    haulsize=500;
end
if ~exist('detail')
    detail=0;
end
if ~exist('lag')
    lag=2;
end

if detail==1;
    for i=1:length(scalars)
    scalar=scalars{i};
    eval([scalar 'list={};'])
    end
end

    
if merge==0;

    if storeM==1;
        filename=fname;
%         filename=sprintf('%s%d',fname,haulnumber);
%         m=matfile(filename,'Writable',true);
%         m.stats=zeros(rulemax,outnumber*2);
%         m.rulelist=cell(rulemax,2);
    end

    if store==1
        stats=zeros(haulsize*soupmax,outnumber);
        rulelist=cell(haulsize*soupmax,2);
        souplist=cell(haulsize*soupmax,1);
    end

end

preallocate

%% main event loop
figure(1)
while ~interrupt;
%list_reset
stop= 0; %wait for a quit button push
for i=1:length(scalars)
    eval([scalars{i} 'o=[];'])
end

random=1;
record=0;

    if soupnumber==soupmax || start==1;
        start=0;
    if soupnumber==soupmax
        last=[mean(temp,1,'omitnan') max(temp,[],1,'omitnan')-min(temp,[],1,'omitnan')];
        str=num2str(last);
        disp(str)   
        if store==1 
        outind=(rulenumber-1)*soupmax+1:rulenumber*soupmax;
        stats(outind,:)=temp;
        rulelist(outind,:)=repmat({ruleB,ruleS},soupmax,1);
        souplist(outind,:)=souptemp;        
        if rulenumber==haulsize;
            if storeM==1 
                if haulnumber*haulsize*soupmax+1>65000
                    filename=[fname num2str(booknumber)];
                    booknumber=booknumber+1;
                    haulnumber=0;
                end
                corner=['A' num2str(haulnumber*haulsize*soupmax+1)];
           
            if haulnumber==0
                header=true;
            else
                header=false;
            end
            print2file
            haulnumber=haulnumber+1;
            rulenumber=0;
            stats=zeros(haulsize*soupmax,outnumber);
            rulelist=cell(haulsize*soupmax,2);
            souplist=cell(haulsize*soupmax,1);
            end
        end
        end
    end
    temp=zeros(soupmax,outnumber);  
    souptemp=cell(soupmax,1);
    rulenumber=rulenumber+1;
    if rulenumber==rulemax;
    interrupt=true;
    end
    if detail==1
        if rulenumber~=length(rulename)+1;
           update_rule
        else
           return
        end
    else
        update_rule
        cchange=rulecurr(1:elenum.S_input)~=rulecurr(mod((0:elenum.S_input-1)+9,18)+1);

    end
    
    try
   
    str=num2str(stats(rulenumber,:));
    disp(str)
    catch err
    end
        random=1;
    soupnumber=0;
    end
    
%%% randomise the soup, record to souptemp
    soupnumber=soupnumber+1;
    random=1;
    if (random==1)
        if soupnumber/soupmax<0.5
        cells = (rand(n+2,n+2)<0.5) ;
        else
        cells(5:25,5:25)=(rand(21,21)<0.5);
        end
        if exist('glider','var')
            if glider
            % put some gliders
            cells=zeros(n+2,n+2);
%             m=[1 7 16];
%             o=[6 12 18];
            m=2;o=6;
            for a=m
                for b =o
            cells([2 3]+a,2+b)=1;
            cells([2 4]+a,3+b)=1;
            cells(2+a,4+b)=1;
                end
            end
            end
        end
        cells(1,:)=cells(x(end),:);
        cells(end,:)=cells(x(1),:);
        cells(:,1)=cells(:,y(end));
        cells(:,end)=cells(:,y(1));
        window=zeros(n+2,n+2,width);
        stepnumber=0;
        record=0;
        random=0;
        s=0;
        soup=mat2rle(cells(x,y));
        souptemp{soupnumber}=soup;
    end
    
DS_input_old=zeros(n,n,lag+1,elenum.S_input);
S_input_old=zeros(n,n,lag+1);
S_pop_old=zeros(n,n,lag+1);
S_linput_old=zeros(n,n-1,lag+1);
S_pop_lagged=zeros(n,n);
S_pop_lagged_old=zeros(n,n);
S_linput=zeros(n,n-1);
SD_Tlinput=zeros(n,n);
LS_input_old=zeros(lag+1,n^2);
MIs=zeros(1,lag+1);
S_defect=ones(n,n);
S_flow=ones(n,n);


cchange=rulecurr(1:elenum.S_input)~=rulecurr(mod((0:elenum.S_input-1)+9,18)+1);

while (stop==0) 

    if (run==1)
        %nearest neighbor sum
        S_pop=cells(x,y);
        LS_pop=S_pop(:)';
%         D_pop=hist(S_pop,0:1);
%         D_pop=D_pop/n^2;
%         H_pop=sum(-D_pop.*log2(D_pop),'omitnan');
         S_popND=findneighbor(S_pop,neighbor.S_pop);

%     S_input=sum(S_popND,3)+9*S_pop;
        S_inputT=conv2(single(cells),repmat(FIR.S_input,1,1,size(cells,3)),'same');
        S_input=S_inputT(x,y);
        D_input=hist(S_input(:),0:sum(FIR.S_input(:)));
        HD_input=entropise(D_input);
        H_input=sum(HD_input);
        LS_input=S_input(:)';
        LS_input_old=cat(1,LS_input,LS_input_old(1:lag,:));
    DS_input=bsxfun(@eq,S_input,permute((1:elenum.S_input)-1,[1 3 4 2]));
    DS_input_old=cat(3,DS_input,DS_input_old(:,:,1:lag,:));
    S_vic=convn(DS_input_old,FIR.outtot_3da,'same');
    SD_vic=bsxfun(@times,DS_input_old,permute(S_vic,[1 2 3 5 4]));
    SD_vic=reshape(SD_vic,[],elenum.S_input,elenum.S_input);
    SD_vic=squeeze(sum(SD_vic,1))';
    SD_vic=(SD_vic+SD_vic')/2;
%     SD_Tinput=reshape(D_Tinput,elenum.S_input,[]);
    SP_vic=SD_vic/sum(SD_vic(:));
    SP_vic1=sum(SP_vic,2);
    SP_vic2=sum(SP_vic,1);
%     H_vic=sum(entropise(SP_vic(:)));
    H_vic1=sum(entropise(SP_vic1));
    H_vic2=sum(entropise(SP_vic2));
%     MI_vic=H_vic1+H_vic2-H_vic;
    SMI_vic=-log2((SP_vic1*SP_vic2)./SP_vic);
    SMI_vic(isnan(SMI_vic)|isinf(SMI_vic))=0;

%     neighbor.S_flow=neighbor.S_pop;
    % S_flow=ones(n,n);
    S_inputND=findneighbor(S_input,neighbor.S_flow);
%     S_flowND=findneighbor(S_flow,neighbor.S_flow);
    % nSMI_vic=SMI_vic/MI_vic;
    % bonddict=logsig(SMI_vic);
    bonddict=SMI_vic;
    id=sub2ind(size(bonddict),S_inputND+1,repmat(S_input+1,1,1,length(neighbor.S_flow)));
    bondND=bonddict(id);
    bondNND=findneighbor(bondND,-neighbor.S_flow);
    bondN=repmat(sum(bondNND,3),1,1,length(neighbor.S_flow));

    lossexponent=-bondN(:,:,1);
    std_LE=std(lossexponent(:));

%     S_flowloss=2.^(-bondN(:,:,1));


%     change0=rulecurr(S_input+1)~=rulecurr(max(S_input,1));
%     change1=rulecurr(S_input+1)~=rulecurr(min(S_input+2,18));
%     change01=cat(4,change0,change1);
%     
%     S_popNDV=bsxfun(@eq,S_popND,reshape(0:1,1,1,1,2));
%     S_change=sum(bsxfun(@times,S_popNDV,change01),4);
%     
%     S_defectND=findneighbor(S_defect,neighbor.S_pop);
%     S_defect=sum(cat(3,cchange(S_input+1),S_change).*cat(3,S_defect,S_defectND),3);
%     
% %     N_defect=2*sum(S_defect(:)>0)/sum(reshape(S_pop_old(:,:,1:2),1,[]));
%     P_defect=sum(S_defect(:)>0)/900;
%     H_defect=-log2(P_defect);
%     H_defect(isnan(H_defect)|isinf(H_defect))=0;
%     
%     HD_input=entropise(D_input);
%     H_input=sum(HD_input);
%     P_input=2^(-H_input);
    
    
        
%         HD_input=-D_input.*log2(D_input);
%         HD_input(isnan(HD_input))=0;     
%         sH_input=sH_input-H_input^2;
%         H_input=sum(HD_input);
%         sH_input=sH_input-H_input^2;
%         H_input=sum(-D_input.*log2(D_input),'omitnan');
%         sHD_input=D_input.*log2(D_input).^2;
%         sHD_input(isnan(sHD_input))=0;
%         sH_input=sum(sHD_input);
    
%         S_input_old=cat(3,S_input,S_input_old(:,:,1:lag));   
        
%         DS_input=bsxfun(@eq,input,permute((1:elenum.S_input)-1,[1 3 4 2]));
%         DS_input_old=cat(3,DS_input,DS_input_old(:,:,1:lag,:));
%         D_input_old=shiftdim(sum(reshape(DS_input_old(x,y,:,:),[],elenum.S_input),1),1);
%         H_input_old=sum(entropise(D_input_old(:)));
%         
%         S_vic=convn(DS_input_old,FIR.outtot_3d,'same');
%         S_vic=S_vic(x,y,:,:);
%         SD_vic=bsxfun(@times,DS_input_old(x,y,:,:),permute(S_vic,[1 2 3 5 4]));
%         SD_vic=reshape(SD_vic,[],elenum.S_input,elenum.S_input);
%         SSD_vic=shiftdim(sum(SD_vic,1),1)/2;
%         HSD_vic=entropise(SSD_vic);
%         HS_vic=sum(HSD_vic(:));
        
%         MIs(1)=H_input;
%     for i=1:lag
%         MIs(i+1)=mutualinfo(LS_input_old(1:i+1,:)',1/n^2.*ones(n^2,1));
%     end
%         S_pop_old=cat(3,S_pop,S_pop_old(:,:,1:lag));    
 


% 
%         MI_vic=2*H_input_old-HS_vic;
%         MI_vic2=2*H_input_old-HS_vic2;
       
        
        
%         S_pop_lagged=sum(S_pop_old.*repmat(pop_elenumV,n,n,1),3);
%         D_pop_lagged=hist(S_pop_lagged(:),(1:elenum.S_pop_lagged)-1);
%         D_pop_lagged=D_pop_lagged./sum(D_pop_lagged);
%         H_pop_lagged=sum(entropise(D_pop_lagged),'omitnan');
    
%     S_linput=conv2(S_pop_lagged,FIR.S_linput,'valid');
%     D_linput=hist(S_linput(:),(1:elenum.S_linput)-1);
%     D_linput=D_linput/sum(D_linput);
%     H_linput=sum(entropise(D_linput));
%     MI=2*H_pop_lagged-H_linput;
%     
        
        if stepnumber>1+lag   
%     Tinput=S_input+S_input_old(:,:,lag).*elenum.S_input;
%     D_Tinput=hist(Tinput(:),(1:elenum.S_input^2)-1);
%     SD_Tinput=reshape(D_Tinput,elenum.S_input,[]);
% 
%     H_Tinput=sum(entropise(D_Tinput));    
%     sHV=sum(entropise(sum(SD_Tinput,1)));
%     sHH=sum(entropise(sum(SD_Tinput,2)));
%     MI_input=sHV+sHH-H_Tinput;
%     
%     VSD_Tinput=SD_Tinput./repmat(sum(SD_Tinput,1),length(SD_Tinput),1);
%     VSD_Tinput=-VSD_Tinput.*log2(VSD_Tinput);
%     VSD_Tinput(isnan(VSD_Tinput))=0;
%     
%     HSD_Tinput=SD_Tinput./repmat(sum(SD_Tinput,2),1,length(SD_Tinput));
%     HSD_Tinput=-HSD_Tinput.*log2(HSD_Tinput);
%     HSD_Tinput(isnan(HSD_Tinput))=0;
% 
%     HV=sum(VSD_Tinput,1)';
%     HH=sum(HSD_Tinput,2);
%     
%     V_input=[HD_input' HV+HH HV-HH  (1:18)'];
%     V=V_input(:,[2,3 ,1])';
%     VV=combvec(V,V);
%     dVV=[VV(1:3,:)-VV(4:6,:)];        
%     wt=D_Tinput;
%     wt=wt/sum(wt);
%     wtdVV=dVV(1:3,:).*repmat(wt,3,1);
%     v2=mean(abs(wtdVV(1,:)));
%     v3=mean(abs(wtdVV(2,:)));
% %     v4=mean(wtdVV,2)';
    
    phaseX=H_input;
    phaseY=std_LE;
    phaseZ=std_LE;
%     phaseZ=MI_input;

%     NFlux=-cross(VVr,VVc);
%     NFlux=NFlux(3,:);
%     NFlux=NFlux.*D_Tinput;
%     NFlux=sum(NFlux);
%     ANFlux=sum(abs(NFlux));
% %     ANFlux=log2(ANFlux);
% %     NFlux=log2(abs(NFlux)+1)*sign(NFlux);
%     RFlux=sum(VV.^2,1)/2;
%     RFlux=RFlux.*D_Tinput;
%     RFlux=sum(RFlux);
%     dmat_V_input=squareform(pdist(V_input));
%     E_dmat_V_input=[mean(V_input(:,[2 3 1]),1) eigs(dmat_V_input,3)'];
%         H_bool=(V_input(:,2)>1)+2*(V_input(:,3)>1);
% D=zeros(2,2);
% % A=zeros(4,4);
% for i=[0 1]
%     for j=[0 1]
%         D(i+1,j+1)=sum(V_input(H_bool==i+2*j),1);
%     end
% end
% % for i=0:3;
% %     for j=0:3;
% %         r=V_input(H_bool==i,4);
% %         c=V_input(H_bool==j,4);        
% %         [r,c]=meshgrid(r,c);
% %         cmb=cat(2,r,c);
% %         cmb=reshape(cmb,[],2);
% %         indices=sub2ind(size(adj), cmb(:,1), cmb(:,2));
% %         A(i+1,j+1)=sum(adj(indices));
% %     end
% % end
% % L_A=A(:)'/sum(A(:));
% L_D=D(:)'/sum(D(:));
% DD=D(1,1)+D(2,2)-D(1,2)-D(2,1);
% H_D=sum(-L_D.*log2(L_D),'omitnan');
% % H_D(isnan(H_D))=0;

% 
%          MSD_Tinput=HSD_Tinput(1:elenum/2,:)+HSD_Tinput(elenum/2+1:end,:);
%         Ents=[HD_input; sum(MSD_Tinput,1)]';
%         St1_Ent=sum(dot(Ents,Ents,2));
%         St2_Ent=sum(Ents(:));
%     


        
        
        end
    
        for i=1:length(scalars);
            scalar=scalars{i};
            try
            eval(sprintf('%so=[%so;%s];',scalar,scalar,scalar))
            catch
            eval(sprintf('%so=[%so;0];',scalar,scalar))
            end
        end
        
        
        if record==1; %detect entropy periodicity
                hnw=H_inputo(end-1:-1:end-2*width);
                [bool,ind]=ismember(H_input,hnw);
                if bool && ind<=width; 
                    l=hnw(1:ind-1)==hnw(ind+1:2*ind-1);
                   s=sum(l)==length(l);
                end
                if s==1
                    stop=1;
                    for i=1:length(scalars)
                        scalar=scalars{i};
                        eval(sprintf('%so=%so(1:end-(2*ind-1),:);',scalar,scalar));

                    end
                     s=0;
                end
          
            
        end
%         try
%         S_input_oold=S_input_old;
%         end
        S_input_old=cat(3,S_input,S_input_old(:,:,1:end-1));
%         FIR=[1 2 3;4 5 6;7 8 9]*(sum(FIR(:))+1);
%         iinput=conv2(input,repmat(FIR,1,1,size(cells,3)));
%         S_iinput=reshape(iinput(x,y),1,[]);
        
%         FIR=[1,2,3;4,5,6;7,8,9]*(sum(FIR(:))+1);
%         SD_iinput=conv2(SD_input,repmat(FIR,1,1,size(SD_input,3)),'same');
%         S_iinput=reshape(SD_iinput(x,y),1,[]);
%         D_iinput=hist(S_iinput,sum(FIR(:)+1));
%         D_iinput=D_iinput/n^2;     
%         H_iinput=sum(-D_iinput.*log2(D_iinput),'omitnan');
%        
        
       
        cells=rulecurr(S_inputT+1);
        cells(1,:)=cells(x(end),:);
        cells(end,:)=cells(x(1),:);
        cells(:,1)=cells(:,y(end));
        cells(:,end)=cells(:,y(1));
%         c=cells(x,y);
        
        stepnumber = 1 + stepnumber;
        set(number,'string',num2str(stepnumber))

%       H_inputo=1;
        
        
        if stepnumber==42;
        record=1;
        end

        if update==1
%         hn=[hn(2:4*width) hnt zeros(1,width)];
        
        if ~isempty(KDE);
        S_input=window_part(:);
        [f,xi] = ksdensity(S_input(:));
        if exist('stdb')==0
            stdb=zeros(1,5*width);
        end
        
        stdb=[stdb(2:4*width) kurtosis(S_input) zeros(1,width)];
        end
            
        if exist('hn')==0
        hn=zeros(1,5*width);
        end
            set(f7,'YData',hn,'XData',1:5*width);
            set(KDE,'XData',xi,'YData',f)
            set(h2,'YData',stdb,'XData',1:5*width)
            set(imh, 'cdata', cat(3,cells(x,y),z(x,y),z(x,y)) )
            addpoints(f8,hn(4*width),hnt);
        end
        
        if stepnumber==stepmax;
            stop=1;
            stepnumber=stepnumber-(2*ind-1);
        end
        drawnow
        
    end
    %need this in the loop for controls to work
end
estep=stepnumber;

if detail==1 ;
     for i=1:length(scalars)
        scalar=scalars{i};
    if estep>40
                eval(sprintf('%slist=[%slist;{%so(20:end,:)}];',scalar,scalar,scalar))
        else
        eval(sprintf('%slist=[%slist;{repmat(0,[142 1])}];',scalar,scalar))
        end
    end
end
process_outputs

end
            %         if ~isempty(freq)
%         for s=1:v
%             freq(x,y,s)=sum(window(x,y,:)==s-1,3);
%         end
%         hft=hft/sum(hft);
%         hft=sum(-hft.*log2(hft),'omitnan');
%         hf(1)=[];
%         hf(end-width)=hft;
%         
%         hstda(1:end-1)=hstda(2:end);
% %         hstda(end-width)=sum(abs(C(:)));
%         hstda(end-width)=std(C(:));
%             if std(hno(end-20:end))<
%                     mm=maxhstda-minhstda;
%             mms=num2str(mm, '%10.5e\n')
%         end
%         set(a,'YData',hstda)
%         end
        


        

        
%  try
% %         refreshdata(5)
%         refreshdata(6)
%         set(h2,'YData',hf)
% %       imt(x,y)=cmap(freq(x,y,1)+1,:);
%          set(dC_display,'string',mms);
%         
%         set(imf,'CData',freq(x,y,2));
%  end



    
  %%
%  figure(9)
%  autocorr(diff(hno),min(160,length(hno)-10));


%%% start processing outputs

%  a=autocorr((hno),min(160,length(hno)-1));
%  b=autocorr(popo,min(160,length(popo)-1));
% Fs = 1;            % Sampling frequency
% T = 1/Fs;             % Sampling period
% L = length(a);             % Length of signal
% t = (0:L-1)*T;        % Time vector
% Y = fft(a);
% P2 = abs(Y/L);
% P1 = P2(1:L/2+1);
% P1(2:end-1) = 2*P1(2:end-1);
% f = Fs*(0:(L/2))/L;
%calculate sum of concentration/partial correlation
% input=diff(hno);

 
% d=pdist(apc);
% d=squareform(d);
% rp=d<thres;

% c=cov(apc);
% ic=pinv(c);
% e=entropyfilt(c);

% sharp=sharpness(squareform(pdist(hno')));
% sumconc=sum((e(:)));
% we=sum(abs(ic(:)));
% maxamp=max(P1);

% lfp=sum(P1(find(f<0.1)));
% temp1(soupnumber+1)=mm;
% temp2(soupnumber+1)=rr;
% mmtemp(soupnumber+1)=mm;
% snrtemp(soupnumber+1)=sumconc;
% apclist=cat(3,apclist,apc);

%         subplot(3,1,1)
%         autocorr(hno,min(160,length(hno)-1));
%         subplot(3,1,2)
%         plot(f,(P1))
%         title(  sprintf('entropy of spectrum=%f \n, max amp=%f, sum amp=%f',we,maxamp,sum(P1)))
%         xlabel('f (Hz)')
%         ylabel('|P1(f)|')
%         subplot(3,1,3)
%         Y = fft(b);
%         P2 = abs(Y/L);
%         P1 = P2(1:L/2+1);
%         P1(2:end-1) = 2*P1(2:end-1);
%         plot(f,(P1))
        %autocorr(a,min(160,length(a)-10));
%         plot(sort(P1,'descend'))

% pause




 %%% store the average

% pause



% ruleB=randi([0 8],[randi([0 8],1) 1]);
% ruleS=ruleB

% 
% %%
% mat.welist=[mat.welist;welist];
% mat.maxamplist=[mat.maxamplist;maxamplist];
% mat.mmlist=[mat.mmlist;mmlist];
% mat.rulelist=[mat.rulelist;rulelist];
% welist=[mat.welist;welist];
% maxamplist=[mat.maxamplist;maxamplist];
% mmlist=[mat.mmlist;mmlist];
% rulelist=[mat.rulelist;rulelist];

%%
% figure(10)
% % ind=1:length(snrlist);
% [~,ind]=sort(eval(slist),'descend');
% % welist=welist(ind);
% % lfplist=lfplist(ind);
% mmlist=mmlist(ind);
% sharplist=sharplist(ind);
% % rulelist=rulelist(ind,:);
% % snrlist=snrlist(ind);
% 
% scatter3(sharplist,mmlist,mmlist);
% s = [1:length(mmlist)]';s = num2str(s); s = cellstr(s);
% text(sharplist,mmlist,mmlist,s);
% 
% id=find(isnan(eval(list))==0);
% xlabel('wavelet entropy of fourier-transformed autocorrelated entropy');
% ylabel('snr of')
% zlabel('max-min entropy')