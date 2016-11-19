pre_init
preallocate
% ruletemp=importrule(rulename);
kmax=length(rulename)*soupmax;
if randrule
    kmax=100000000;
end

for k=k0:kstep:kmax

    D=zeros(stepmax,stepmax);
    
    for soup=1:soups
        initialise

    if sum(S_popo(:))==0
            for i=1:length(lists)
               list=lists{i};
               eval(sprintf('%ss=[%ss;zeros(1,size(%ss,2))];',list,list,list))
            end

            continue
    end

    getrule
%         rind=floor((k-1)/soupmax)+1;
        if update2==1
        figure(1)
        subplot(2,2,1)
        f1=imagesc(cells(x,y),'AlphaData',1);
        caxis([0 1])
        subplot(2,2,3)
        f2=imagesc(-1+cells(x,y),'AlphaData',1);
        subplot(2,2,2)
        f3s=scatter([],[],3,'x');
        title(sprintf('dm=%d', dm) )
        xlim([0 1])
        ylim([0 1])
        subplot(2,2,4)
                f4s=scatter([],[],3,'x');
%         title(sprintf('dm=%d', dm) )
        xlim([0 1]*4)
        ylim([0 1]*4)

%         f4h=histogram([],0:0.01:1);
%         set(gca,'Yscale','log')
%         ylim([1 n^2*2])
%         refline([1 0])
%         colorbar
        suptitle(['mean-field p map' ' rulenumber='  num2str(rind) ' ' rulename{rind}])
        
        end


        for step=1:stepmax
        S_inputT=conv2(single(cells),FIR.S_input,'same');
        S_input=S_inputT(x,y);

        cells(x,y)=rulecurr(S_input+1);
    %     [xs,ys]=ind2sub([n,n],);
 
    %     if ismember(mean(S_pop(:)),[0,1])  
        if mod(step,rstep)==0 || ismember(mean(mean(cells(x,y))),[0 1])
            cells(x,y)=rand(n,n)<rand(1,1);
        end
               ind=randperm((n/dm)^2);
        S_pop=cells(x,y);

if mod(step,sstep)==0
    
        if dm>1
        S_convT=conv2(cells,FIR.dm,'same');
        S_conv=S_convT(x,y);
        [m1,m2]=size(S_conv);
        S_convt=reshape(S_conv,dm,m1/dm,dm,m2/dm);
        

        S_convtt=S_convt(dmm,:,dmm,:);       
        S_convt(dmm,:,dmm,:)=reshape(S_convtt(ind),1,n/dm,1,n/dm);
        for x1=1:m1/dm;
            for x2=1:m2/dm;
                S_convt(:,x1,:,x2)=reshape(str2num(dec2base(S_convt(dmm,x1,dmm,x2),2,dm^2)')',dm,dm);
            end
        end
        S_conv=reshape(S_convt,n,n);
        cells(x,y)=S_conv;

        else 
        if dm==1
        cells(x,y)=reshape(S_pop(ind),n/dm,n/dm);
        end
        end  

        
end

        cells(1,:)=cells(x(end),:);
        cells(end,:)=cells(x(1),:);
        cells(:,1)=cells(:,y(end));
        cells(:,end)=cells(:,y(1)); 
        LS_pop=reshape(cells(x,y),1,[]);
        cells_old=cat(3,cells,cells_old(:,:,1:lag));  
        
        
    %     df=sum(cells_old(x,y,1:lag),3)-sum(cells_old(x,y,2:lag+1),3);
%         df=mean(diff(cells_old(x,y,:),[],3),3);
    %     df=mean(cells_old(x,y,1:lag),3);
        %     df=cells_old(x,y,1)-cells_old(x,y,3);
%         df=abs(df);
        D_input=hist(S_input(:),0:elenum.S_input-1);
        H_input=sum(entropise(D_input));
%         H=sum(entropise(S_convtt(:)));
        rec=[mean(mean(cells(x,y))) H_input];
% rec=[H mean(df(:))];
%         rec=[std(df(:)) mean(df(:))];
         for i=1:length(scalars); 
            scalar=scalars{i};
            eval(sprintf('%so=[%so;%s];',scalar,scalar,scalar))
         end
%       
        if step>lag+1
            if update2==1
           set(f1,'CData',df);
            set(f2,'CData',cells(x,y))
             set(f3s,'XData',[f3s.XData reco(step-2,1)],'YData',[f3s.YData,reco(step-1,1)]);
             set(f4s,'XData',[f4s.XData reco(step-2,2)],'YData',[f4s.YData,reco(step-1,2)]);

             if step>lag+2
%         set(f3s,'XData',[f3s.XData reco(step-lag-2,1)],'YData',[f3s.YData,reco(step-lag-1,2)]);
%         set(f4h,'Data',df(:))
        end
            end
            if mod(step,2)==0
        drawnow
             end
        end

        step
    %     pause(0.2)
        end
soup
%         subplot(2,2,4)
%         histogram(reco(:,1),0:0.01:1)

% rind=ceil((k)/soupmax)    
% if rind>size(rulename,1);
%     rind=size(rulename,1)
%     rind=1;
% end
%     figure(2)   
    D=D+squareform(pdist(LS_popo,'hamming'));
    
    end
    D=D/soups;
    if update2==1
     drawnow
% 
%     figure(2)
%     subplot(2,2,1)
%     imagesc(D)
% %     caxis([0 1])
%     subplot(2,2,2)
%     plot(D(:,1))
%     ylim([0 1])
%     subplot(2,2,3)
%     [x1 x2]=size(D);
%     A=reshape(D,2,x1/2,2,x2/2);
%     C=squeeze(A(1,:,1,:))
%     imagesc(C)
    suptitle(['mean-field p map' ' rulenumber='  num2str(rind) ' ' rulename{rind}])
    end
    if store
            makeprint
        end

    
end

%%
% [m n ]=size(D);
% A=reshape(D,2,m/2,2,n/2);
% C=squeeze(A(1,:,1,:))