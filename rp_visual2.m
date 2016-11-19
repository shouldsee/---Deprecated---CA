

list=S_poplist;
m=length(list);
timpms=[];
for rind=1:floor((m)/soupmax)+1
    figure(8)
    clf
%     fprintf('mmin=%.5f rr=%.5f,rulenumber=%d',rr1,rr2,rind);
% rind=floor((k)/soupmax)+1    
suptitle(['MISC' ' rulenumber='  num2str(rind) ' ' rulename{rind,:}])
    cellos=[];
    cs=zeros(4,4);
    cs=[];
    bs=[];
    timps=[];
    for k=(soupmax*(rind-1)+1):soupmax*rind
%         for k=soupmax*(rind-1)+1
    S_pop=S_poplist{k}';
    if sum(S_pop(:))~=0
        x=movmean(H_poplist{k}',6);
        y=movmean(H_inputlist{k}',6);
        z=movmean(H_iinputlist{k}',6);
        a=movmean(H_iinputlist{k}',6);

        x=double(x);
        y=double(y);
        z=double(z);
        a=double(a);
        if isempty(x);
            x=zeros(1,144);
            y=zeros(1,144);
            z=zeros(1,144);
            a=zeros(1,144);
            S_pop=zeros(900,144);
        end
        t=1:length(x);
        
        dx=diff(x);
        dy=diff(y);
        dz=diff(z);
        da=diff(a);
        dt=diff(t);
if ~exist('vdim')
    vdim=2;
end
switch(vdim)
    case 2
        subplot(2,3,1)
        quiver(x,y,[dx 0],[dy 0],'MaxHeadSize',5)
        text(x(1),y(1),num2str(k))
%         line(x,y,z,'Linewidth',1)
        hold on
        subplot(2,3,2)
        quiver(x,y,[dx 0],[dy 0],'MaxHeadSize',5)
        text(x(1),y(1),num2str(k))
        hold on
    case 3
        subplot(2,3,1)
%         scatter3(x,y,z,'Linewidth',1)
%         quiver(x,y,[dx 0],[dy 0],'MaxHeadSize',5)
        plot(x,y)
        text(x(1),y(1),num2str(k))
%          quiver3(x,y,z,[dx 0],[dy 0],[dz 0],'MaxHeadSize',5)
%          text(x(1),y(1),z(1),num2str(k))
        hold on
        subplot(2,3,2)
        plot3(x,y,z)
%         line(x,y,z)
%          quiver3(x,y,z,[dx 0],[dy 0],[dz 0],'MaxHeadSize',5)
%          text(x(1),y(1),z(1),num2str(k))
             hold on
end

subplot(2,3,4)
% plot(real(fft(x)))
hold on
fx=abs(fft(x));
fx(1)=[];
plot(fx)

subplot(2,3,5)
% plot(real(fft(y)))
hold on 
fy=abs(fft(y));
fy(1)=[];
plot(abs(fy))
subplot(2,3,6)
% plot(real(fft(z)))
hold on
fz=abs(fft(z));
fz(1)=[];
plot(abs(fz))
%         scatter3(dx,[dx(2:end) 0],(3*k)*ones(1,length(dx)));
%         scatter3(dy,[dy(2:end) 0],(3*k+1)*ones(1,length(dx)));
%         scatter3(dz,[dz(2:end) 0],(3*k+2)*ones(1,length(dx)));
                 

            cellos=[cellos;S_pop'];

        dxt=dx/sum(dx);
        dxh=sum(-dxt.*log2(dxt),'omitnan');
        
%         hist(dxt);
         c=cov([dx;dy;dz;da;t(1:end-1)]');
        invc=diag(diag(c))^-0.5;
        c=invc*c*invc;
        
%         cs=(cs+c);

        cs=[cs;c];
        b=(xcorr([x;y;z]'));
        bs=[bs;b];

        fprintf('\n %d ',k)
        fprintf(['dxh=' num2str(mean(dx))])

        data=[dx;dy;dz;da]';
        dtdata=detrend(data);
        
        timp=1-var(dtdata,1)./var(data,1);
        timps=[timps;mean(fx),mean(fy),mean(fz)];
        end
        
    end
    timpms=[timpms;mean(timps,1)];
    subplot(2,3,1)
    hold off
    subplot(2,3,2)
        axis([0 1 0 10 0 10])
    hold off
    subplot(2,3,3)
    imagesc(cs)
%     subplot(2,1,2)
    imagesc(cellos)
    
%     waterfall(timps')
    zlim([0 0.5])
    sumcov=mean(cs(:),'omitnan');
    title(sprintf('sum-cov=%f',sumcov))
    caxis([0 1])
    colorbar
    subplot(2,3,4)
%     plot(fft(H_pop))
    ylim([-10 10])
    hold off
    subplot(2,3,5)
    ylim([-100 100])
%     plot(fft(H_input))
    hold off
    subplot(2,3,6)
    ylim([-100 100])
%     plot(fft(H_iinput))
    hold off
%     figure(10)
%     [cr,lgs] = xcorr([dx;dy;dz]','coeff');

%     for row = 1:3
%         for col = 1:3
%             nm = 3*(row-1)+col;
%             subplot(3,3,nm)
%             stem(lgs,cr(:,nm),'.')
%             title(sprintf('c_{%d%d}',row,col))
%             ylim([-1 1])
%         end
%     end
%     subplot(2,3,4);
%     autocorr(dx)
%     subplot(2,3,5)
%     autocorr(dy)
%     subplot(2,3,6)
%     autocorr(dz)
%      spy(cellos);
%     imagesc(cellos);
%     d=squareform(pdist([x;y;z]'));
%     subplot(2,3,2)
%     imagesc(d)
%     imagesc(bs)

    
    figure (9)
    clf
    cs=[];
    imagesc(cellos); 
    suptitle([listname ' rulenumber='  num2str(rind) ' ' rulename{rind,:}])
%      for k=(soupmax*(rind-1)+1):soupmax*rind
%         S_pop=cellolist{k}';
%         if sum(S_pop(:))~=0
%             cellos=[cellos;S_pop];
%         end
%          
%         c=cov([x;y;z;a]');
%         invc=diag(diag(c))^-0.5;
%         c=invc*c*invc;
%         
% %         cs=(cs+c);
%         cs=[cs;c];
%         disp(k)
%     end
%    
%     pause
end

%%
figure(11)
waterfall(timpms);
