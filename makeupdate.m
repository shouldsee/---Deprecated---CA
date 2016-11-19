if ~ishandle(2)
    figure(2)
%     clf
    subplot(2,2,1)
    f1=imagesc(S_input+1);
    caxis([0 1])
    view(-90,90)
    colorbar      
    caxis([1 elenum.S_input+1])
%     view([90 -90])
%     legend(cellstr(num2str((1:elenum)')))
    subplot(2,2,3)
%     f3h=histogram2(HIVxs,HIVys)
%     ,[m1 20]);
    
    zlim([1 n^2])
    set(gca,'ZScale','log')
%     colormap(gca,'colorcube')
%     f3q=quiver(pdata(1,:),pdata(3,:),pdata(4,:),pdata(6,:),0);
%     f3q=quiver3(pdata(1,:),pdata(2,:),pdata(3,:),pdata(4,:),pdata(5,:),pdata(6,:),0);
%     markers = {'+','o','*','.','x','s','d','^','v','>','<','p','h'};
%     marker=cell2mat(markers([1:9 1:9]));
%     f3=gca;
%     xlim([0,4]) 

%     ylim([0 0.75])
%     ylim([-4 4])
%     zlim([0 0.75])

%     f3s=scatter(f3p(:,1),f3p(:,2),3,'x');
%         xlim([0 5])
%     ylim([0 1]*5)

%     f3i=imagesc(S_flow);
%     view(-90,90)
%     caxis([0 15])

%     hold on
%     f3g=gscatter(V(1,:),V(2,:),1:elenum.plot,colormap(gca),marker);
%     legend(cellstr(num2str((0:elenum.plot)')),'Location','bestoutside')
    hold off
%     refline([0 0])

    subplot(2,2,4)
%     f4i=imagesc(VSD_Tinput);
%     f4i=imagesc(ldCS_input);
%     f4i=imagesc(fHSD_vic);
%     f4i=imagesc(St_input_MI);
%     f4i=imagesc(SI_input);
%     f4i=imagesc(S_defect);
        view(-90,90)
        colorbar
%     caxis([0,5])
%     ax4=gca;
%     ax4.XTickLabels=num2str((ax4.XTick-1)')';
    
    subplot(2,2,2)
%     ys=std(H_inputV_scale,1);
%     xs=1:length(ys);
%     f2s=scatter(xs,ys);
%     f2l=line(xs,ys);
    ylim([0 0.5])
    xlim([0 m1])
%     f2h=histogram(S_flow(:),hedge);
%     set(gca,'Yscale','log')
%     set(gca,'Xscale','log')
%     ylim([1 600])
%     f2i=imagesc(SMI_vic);
%     f2i=imagesc(Hy_x);
%     f2i=imagesc(SI_dinput);
%     view(-90,90)
%     caxis([0,6])
%     f2=scatter3(V_input(:,2),V_input(:,3) ,V_input(:,1),5,(1:elenum.plot)*3,'x');
%     
%     colormap(gca,'colorcube')
%     hold on
%     axis([0 4 0 4 0 0.5 ]);
%     view(50,90);
% 
%     xlabel('Forward uncertainty')
%     ylabel('Backward uncertainty')
%     zlabel('Entropic weight')
    rind=floor((k-1)/soupmax)+1;
    suptitle(['MISC' ' rulenumber='  num2str(rind) ' ' rulename{rind,:}])
    elseif mod(stepnumber,2)==0 
        set(f1,'CData',cellsM(x,y,1))
% %         set(f4i,'CData',bonddict)
%         set(f4i,'CData',lossexponent)
% %         set(f2i,'CData',SMI_vic);
% %         set(f2h,'Data',lossexponent(:));;
%         set(f2s,'YData',stds);
%         set(f2l,'YData',stds);
%         
% %         set(f3i,'CData',S_flow);
%           set(f3h,'Data',[HIVxs,HIVys]);
%         set(f1,'CData',S_linput);
%         set(f3s,'XData',f3p(:,1),'YData',f3p(:,2))
%         set(f3q,'XData',pdata(1,:),'YData',pdata(3,:),'UData',pdata(4,:),'VData',pdata(6,:));
%         set(f3q,'XData',pdata(1,:),'YData',pdata(2,:),'ZData',pdata(3,:),'UData',pdata(4,:),'VData',pdata(5,:),'WData',pdata(6,:));
%         set(n,ni,'CData',St_input_MI);
%         set(f4i,'CData',MI_vic1(S_input+1))
%         hist(S_flow(:),0.0.1)
%         for i=1:length(f3g);
% %             set(f3g(i),'XData',V(1,i),'YData',V(3,i));
%             set(f3g(i),'XData',V(1,i),'YData',V(2,i),'ZData',V_input(i,1));
%         end
%         subplot(2,2,2)
%         scatter3(V_input(:,2),V_input(:,3) ,V_input(:,1),5,(1:elenum.plot)*3,'x');
% %                                 view(50,90);
%                     colormap(gca,'colorcube')
%                            view(50,90);
%     hold on
    end
    drawnow
