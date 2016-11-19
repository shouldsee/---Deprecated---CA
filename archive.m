%%archive

%         cello=[cello cellt(:)];
%         part(x,y)=total(x,y);
%         x=2:n+1;y=2:n+1;
%         pt2= conv2(single(cells),[1 1 ; 1 1],'same');
%         pt2=pt2(2*x(1:floor(n/2)),2*y(1:floor(n/2)));
%         pt3= conv2(single(cells),[1 1 1;1 1 1;1 1 1],'same');
%         pt3=pt3(3*x(1:floor(n/3)-1),3*y(1:floor(n/3))-1);
%         pt4= conv2(single(cells),[1 1 1 1;1 1 1 1; 1 1 1 1; 1 1 1 1],'same');
%         pt4=pt4(4*x(1:floor(n/4)),4*y(1:floor(n/4)));
%         ptt=input(x,y);
%         cellt=cells(x,y);
%         ptt=ptt+cellt*9;
%         
%         hnt=hist(ptt(:),0:17);
%         hnt=hnt/sum(hnt);
%         hnt=sum(-hnt.*log2(hnt),'omitnan');
%         hno=[hno hnt];
%         
%         popt=hist(cellt(:),0:v-1);
%         popt=popt/sum(popt);
%         popt=sum(-popt.*log2(popt),'omitnan');
%         popo=[popo popt];
%         
%         popit=cellt(:)/sum(cellt(:));
%         popit=sum(popit.*log2(popit),'omitnan');
%         popio=[popio popit];
 
        
        
%         pt2it=pt2(:)/sum(pt2(:));
%         pt2it=sum(-pt2it.*log2(pt2it),'omitnan');
%         pt2io=[pt2io pt2it];
%         
%         pt4it=pt4(:)/sum(pt4(:));
%         pt4it=sum(-pt4it.*log2(pt4it),'omitnan');
%         pt4io=[pt4io pt4it];
%         
%         pt4t=hist(pt4(:),0:17);
%         pt4t=pt4t/sum(pt4t);
%         pt4t=sum(-pt4t.*log2(pt4t),'omitnan');
%         pt4o=[pt4o pt4t];
%         
%         
%         pt2t=hist(pt2(:),0:9);
%         pt2t=pt2t/sum(pt2t);
%         pt2t=sum(-pt2t.*log2(pt2t),'omitnan');
%         pt2o=[pt2o pt2t];

%         old update function
%         cells(x,y-1) +cells(x,y+1) + ...
%             cells(x-1, y) + cells(x+1,y) + ...
%             cells(x-1,y-1) + cells(x-1,y+1) + ...
%             cells(x+1,y-1) + cells(x+1,y+1);
%         % The CA rule
%         cells = ismember(part,rule) | (part==2 & cells);       

        %draw the new image

        %update the step number diaplay


        
%         popo=[popo mean(mean(cells(x,y)))];