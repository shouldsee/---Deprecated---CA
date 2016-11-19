
    DS_input=bsxfun(@eq,S_input,permute((1:elenum.S_input)-1,[1 3 4 2]));
    DS_input_old=cat(3,DS_input,DS_input_old(:,:,1:lag,:));
%%    
%     S_inputV=bsxfun(@eq,S_input,shiftdim((1:elenum.S_input)-1,-1));
    S_inputV=bsxfun(@eq,S_pop,shiftdim(0:1,-1));
    repS_inputV=repmat(S_inputV,3,3);
    m1=size(S_inputV,1);
    m2=size(S_inputV,2);
    m1=8;
%     cx=m1+1:2*m1;
    cy=m2+1:2*m2;
    cx=cy;
    
    S_inputV_scale=[];
    for i=1:min(m1,m2);
        FIR.temp=ones(i,i);
        conv_out=convn(repS_inputV,FIR.temp,'same');
        c_conv_out=conv_out(cx,cy,:);
        mcy=find(mod(1:m2,i)~=0);
        mcx=mcy;
%         mcy
        c_conv_out(mcx,mcy,:)=nan;
        
        S_inputV_scale(:,:,:,i)=c_conv_out;
    end
    %%
    HS_inputV_scale=bsxfun(@rdivide,S_inputV_scale,shiftdim((1:min(m1,m2)).^2,-2));
    HD_input_V_scale=HS_inputV_scale.*-cslog2(HS_inputV_scale);
%     HD_input_V_scale(isnan(HD_input_V_scale))=0;
%     imagesc(H_input_V_scale)
%     H_inputV_scale=reshape(sum(HD_input_V_scale,3),[],min(m1,m2));
    H_inputV_scale=reshape(bsxfun(@rdivide,S_inputV_scale(:,:,1,:),sum(S_inputV_scale,3)),[],min(m1,m2));
%     ys=H_input_V_scale;
    HIVxs=repmat(1:size(H_inputV_scale,2),n^2,1);
    HIVys=reshape(H_inputV_scale,[],1);
    HIVxs=HIVxs(:);
    HIVys=HIVys(:);

    xs=1:n;
    ys=std(H_inputV_scale,1,'omitnan');
    std_H_inputV_scale=std(H_inputV_scale,1,'omitnan');
    stds=std_H_inputV_scale;
%     means=mean(H_inputV_scale)
%     scatter(xs,ys);
%     line(xs,ys)
% plot(xs,ys)
%%
%     S_vic=convn(DS_input_old,FIR.outtot_3d,'same');
%     SD_vic=bsxfun(@times,DS_input_old,permute(S_vic,[1 2 3 5 4]));
%     SD_vic=reshape(SD_vic,[],elenum.S_input,elenum.S_input);
%     SSD_vic=shiftdim(sum(SD_vic,1),1)/2;
% %     SSD_vic=(SSD_vic+SSD_vic')/2;
%     HSD_vic=entropise(SSD_vic);
%     HS_vic=sum(HSD_vic(:));
%     fHSD_vic=triu(HSD_vic);
    
    S_vic=convn(DS_input_old,FIR.outtot_3da,'same');
    SD_vic=bsxfun(@times,DS_input_old,permute(S_vic,[1 2 3 5 4]));
    SD_vic=reshape(SD_vic,[],elenum.S_input,elenum.S_input);
    SD_vic=squeeze(sum(SD_vic,1))';
    SD_vic=(SD_vic+SD_vic')/2;
%     SD_Tinput=reshape(D_Tinput,elenum.S_input,[]);
    SP_vic=SD_vic/sum(SD_vic(:));
    SP_vic1=sum(SP_vic,2);
    SP_vic2=sum(SP_vic,1);
    H_vic=sum(entropise(SP_vic(:)));
    H_vic1=sum(entropise(SP_vic1));
    H_vic2=sum(entropise(SP_vic2));
    MI_vic=H_vic1+H_vic2-H_vic;
%     MI_vic1
%     MI_vic=H_vic1
    SMI_vic=-log2((SP_vic1*SP_vic2)./SP_vic);
    SMI_vic(isnan(SMI_vic)|isinf(SMI_vic))=0;
%     ESMI_vic=sum(sum(SMI_vic.*SP_vic));
    
    MI_vic1=sum(SMI_vic.*SP_vic,2)./SP_vic1;
    MI_vic1(isnan(MI_vic1))=0;
    SSMI_vic1T=MI_vic1(S_inputT+1);
    SSMI_vic1=SSMI_vic1T(x,y);
    SSMI_vic1_inputT=conv2(SSMI_vic1T,FIR.outtot,'same');
    SSMI_vic1_input=SSMI_vic1_inputT(x,y);
%     scatter(SSMI_vic1,SSMI_vic1_input,3,'x');

