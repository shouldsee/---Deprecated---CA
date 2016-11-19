 function[VDV,udiF_WtdVDV]=calcflux(V_input,SD_Tinput)
elenum=length(V_input);
 X=V_input(:,1:3);
% distfun = @(XI,XJ) sum(bsxfun(@minus,XI,XJ),2); % swap XI,XJ if needed
% D = squareform(pdist(X,distfun));
V=V_input(:,[2,3,1])';
VV=combvec(V,V);
VV=VV(4:6,:)-VV(1:3,:);
VDV=reshape(VV,3,elenum,elenum);
% SD_Tinput=reshape(D_Tinput,[],elenum);
Wt=repmat(permute(SD_Tinput,[3 1 2]),3,1,1);
Wt=2.^(Wt);
WtdVDV=Wt.*VDV;

% diF_WtdVDV=sum(sum(WtdVDV,3),2)';
Abs_VDV=squeeze(sum(VDV(1:2,:,:).^2,1).^0.5);
udiF_WtdVDV=Abs_VDV.*SD_Tinput;

% [SD_Tinput(6,16) Wt(1,6,16)]
% Wt=repmat(SD_Tinput,)
% imagesc(SD_Tinput);
end