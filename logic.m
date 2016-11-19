

H_bool=V_input(:,2)>0.6+2*V_input(:,3)>0.6;
D=zeros(2,2);
for i=[0 1]
    for j=[0 1]
        D(i+1,j+1)=sum(V_input(H_bool==i+2*j),1);
    end
end
% adj;
A=zeros(4,4);
for i=1:max(H_bool);
    for j=1:max(H_bool);
        r=V_input(H_bool==i-1,4);
        c=V_input(H_bool==j-1,4);
        
        [r,c]=meshgrid(r,c);
        cmb=cat(2,r,c);
        cmb=reshape(cmb,[],2);

        indices=sub2ind(size(adj), cmb(:,1), cmb(:,2));

        A(i,j)=sum(adj(indices));
    end
end

figure(9)
subplot(2,2,1)
imagesc(D)
subplot(2,2,2)
imagesc(A)