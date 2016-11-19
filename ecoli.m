keyset={'A','T','C','G'};
valueset=[1,2,3,4];
mapObj=containers.Map(keyset,valueset);
%%
input=ecoliG;
%%
out=arrayfun(@(a) mapObj(a),ecoliG);
ecoliGN=out;
%%
pwmax=6;
wd=4^pw+pw-1;

% genome=randsample(ecoliGN,length(ecoliGN));
% genome=ecoliGN(1:20000);
genome=ecoliGN;
input=conv(genome-1,4.^(0:pw-1),'same');
a=length(input)-wd;
H=zeros(pw,a);

for pw=1:pwmax
    wd=4^pw+pw-1;
    input=conv(genome-1,4.^(0:pw-1),'same');
    a=length(input)-wd;
    for i=1:a; %% rolling window
        D=hist(input(i:i+wd-1),0:(4^pw));
    %   HD=hist(input(a:a+wd),1:4);
        ct=sum(D~=0);
        H(pw,i)=log(ct)/log(4^pw);
        fprintf('%d of %d \n',i,a);
    end
end
%%
imagesc(H)
%%
waterfall(H)
