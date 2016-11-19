% %%
% rb=ruletemp{rind,1};
% rs=ruletemp{rind,2};
% rulecurr=[rb,rs]
preallocate

%%
% inputdict=squeeze(convn(ndict,FIR.S_input,'valid'));
% outdict=rulecurr(inputdict+1)';
% syms p v
% 
% pfunc=@(v)p^v*(1-p)^(9-v);
%%
mfunc=0;
I=size(outdict,1);
for i=1:I
    mfunc=mfunc+outdict(i)*pfunc(compdict(i));
end
%%

% mpfunc=matlabFunction(pfunc);
mmfunc=matlabFunction(mfunc);
%%
xs=linspace(0,1,1000);
% ys=arrayfun(@(p,v,out)mpfunc(m,p)*rulecurr(out+1),xs,compdict,outdict);

%%
plot(xs,mmfunc(xs));
refline([1 0])