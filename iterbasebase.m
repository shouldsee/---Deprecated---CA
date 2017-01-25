
function[cells]=iterbasebase(to,cells,cell2transi,async)
if to;
cells=torus(cells);
end

%% generate transition prob.
transi=cell2transi(cells);
% S_input=conv2((cells-0.5)*2,FIR,'same');


%% select cells for dicriminative update
    
uind=rand(size(cells))>async;
        
% if mean(async(:))==0;
%     uind=ones(size(cells));
%     elseif mean(async)==1;
%     uind=randi(1,numel(cells));
%     end
    pmat=rand(sum(uind(:)),1);
    cells(uind)=cells(uind).*(transi(uind)<pmat)+(1-cells(uind)).*(transi(uind)>pmat);
%     proj(S_input(uind)+1)>pmat;


end