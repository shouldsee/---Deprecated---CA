
function[cells]=update(cells)
global rulecurr async n a b
cells=torus(cells);
S_input=conv2(cells,[1 1 1 ; 1 9 1 ; 1 1 1],'same');
if async==0;
    cells=rulecurr(S_input+1);
elseif async==1;
    uind=rand(size(cells))>ndgrid(linspace(-0.5,1.5,a),linspace(-0.5,1.5,b));
    cells(uind)=rulecurr(S_input(uind)+1);

else
    uind=rand(size(cells))>async;
    cells(uind)=rulecurr(S_input(uind)+1);
%     cells=1-cells;
end

end