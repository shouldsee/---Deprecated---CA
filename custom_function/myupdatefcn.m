function txt = myupdatefcn(~,event_obj,r,id)
% Customizes text of data tips
pos = get(event_obj,'Position');
try
%     X=num2str(pos(1));
    Y=pos(2);

    I = get(event_obj, 'DataIndex');
txt = {['X: ',num2str(pos(1))],...
       ['Y: ',num2str(pos(2))],...
       ['I: ',num2str(I)],...
        ['Rule: ',r{eval(id)}]};
    %        ['Rule:',r{eval(I)}]};

    %        ['Z: ',num2str(pos(3))],...

catch err
    txt={num2str(pos(1))};
    txt=num2str(eval('I'))
end
%        ['T: ',num2str(t(I))]};

