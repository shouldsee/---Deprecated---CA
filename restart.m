
%% modify rule list
% rulename=repmat({'B3/S23'},10,1);
% ruletemp=importrule(rulename);


%%
rrule=0;
import_rule=1;
haulsize=1000;
storeM=1;
store=1;
fname='man_test1';

%%
clear welist
rulenumber=0;
soupnumber=0;
start=1;
k=8;
v=2;
rrule=0;
import_rule=1;
detail=1;
soupmax=6;
% lists={'list1','list2','list3','list4','list5'};
% outnumber=13;
% % scalars={'E_dmat_V_input','H_input','S_pop'};
% scalars={'S_pop','H_input'};
% % scalars={'S_pop','H_input','D_input'};
% % loadrule
% % outnumber=8;
% % scalars={'NFlux','RFlux','ANFlux','H_input'};
% scalars={'phaseX','phaseY','phaseZ','H_input','S_pop','S_age'};
% scalars={'phaseX','phaseY','H_input','S_pop'};
% outnumber=8;
% scalars={'phaseX','phaseY','phaseZ','H_input','HD_linput'};
% outnumber=7;
% scalars={'H_pop_lagged','MI','H_input'};
% outnumber=5;
% scalars={'LS_pop','phaseX','phaseY','H_input'};
outnumber=7;
scalars={'LS_pop','phaseX','phaseY','phaseZ','H_input'};


% CA_console_opt
% CA_mc
CA_console_opt2
