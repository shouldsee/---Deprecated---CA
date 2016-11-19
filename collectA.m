% lists={'list1','list2','list3','list4','list5'};
% scalars={'St1_Ent','St2_Ent','H_input'};
% outnumber=8;
% scalars={'phaseX','phaseY','phaseZ','H_input','HD_linput'};
% outnumber=7;
% scalars={'H_pop_lagged','MI','H_input'};
% scalars={'NFlux','RFlux','ANFlux','H_input'};
% scalars={'DD','V_input','H_D','L_D','H_input'};
outnumber=13;
scalars={'phaseX','phaseY','phaseZ','H_input'};

% fname='a';
%%
storeM=1;
rrule=1;
store=1;
import=0;
detail=0;
fname='pcami_input_A'; %
v=2;
k=8;


haulsize=50; %increse haulsize to achieve less frequent I/O to speed up.
haulnumber=0;


CA_console_opt2;