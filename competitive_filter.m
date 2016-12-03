W=[0 -1 0 0;
   1 0 0 -1;
   0 0 0 -1;
   0 -1 1 0];
u=rand(4,1);
uall=[];
for epoch=1:10
u=(logsig(W*u));
uall=[uall,u];
end
imagesc(uall)