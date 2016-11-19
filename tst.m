%Example:
datanum = 1024;
outputnum = 4;
hiddennum = 256;
inputnum = 16;

inputdata = rand(datanum, inputnum,inputnum);
% outputdata=convn(inputdata,shiftdim(FIR.Sinput,-1),'valid');
% inputdata=reshape(inputdata,datanum,[]);
% outputdata=reshape(outputdata,datanum,[]);
% inputnum=size(inputdata,2);
% outputnum=size(outputdata,2);
inputdata=rand(datanum,inputnum);
outputdata = rand(datanum, outputnum);

dbn = randDBN([inputnum, outputnum]);
dbn = pretrainDBN( dbn, inputdata );
dbn = SetLinearMapping( dbn, inputdata, outputdata );
dbn = trainDBN( dbn, inputdata, outputdata );

% estimate = v2h( dbn, inputdata );
estout=v2h(dbn,inputdata);

imagesc(estout)
%%
imagesc(h2v(dbn,estout));

%%
imagesc(inputdata)
%%
% Example:
datanum = 1024;
outputnum = 16;
inputnum = 4;
inputdata = rand(datanum, outputnum);
dnn = randRBM( inputnum, outputnum );

% outputdata = h2v( dnn, inputdata );