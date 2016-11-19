function [ rle ] = mat2rle( mat )
% tic
[m,n]=size(mat);
old=-1;
% ones(1,n)*-1
rle='';

dict='$bo';
count=1;
% rle=repmat(' ',1,m*n*2);
pos=1;
for i=1:m
%     row=mat(i,:);
    for j=1:n
        
        curr=mat(i,j);
        if curr==old;
            count=count+1;
        else
            rle=[rle num2str(count) dict(old+2)];
%             num=num2str(count);
%             rle(pos:pos)=;
%             rle(pos+1)=dict(old+2);
%             pos=pos+2;

            count=1;
        end
        old=curr;

    end
  rle=[rle num2str(count) dict(old+2)];
%     rle(pos)=num2str(count);
%     rle(pos+1)=dict(old+2);
%     pos=pos+2;
    count=1;
    old=-1;
end
rle(pos)='!';
% f
% toc

% rle='oboboboboob$bbbo!';

%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here


        end


