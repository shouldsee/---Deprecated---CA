function fractal_growth
%����
%ÿ��0ϸ����1-��1-p)^n�ĸ��ʱ��2ϸ������p��һ����(�������0.245)��nΪ��Χ2ϸ��������0ϸ���0ϸ��
%���ο������еļ�����ģ�ͣ�
%2ϸ�����������1ϸ��
%1ϸ������������1ϸ��
a=zeros(200,200);
a(98:102,98:102)=2;%����2ϸ��䵱���ӵĽ�ɫ
k=[1 1 1
   1 0 1
   1 1 1];
p=1-(1-0.245).^(0:8);%Ԥ�����1-��1-p)^n
tic
while(1)
    a1=conv2(double(a==2),k,'same');
    a=double(a>=1)+... %��һ�֤(1,2)->1
      and(rand(200)<p(a1+1),a==0)*2; %��һ������0->2
    if true
  %     if toc>0.1
       tic
       image(a*24);axis equal
       drawnow
    end
end
