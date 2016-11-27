fa=5;fb=5;
a=10;
b=10;
t=0;
% t=1;
h=1.0E-2;
as=linspace(0,10,1/h);
bs=linspace(0,10,1/h);
[x,y]=ndgrid(as,bs);
z=(sin(fa*pi*x/a).*sin(fb*pi*y/b).*exp(i*t));
subplot(2,1,1)
sf1=surf(x,y,real(z));
zlim([-1 1])
caxis([-1 1])
subplot(2,1,2);
sf2=surf(x,y,imag(z));
zlim([-1 1])
caxis([-1 1])
tstep=1E-1;
while true
    z=(sin(fa*pi*x/a).*sin(fb*pi*y/b).*exp(i*t));
    set(sf1,'CData',real(z),'ZData',real(z));
    set(sf2,'CData',imag(z),'ZData',imag(z));
    t=t+tstep;
    pause(0.02)
end

% shading interp