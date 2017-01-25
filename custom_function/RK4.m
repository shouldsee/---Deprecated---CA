function[incr]=RK4(eof,dof,dt);
wt=[1 2 2 1]/6;
k1=eof(dof);
k2=eof(dof+dt/2*k1);
k3=eof(dof+dt/2*k2);
k4=eof(dof+dt*k3);
incr=dt*sum(bsxfun(@times,wt,[k1 k2 k3 k4]),2);

end