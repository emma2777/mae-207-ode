function NR_Burgers_IMEXRKCB3e_1
% Simulate the 1D Burgers on 0<x<L with homogeneous Dirichlet BCs using IMEXRKCB3e in time
% (explicit on nonlinear terms, implicit on linear terms) & 2nd-order central FD in space.

%%%%%%%%%%%%%%%%%%%% Initialize the simulation paramters (user input) %%%%%%%%%%%%%%%%%%%%
L=100; Tmax=50; N=100; dt=0.5; PlotInterval=1;
dx=L/N; x=(0:N)*dx; y=-sin(pi*x/L)-sin(2*pi*x/L)+sin(6*pi*x/L); NR_PlotXY(x,y,0,0,L,-3,3)
%%%%%%%%%%%% Precalculate the time-stepping coefficients used in the simulation %%%%%%%%%%
c = [0,1/3,1,1];
b = [0,3/4,-1/4,1/2];
a_im = [0,0,0,0;
        0,1/3,0,0;
        0,1/2,1/2,0;
        0,3/4,-1/4,1/2];
a_ex = [0,0,0,0;
        1/3,0,0,0;
        0,1,0,0;
        0,3/4,1/4,0];
A = diag(-2*ones(1,N-1)) + diag(1*ones(1,N-2),1) + diag(1*ones(1,N-2),-1);
for k=1:Tmax/dt
  for rk=1:4 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ALL RK SUBSTEPS %%%%%%%%%%%%%%%%%%%%%%%%%%%% 
      if (rk==1)
              r = y;
        else
              r(2:N) = y(2:N)+(a_im(rk,rk-1)-b(rk-1))*dt*z+(a_ex(rk,rk-1)-b(rk-1))*dt*r(2:N);
      end
      z = ((eye(size(A))-a_im(rk,rk)*dt*A)\A*r(2:N)')';
        r(2:N) = r(2:N)+a_im(rk,rk)*dt*z;
        r(2:N) = -r(2:N).*(r(3:N+1)-r(1:N-1))/2/dx;
        y(2:N) = y(2:N)+b(rk)*dt*z+b(rk)*dt*r(2:N);

  end %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% END OF RK LOOP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if (mod(k,PlotInterval)==0) NR_PlotXY(x,y,k*dt,0,L,-3,3); end            
end                                                                     
end