function [Phi_Phi, Phi_F, Phi_R, eyeNc] = mpc_model(num_c, den_c, Delta_t, Nc, Np, display)
[Ap,Bp,Cp,Dp]=tf2ss(num_c,den_c);
Delta_t=0.01;
[Ad,Bd,Cd,Dd]=c2dm(Ap,Bp,Cp,Dp,Delta_t);

[m1,n1]=size(Cd);
[n1,n_in]=size(Bd);
A_e=eye(n1+m1,n1+m1);
A_e(1:n1,1:n1)=Ad;
A_e(n1+1:n1+m1,1:n1)=Cd*Ad;
B_e=zeros(n1+m1,n_in);
B_e(1:n1,:)=Bd;
B_e(n1+1:n1+m1,:)=Cd*Bd;
C_e=zeros(m1,n1+m1);
C_e(:,n1+1:n1+m1)=eye(m1,m1);

n=n1+m1;
h(1,:)=C_e;
F(1,:)=C_e*A_e;
for kk=2:Np
    h(kk,:)=h(kk-1,:)*A_e;
    F(kk,:)= F(kk-1,:)*A_e;
end
v=h*B_e;
Phi=zeros(Np,Nc); % declare the dimension of Phi
Phi(:,1)=v; % first column of Phi
for i=2:Nc
    Phi(:,i)=[zeros(i-1,1);v(1:Np-i+1,1)]; %Toeplitz matrix
end
BarRs=ones(Np,1);
Phi_Phi= Phi'*Phi;
Phi_F= Phi'*F;
Phi_R=Phi'*BarRs;
eyeNc = eye(Nc,Nc);

if display=='YES'
[n,n_in]=size(B_e);
xm=zeros(n1,1);
Xf=zeros(n,1);
N_sim=100;
r=ones(N_sim,1);
u=0;
y=0;
for kk=1:N_sim;
DeltaU=inv(Phi_Phi+0.1*eye(Nc,Nc))*(Phi_R*r(kk)-Phi_F*Xf);
deltau=DeltaU(1,1);
u=u+deltau;
u1(kk)=u;
y1(kk)=y;
xm_old=xm;
xm=Ad.*xm+Bd.*u;
y=Cd*xm+Dd*u;
Xf=[xm-xm_old;y];
end

for i = 0:99
    k(i+1) = i * Delta_t;
end
figure
plot(k,y1)
grid;
end
end
