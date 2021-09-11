function [ D1_0] = Newmark_linear_shear_frame_structure( num, g )
% NEWMARK_LINEAR_SHEAR_FRAME_STRUCTURE 此处显示有关此函数的摘要
%   此处显示详细说明
mg = g;
storey=12;
m=[3.487	3.852	3.225	2.887	2.667	2.558	2.558	2.558	2.558	2.558	2.558	2.558]'.*1e5;
k=[1.618     1.562  1.454    1.258    0.975   0.882	  0.882	 0.882	 0.882	  0.882   0.882   0.882]'.*1e8;
% cov(1:storey,1)=0.15;
% [m_mu,m_sig]=lognormal_permaters(mm,cov);
% [k_mu,k_sig]=lognormal_permaters(kk,cov);
% for aa=1:storey
% m(aa,:)=lognrnd(m_mu(aa),m_sig(aa),[1,num]);
% k(aa,:)=lognrnd(k_mu(aa),k_sig(aa),[1,num]);
% end

    %形成质量矩阵
mass=zeros(storey,storey);
for i=1:storey
   mass(i,i)=m(i);
end

%形成刚度矩阵
stiffness=zeros(storey,storey);
for i=1:storey-1
  stiffness(i,i)=k(i)+k(i+1);
  stiffness(i,i+1)=-k(i+1);
  stiffness(i+1,i)=-k(i+1);
end
stiffness(storey,storey)=k(storey);

dis(:,1)=zeros(storey,1);
vel(:,1)=zeros(storey,1);
acc(:,1)=zeros(storey,1);
dt=0.02;
t=0:dt:20;
%%%读取地震   



%%%%%
gama=1/2;
beta=1/6;
% damp= [2086072.06000000,-681600,0;-681600,1400137.54000000,-631900;0,-631900,712670.080000000];

for iii=1:num
    iii
    


%求解结构动力特性
 [x,d]=eig(stiffness,mass);
d=diag(sqrt(d));
for i=1:storey
[d1(i),j]=min(d);
xgd(:,i)=x(:,j);
d(j)=max(d)+1;
end
w=d1;
x=xgd;

%求解阻尼系数
a1=2*w(1)*w(2)*(0.05*w(2)-0.05*w(1))/(w(2)^2-w(1)^2);
a2=2*(0.05*w(2)-0.05*w(1))/(w(2)^2-w(1)^2);

%形成阻尼矩阵
damp=a1.*mass+a2.*stiffness;
    
    
    
stiffness_=stiffness+gama/(beta*dt).*damp+1/(beta*dt^2).*mass;
a=1/(beta*dt).*mass+gama/beta.*damp;
b=1/(2*beta).*mass+dt*(gama/(2*beta)-1).*damp;
for ii=2:1001
    ff_(:,ii)=-mass*ones(storey,1).*mg(ii,iii)+mass*ones(storey,1).*mg(ii-1,iii)+a*vel(:,ii-1)+b*acc(:,ii-1);
    deltadis(:,ii)= inv(stiffness_)*ff_(:,ii);
    deltavel(:,ii)=gama/(beta*dt).*deltadis(:,ii)-gama/(beta).*vel(:,ii-1)+dt*(1-gama/(2*beta)).*acc(:,ii-1);
    deltaacc(:,ii)=1/(beta*dt^2).*deltadis(:,ii)-1/(beta*dt).*vel(:,ii-1)-1/(2*beta).*acc(:,ii-1);
    dis(:,ii)=dis(:,ii-1)+deltadis(:,ii);
    vel(:,ii)=vel(:,ii-1)+deltavel(:,ii);
    acc(:,ii)=acc(:,ii-1)+deltaacc(:,ii);
end
  D1_0(iii)=10*max(abs(dis(1,:)));
  D2_1(iii)=10*max(abs(dis(2,:)-dis(1,:)));
  D3_2(iii)=10*max(abs(dis(3,:)-dis(2,:)));
  D4_3(iii)=10*max(abs(dis(4,:)-dis(3,:)));
  D5_4(iii)=10*max(abs(dis(5,:)-dis(4,:)));
  D6_5(iii)=10*max(abs(dis(6,:)-dis(5,:)));
  D7_6(iii)=10*max(abs(dis(7,:)-dis(6,:)));
  D8_7(iii)=10*max(abs(dis(8,:)-dis(7,:)));
  D9_8(iii)=10*max(abs(dis(9,:)-dis(8,:)));
  D10_9(iii)=10*max(abs(dis(10,:)-dis(9,:)));
  D11_10(iii)=10*max(abs(dis(11,:)-dis(10,:)));
  D12_11(iii)=10*max(abs(dis(12,:)-dis(11,:)));
end


end

