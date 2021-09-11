function [ gg ] = Non_stationary_seismic_motoins_mcs( num )
%NON_STATIONARY_SEISMIC_MOTOINS 此处显示有关此函数的摘要
%   此处显示详细说明


wg=5*pi;
sg=0.60;
wf=0.5*pi;
sf=0.60;

amax=200;                     %峰值加速度
gama=2.8;
s_=amax^2/((gama^2)*(pi*wg*(2*sg+1/(2*sg))));            %谱强度因子


N=1000; 
i=1:N;
wu=100; %圆频率上限值
dw = 0.1;
wl=wu-N*dw;                         %圆频率下限值
% dw=(wu-wl)./N;
w=wl+i.*dw;                   %频率划分
t1=0.5;                            %包络函数参数
t2=10;
c=0.45;
dt=0.02;
T=20;
t=0:dt:T;                       %时间划分

%%%% Clough-Penzien SDF %%%%
sw=2.*(wg.^4+(2.*sg.*wg.*w).^2).*(w.^4).*s_./(((w.^2-wg.^2).^2+(2.*sg.*wg.*w).^2).*((w.^2-wf.^2).^2+(2.*sf.*wf.*w).^2));

%%%%  Amplitude calculation %%%%
A=sqrt(2*pi.*sw.*dw./pi);

%%%% Random Phase Angle MCS　%%%%
% fai=unifrnd(0,2.*pi,[N,num]);           %%%%均匀分布的随机变量

%%%% Envelope Function %%%%
a=(t./t1).^2.*(t<=t1)+1.*(t>t1&t<t2)+exp(-c.*(t-t2)).*(t>=t2);

%%%% random non-stationary  seismic acceleration (only amplitude is non-stationary) 

x=zeros(1001,num);
for j=1:num 
    j
    fai=unifrnd(0,2.*pi,[1,N]);           %%%%均匀分布的随机变量
for i=1:N
    x(:,j)=x(:,j)+(A(i)*(cos(w(i).*t+fai(i))))';
end
    gg(:,j)=a'.*x(:,j);
end


end

