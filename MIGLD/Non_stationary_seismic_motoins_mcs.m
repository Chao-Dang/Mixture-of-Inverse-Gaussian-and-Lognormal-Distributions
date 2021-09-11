function [ gg ] = Non_stationary_seismic_motoins_mcs( num )
%NON_STATIONARY_SEISMIC_MOTOINS �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��


wg=5*pi;
sg=0.60;
wf=0.5*pi;
sf=0.60;

amax=200;                     %��ֵ���ٶ�
gama=2.8;
s_=amax^2/((gama^2)*(pi*wg*(2*sg+1/(2*sg))));            %��ǿ������


N=1000; 
i=1:N;
wu=100; %ԲƵ������ֵ
dw = 0.1;
wl=wu-N*dw;                         %ԲƵ������ֵ
% dw=(wu-wl)./N;
w=wl+i.*dw;                   %Ƶ�ʻ���
t1=0.5;                            %���纯������
t2=10;
c=0.45;
dt=0.02;
T=20;
t=0:dt:T;                       %ʱ�仮��

%%%% Clough-Penzien SDF %%%%
sw=2.*(wg.^4+(2.*sg.*wg.*w).^2).*(w.^4).*s_./(((w.^2-wg.^2).^2+(2.*sg.*wg.*w).^2).*((w.^2-wf.^2).^2+(2.*sf.*wf.*w).^2));

%%%%  Amplitude calculation %%%%
A=sqrt(2*pi.*sw.*dw./pi);

%%%% Random Phase Angle MCS��%%%%
% fai=unifrnd(0,2.*pi,[N,num]);           %%%%���ȷֲ����������

%%%% Envelope Function %%%%
a=(t./t1).^2.*(t<=t1)+1.*(t>t1&t<t2)+exp(-c.*(t-t2)).*(t>=t2);

%%%% random non-stationary  seismic acceleration (only amplitude is non-stationary) 

x=zeros(1001,num);
for j=1:num 
    j
    fai=unifrnd(0,2.*pi,[1,N]);           %%%%���ȷֲ����������
for i=1:N
    x(:,j)=x(:,j)+(A(i)*(cos(w(i).*t+fai(i))))';
end
    gg(:,j)=a'.*x(:,j);
end


end

