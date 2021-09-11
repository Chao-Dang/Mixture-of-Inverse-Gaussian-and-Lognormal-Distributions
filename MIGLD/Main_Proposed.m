% This script aims to run the proposed methtod based on a  mixture distribuion with fractional
% moments for EVD reconstruction

clear all; close all;clc
%% Proposed method
 rng(0) % for reproducibility
num  = 25^2;
[ gg ] = Non_stationary_seismic_motions_lpss( num ); 
[ L1_625] = Newmark_linear_shear_frame_structure( num, gg );

%%

scale = 4;

load D1_0.mat
G=D1_0./scale;
GX=L1_625'./scale;


weights = ones(1,length(GX))./length(GX);

mu = mean(G);
st = std(G);

gmin = mu-4*st;
gmax = mu+10*st;
g = gmin:(gmax-gmin)/100:gmax;

k = [0.1,0.2,0.3,0.4,0.5]*4;

 [xx ,Pdf,Cdf ] = mixture_of_inverse_gaussian_and_lognormal_distribution( k,weights,GX,g);

 

 

PDF1 = xx(1).*pdf('InverseGaussian',g,xx(2),xx(3));
CDF1 = xx(1).*cdf('InverseGaussian',g,xx(2),xx(3));
PDF2 = (1-xx(1))*lognpdf(g,xx(4),xx(5));
CDF2 = (1-xx(1))*logncdf(g,xx(4),xx(5));


g = scale.*g;
G = scale.*G;




%%
figure(1)
[j,i]=hist(G,50);
j=j/length(G)/mean(diff(i));
b=bar(i,j,1);
hold on;
plot(g,Pdf./scale,'g--','LineWidth',2)
plot(g,PDF1./scale,'r:','LineWidth',2)
plot(g,PDF2./scale,'m-.','LineWidth',2)
%  plot(g,PDF,'m--','LineWidth',2)
xlim([0 70])
h=legend('MCS','Proposed method','Component A: inverse Gaussian','Component B: lognormal');
set(h,'Interpreter','latex','FontSize',12)
xlabel('$\rm Extreme~value, (mm)$','interpreter','latex','FontSize',12)
ylabel('$\rm PDF$','interpreter','latex','FontSize',12)
set(gca,'FontSize',12);
set(gca,'FontName','Timesnewroman');


figure(2)
gg = min(G):0.01:max(G);
h_mcs = hist(G,gg);
cdf_mcs = cumsum(h_mcs)/sum(h_mcs);
semilogy(gg,cdf_mcs,'b-','LineWidth',2)
hold on
semilogy(g,Cdf,'g--','LineWidth',2)
% semilogy(g,CDF1,'r:','LineWidth',2)
% semilogy(g,CDF2,'m-.','LineWidth',2)
% semilogy(g,1-CDF,'m--','LineWidth',2)
ylim([1e-6 1])
grid on
h=legend('MCS','Proposed method');
set(h,'Interpreter','latex','FontSize',12)
xlabel('$\rm Extreme~value, (mm)$','interpreter','latex','FontSize',12)
ylabel('$\rm CDF(log~scale)$','interpreter','latex','FontSize',12)
set(gca,'FontSize',12);
set(gca,'FontName','Timesnewroman');

figure(3)
gg = min(G):0.01:max(G);
h_mcs = hist(G,gg);
cdf_mcs = cumsum(h_mcs)/sum(h_mcs);
semilogy(gg,1-cdf_mcs,'b-','LineWidth',2)
hold on
semilogy(g,1-Cdf,'g--','LineWidth',2)
% semilogy(g,1-CDF1,'r:','LineWidth',2)
% semilogy(g,1-CDF2,'m-.','LineWidth',2)
% semilogy(g,1-CDF,'m--','LineWidth',2)

hold on 

% errorbar(g,mu_Poe,std_Poe)


ylim([1e-6 1])
grid on
h=legend('MCS','Proposed method');
set(h,'Interpreter','latex','FontSize',12)
xlabel('$\rm Extreme~value, (mm)$','interpreter','latex','FontSize',12)
ylabel('$\rm POE(log~scale)$','interpreter','latex','FontSize',12)
set(gca,'FontSize',12);
set(gca,'FontName','Timesnewroman');











