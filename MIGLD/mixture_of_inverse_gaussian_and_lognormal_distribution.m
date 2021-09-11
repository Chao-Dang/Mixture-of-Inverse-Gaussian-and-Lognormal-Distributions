function [xx ,Pdf,Cdf ] = mixture_of_inverse_gaussian_and_lognormal_distribution( k,weights,GX,g)
%MIXTURE_OF_TWO_LOGNORMAL 此处显示有关此函数的摘要
%   此处显示详细说明

% initial value 
mu = weights*GX;
va = weights*GX.^2-mu.^2;

para1 = mu;
para2 = mu^3/va;

para3 = log((mu.^2)./sqrt(va+mu.^2));
para4 = sqrt(log(va./(mu.^2)+1));

x0 = [0.5,para1,para2,para3,para4];


fun = @(x)[ analytical_raw_moments(x,k(1))-estimated_raw_moments(weights,GX,k(1));
    analytical_raw_moments(x,k(2))-estimated_raw_moments(weights,GX,k(2));
    analytical_raw_moments(x,k(3))-estimated_raw_moments(weights,GX,k(3));
    analytical_raw_moments(x,k(4))-estimated_raw_moments(weights,GX,k(4));
    analytical_raw_moments(x,k(5))-estimated_raw_moments(weights,GX,k(5));
];
AlGO = {'trust-region-dogleg','trust-region','levenberg-marquardt'};
opts = optimoptions('fsolve','Algorithm',AlGO{3},'MaxFunctionEvaluations',1e4,'FunctionTolerance',1e-100,'StepTolerance',1e-100,'MaxIterations',1e4,'Display','iter');
xx = fsolve(fun,x0,opts);

Pdf = xx(1).*pdf('InverseGaussian',g,xx(2),xx(3))+ (1-xx(1)).*lognpdf(g,xx(4),xx(5));
Cdf = xx(1).*cdf('InverseGaussian',g,xx(2),xx(3))+ (1-xx(1)).*logncdf(g,xx(4),xx(5));


end

function [raw_a] = analytical_raw_moments(x,k)

raw_a = x(1)*exp(x(3)/x(2))*sqrt(2*x(3)/pi)*x(2)^(k-1/2)*besselk(1/2-k,x(3)/x(2))+(1-x(1))*exp(k*x(4)+k^2*x(5)^2/2);

end

function [raw_e] = estimated_raw_moments(weights,GX,k)

raw_e = weights*GX.^(k);



end
