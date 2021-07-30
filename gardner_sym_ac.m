options = optimset('Display','off');
options.MaxIter = 5000000 ;
options.MaxFunEvals = 5000000 ;
options.FunctionTolerance = 1e-16;
options.OptimalityTolerance = 1e-16;
options.StepTolerance = 1e-16;

k = 0;
iter = 0;
iter = cast(iter,'uint8');
k_max = 2;
delta_k = 0.001;
size = 2000;
alpha_c = zeros(1,size);
alpha_cc = zeros(1,size);

kk = zeros(1,size);

for i = 1:size

    F = @(x) (k-x)*0.5*(1-erf((x-k)/sqrt(2))) + (1/sqrt(2*pi))*exp(-0.5*(x-k)*(x-k)) - x; 
    r = fzero(F,0.1);

    %I = (k-r)*exp(-0.5*(r-k)*(r-k))/sqrt(2*pi) + (0.5*(k-r)*(k-r) + 1)*(1 - erf((r-k)/sqrt(2))) + r*r;
    I = 0.5 + (k-r)*0.398942*exp(-0.5*(k-r)*(k-r)) + 0.5*(k-r)*(k-r) + (0.5 + 0.5*(k-r)*(k-r))*erf(0.707107*(k-r)) + r*r;
    alpha_c(i) = 0.5/I;
    alpha_cc(i) = 1/(0.5 + k*0.398942*exp(-0.5*k*k) + 0.5*k*k + (0.5 + 0.5*k*k)*erf(0.707107*k));
    k = i*delta_k;
    kk(i) = k;
end

xlabel('k');
ylabel('\alpha_c');
hold on
plot(kk,alpha_c)
%plot(kk,0.4*ones(1,size))
%plot(kk,0.3*ones(1,size))
%plot(kk,0.5*ones(1,size))
plot(kk,0.65*ones(1,size))