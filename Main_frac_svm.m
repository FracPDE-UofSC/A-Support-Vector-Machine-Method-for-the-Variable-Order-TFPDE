%% using support vector machine method for solving space fractional equations
clc;clear;
n = 2^10;
lam = 1e7;
sig = 1e-1;
alpha = 1.9;
h = 1/n;
dim = n+1;
theta = 0.9; 
D = 1; 
for i = 1:n+1
    x(i) = 0 + (i-1)*h;
end
% K = 1/(2*cos(pi*alpha/2));
K = 1;
hh = (1/h)^alpha;
fuv = @(u,v) exp(-(u-v)^2/sig^2);%kernel function
z1=@(xx) -10/(2*(cos(alpha*pi/2)));
z2=@(xx) (2/gamma(3-alpha))*(xx.^(2-alpha)+(1-xx).^(2-alpha));
z3=@(xx) (12/gamma(4-alpha))*(xx.^(3-alpha)+(1-xx).^(3-alpha));
z4=@(xx) (24/gamma(5-alpha))*(xx.^(4-alpha)+(1-xx).^(4-alpha));
fr = @(xx) z1(xx)*(z2(xx)-z3(xx)+z4(xx));
fr =  @(xx) - D * 10 * 2*( theta *  (xx.^(2-alpha)/gamma(3-alpha)- 6*xx.^(3-alpha)/gamma(4-alpha) + 12*xx.^(4-alpha)/gamma(5-alpha) ) ...
    +  (1-theta) * ( (1-xx).^(2-alpha)/gamma(3-alpha) - 6*(1-xx).^(3-alpha)/gamma(4-alpha)   ...
    + 12*(1-xx).^(4-alpha)/gamma(5-alpha) ) );
I = eye(n-1); 
f = zeros(n+2,1);
g = g_alpha(n,alpha);
%% Compute  left  fractional differential 
AL = com_AL(g,n);%  dim(AL) = n-1;
%%
AR = AL';% Compute  left  fractional differential 
fuv1 = @(u) exp(-(u-x(1))^2/sig^2);  
fuvn = @(u) exp(-(u-x(n+1))^2/sig^2); 
for i = 2:n 
    Fuv1(i-1) = fuv1(x(i));
    Fuvn(i-1) = fuvn(x(i));
end
for i = 2:n
    for j = 2:n
        Fuv(i-1,j-1) =  fuv(x(i),x(j));
    end
end
k11 = hh*((theta*AL+(1-theta)*AR)*Fuv);
k11 = hh*(theta*AL+(1-theta)*AR)*k11';
k12 = hh*(theta*AL+(1-theta)*AR)*Fuv1';
k13 = hh*(theta*AL+(1-theta)*AR)*Fuvn';

k14 = ones(n-1,1);
k14 = hh*(theta*AL+(1-theta)*AR)*k14;
for i = 2:n
    f(i-1) = fr(x(i));
end

k11 =  K^2*k11;
k12 =  K*k12;
k13 =  K*k13;
k14 = -K*k14;
k21=   k12';
k31 =  k13';
k11 = k11 + I./lam;
k22 = K*fuv(x(1),x(1));
k23 = K*fuv(x(n+1),x(1));
k32 = K*fuv(x(1),x(n+1));
k33 = K*fuv(x(n+1),x(n+1));
k41 = k14';
A = [k11,k12,k13,k14;
    k21,k22,k23,-1;
    k31,k32,k33,-1;
    k41,-1,-1,0];
ha = A\f;
beta0 = ha(1:n-1);
beta1 = ha(n);%
beta2 = ha(n+1);
b = ha(n+2);
[u_svm,u_exact]= deal(zeros(n+1,1));
kk = zeros(n-1,1);
%% svm_solution u_svm
for i = 1:n+1
    for j = 1:n-1
        fuj(j) = fuv(x(j+1),x(i));
    end
    kk = hh*((theta*AL+(1-theta)*AR))*fuj';
    u_svm(i) = b - K*kk'*beta0 - beta1*fuv(x(1),x(i)) - beta2*fuv(x(n+1),x(i));
end
%% True_solution  u_exact
for i = 1:n+1
    u_exact(i) = 10*x(i)^2*(1-x(i))^2; % frac_true solution
end

%% Error
table(u_exact,u_svm);
error = norm(u_exact - u_svm,inf); 

%% Plot 
plot(x,u_exact,'r')
hold on;
plot(x,u_svm,'b--')
legend('Uexact','Usvm')

table(dim,alpha,lam,sig,error)