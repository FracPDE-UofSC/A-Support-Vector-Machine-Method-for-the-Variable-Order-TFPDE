function g = g_alpha(n,alpha)
% n == dim(g)

g=zeros(1,n-1);%”√”⁄¥Ê¥¢ g(alpha)
g(1)=-alpha;
for i=2:n-1
    g(i)=(1-(alpha+1)/i)*g(i-1);
end
g=[1,g];

end