function y = com_AL(g,n)
%% Compute  left  fractional differential
AL=zeros(n-1,n-1);
for i=1:n-1
    for j=1:n-1
        if i==j
            AL(i,j)=g(2);
        elseif j==i-1
            AL(i,j)=g(3);
        elseif i==j-1
            AL(i,j)=g(1);
        elseif j<i-1
            AL(i,j)=g(i-j+2);
        else
            AL(i,j)=0;
        end
    end
end
y = AL;
end
