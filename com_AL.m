function y = com_AL(g,n)
%% Compute  left  fractional differential 
% n 代表区间个数 内点个数为n-1
AL=zeros(n-1,n-1);
for i=1:n-1
    for j=1:n-1
        if i==j
            AL(i,j)=g(2);%主对角线元素
        elseif j==i-1
            AL(i,j)=g(3);%下次对角线元素
        elseif i==j-1
            AL(i,j)=g(1);%上次对角线元素
        elseif j<i-1
            AL(i,j)=g(i-j+2);%左下方元素
        else
            AL(i,j)=0;%其他地方为零
        end
    end
end
y = AL;
end
