function y = com_AL(g,n)
%% Compute  left  fractional differential 
% n ����������� �ڵ����Ϊn-1
AL=zeros(n-1,n-1);
for i=1:n-1
    for j=1:n-1
        if i==j
            AL(i,j)=g(2);%���Խ���Ԫ��
        elseif j==i-1
            AL(i,j)=g(3);%�´ζԽ���Ԫ��
        elseif i==j-1
            AL(i,j)=g(1);%�ϴζԽ���Ԫ��
        elseif j<i-1
            AL(i,j)=g(i-j+2);%���·�Ԫ��
        else
            AL(i,j)=0;%�����ط�Ϊ��
        end
    end
end
y = AL;
end
