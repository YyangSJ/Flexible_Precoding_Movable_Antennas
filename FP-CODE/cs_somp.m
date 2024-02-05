function [pos_array,hat_x]=cs_somp(y,T_Mat,G,L)
% y=T_Mat*x, T_Mat is n-by-m
% y - measurements
% T_Mat - combination of random matrix and sparse representation basis
% m - size of the original signal
% the sparsity is length(y)/4

[m,n]=size(y);
s=L; %  测量值维数

hat_x=zeros(G,n); %  待重构的谱域(变换域)向量                     
Aug_t=[];        %  增量矩阵(初始值为空矩阵)
r_n=y;  %  残差值 

for n=1:s; %  迭代次数(稀疏度是测量的1/4)
    pro=T_Mat'*r_n; 
    for i=1:G
       product(i)=sum(abs(pro(i,:)));
    end
   
    [val,pos]=max(product);   %最大投影系数对应的位置
        x_ind(n) = mod(pos - 1, sqrt(G)) + 1;
    z_ind(n) = ceil(pos / sqrt(G));
    Aug_t=[Aug_t,T_Mat(:,pos)];   %矩阵扩充
   % T_Mat(:,pos)=zeros(m,1); %选中的列置零
    aug_x=(Aug_t'*Aug_t+1e-7*eye(n))^(-1)*Aug_t'*y;   
  %   aug_x=pinv(Aug_t)*y;  
    r_n=y-Aug_t*aug_x;   %残差 
    pos_array(n)=pos;   %纪录最大投影系数的位置
    
end

hat_x(pos_array,:)=aug_x;  %  重构的向量 


end