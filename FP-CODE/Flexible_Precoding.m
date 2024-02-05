function [F,H,x,z]=Flexible_Precoding(y,T_Mat,G,N,alpha,phi,theta,beta,K,L,lambda,I,d,d_min)
 

[m,n]=size(y);
hat_x=zeros(G,n); %  待重构的谱域(变换域)向量
Aug_t=[];        %  增量矩阵(初始值为空矩阵)
r_n=y;  %  残差值
x=[];z=[];
MG=1:G;
for n=1:N
    pro=T_Mat'*r_n;
    product=[];
    for g=1:length(MG)
        product(g)=sum(abs(pro(MG(g),:)));
    end
    
    [val,pos]=max(product);   %最大内积系数对应的位置
    g_opt=MG(pos);
    x_ind(n) = mod(g_opt(1) - 1, sqrt(G)) + 1;
    z_ind(n) = ceil(g_opt(1) / sqrt(G));
    x(n)=(x_ind(n)-1)*d;
    z(n)=(z_ind(n)-1)*d;
    [x(n),z(n)]=OffG_Position_Refinement(x,z,K,L,phi,theta,beta,alpha,lambda,n,I,d_min);
    MG_temp=MG;
    g_cel=[];
    x_setind=[];z_setind=[];
    for gg=1:length(MG_temp)
        x_setind(gg) = mod(MG_temp(gg) - 1, sqrt(G)) + 1;
        z_setind(gg) = ceil(MG_temp(gg) / sqrt(G));
        if sqrt((x(n)-(x_setind(gg)-1)*d)^2+(z(n)-(z_setind(gg)-1)*d)^2)<d_min
            g_cel=[g_cel MG_temp(gg)];
        end
    end
    MG=MG(~ismember(MG, g_cel));
    Aug_t=[Aug_t,position_manifold(x(n),z(n),K,L,phi,theta,beta,lambda)];   %矩阵扩充
    T_Mat(:,g_opt)=zeros(m,1); %选中的列置零
    aug_x=(Aug_t'*Aug_t+alpha*eye(n))^(-1)*Aug_t'*y;
    r_n=y-Aug_t*aug_x;   %残差
    
end

H=Aug_t;
F=aug_x;
end