function [xn,zn] = OffG_Position_Refinement(x,z,K,L,phi,theta,beta,alpha,lambda,n,I,d_min)

Bx=[];
Bz=[];
H_re_t=[];
if n==1
    H_re_t=[];
else
    for t=1:n-1
        H_re_t=[H_re_t position_manifold(x(t),z(t),K,L,phi,theta,beta,lambda)];
    end
end
for iter=1:I
    %         if x(n)==0
    %             x(n)=1e-5;
    %         end
    %         if z(n)==0
    %             z(n)=1e-5;
    %         end
    H_re=[H_re_t position_manifold(x(n),z(n),K,L,phi,theta,beta,lambda)];
    F=inv(H_re'*H_re+alpha*eye(n))*H_re';
 
    F_T=F.';
    %         bx=kron(F_T(:,n),x(n)*par_b_x(x(n),z(n),K,L,phi,theta,beta,lambda));
    %          bz=kron(F_T(:,n),z(n)*par_b_z(x(n),z(n),K,L,phi,theta,beta,lambda));
    bx=kron(F_T(:,n),par_b_x(x(n),z(n),K,L,phi,theta,beta,lambda));
    bz=kron(F_T(:,n),par_b_z(x(n),z(n),K,L,phi,theta,beta,lambda));
    
    R=eye(K)-H_re*F;
    r=R(:);
    if n==1
        gamma_opt=0;
    else
        if Spacing_Constraint(0,bx,bz,r,x,z,n)>=lambda/2
            gamma_opt=0;
        else
            gamma_opt=binarySearch(@Spacing_Constraint,d_min,bx,bz,r,x,z,n);
        end
    end
    % gamma_opt=0;
    eta=bx'*r/(bx'*bx+gamma_opt);
    xi=bz'*r/(bz'*bz+gamma_opt);
    %         x(n) =x(n)+x(n)*real(eta);
    %         z(n) =z(n)+z(n)*real(xi);
    x(n) =x(n)+real(eta);
    z(n) =z(n)+real(xi);
end
xn=x(n);
zn=z(n);


end

