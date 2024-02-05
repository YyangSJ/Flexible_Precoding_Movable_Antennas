function [H_sel] = AS_SEL(H,N,G,K,P_max,sigma_2)
for n=1:N
    ini_a(n)=H(n,:)*H(n,:)';
end
[~,indm_ini]=max(abs(ini_a));
H_sel=H(indm_ini,:);
for n=1:N-1 
    for g=1:G
        gain_a(g)=H(g,:)*inv(eye(K)+P_max/K/sigma_2*H_sel'*H_sel)*H(g,:)';
    end
    [~,indm_iter]=max(abs(gain_a));
    H_sel=[H_sel ;H(indm_iter,:)];
end
%F_SEL= H_sel/(H_sel'*H_sel+alpha*eye(K));
%H_SEL=H_sel';
end

