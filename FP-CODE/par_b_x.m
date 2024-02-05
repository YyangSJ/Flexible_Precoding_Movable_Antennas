function bx = par_b_x(x,z,K,L,phi,theta,beta,lambda)
bx=zeros(K,1);
for k=1:K
    for l=1:L
       bx(k)=bx(k)+sqrt(1/L)*conj(beta(l,k))*(-1i)*2*pi/lambda*phi(l,k)*exp(-1i*2*pi/lambda*(phi(l,k)*x+theta(l,k)*z));
    end 
end
end

