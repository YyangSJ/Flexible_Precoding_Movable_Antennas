function bz = par_b_z(x,z,K,L,phi,theta,beta,lambda)
bz=zeros(K,1);
for k=1:K
    for l=1:L
       bz(k)=bz(k)+sqrt(1/L)*conj(beta(l,k))*(-1i)*2*pi/lambda*theta(l,k)*exp(-1i*2*pi/lambda*(phi(l,k)*x+theta(l,k)*z));
    end 
end
end

