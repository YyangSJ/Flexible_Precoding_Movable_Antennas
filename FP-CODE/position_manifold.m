function b = position_manifold(x,z,K,L,phi,theta,beta,lambda)
b=zeros(K,1);
for k=1:K
    for l=1:L
       b(k)=b(k)+sqrt(1/L)*conj(beta(l,k))* exp(-1i*2*pi/lambda*(phi(l,k)*x+theta(l,k)*z));
    end 
end
end

