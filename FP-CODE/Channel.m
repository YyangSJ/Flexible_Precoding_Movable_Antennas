function h = Channel(beta,phi_r,theta_r,phi_t,theta_t,lambda,d,N,L,x,z)
h=zeros(N,1);
for l=1:L
    h=h+sqrt(1/L)*beta(l)*exp(1i*2*pi/lambda*(x*phi_r(l)+z*theta_r(l)))*PW(theta_t(l),phi_t(l),d,lambda,sqrt(N),sqrt(N));
end

end

