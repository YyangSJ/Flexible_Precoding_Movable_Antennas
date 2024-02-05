function g_s = PW(theta,phi,d,lambda,N,M)
for m=1:M
    for n=1:N 

         r_s(m,n)=-((n-1)*d)*theta-((m-1)*d)*phi;   
        G_s(m,n)=exp(-1i*2*pi/lambda*(r_s(m,n))); 
    end
end
g_s=G_s(:);
end

