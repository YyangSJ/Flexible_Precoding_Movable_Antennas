function fun_g = Spacing_Constraint(gamma,bx,bz,r,x,z,n)
eta=bx'*r/(bx'*bx+gamma);
xi=bz'*r/(bz'*bz+gamma);
    x(n) =x(n)+x(n)*real(eta);
    z(n) =z(n)+z(n)*real(xi); 
for t=1:n-1
    spa(t)=sqrt((x(n)-x(t))^2+(z(n)-z(t))^2);
end
fun_g=min(spa);
end

