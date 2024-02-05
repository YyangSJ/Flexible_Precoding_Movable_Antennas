function [Lambda,Gain] = OLS_MMV(Y,PHI,L)
% Distributed orthogonal least squares for different sensing matrices
% sampling the common subspace.Written by Songjie Yang-9.13.2023.

[P,M]=size(Y);
G=size(PHI,2);
Lambda=[];
for l=1:L   
    for m=1:M
        for g=1:G
%             if ismember(g,Lambda)
%                 g=g+1;
%             end 
            PHI_G=PHI(:,[Lambda,g]);
            Res(g,m)=norm((eye(P)-PHI_G*((PHI_G'*PHI_G)\PHI_G'))*Y(:,m))^2;
        end 
    end 
    Res_com=sum(Res,2);
    [~,indm]=min(Res_com);
    Lambda=[Lambda,indm];
end

    PHI_G=PHI(:,Lambda);
for m=1:M

    Gain(:,m)=((PHI_G'*PHI_G)\PHI_G')*Y(:,m);
end

end
