function H_SWAP = AS_SWAP(H,N,G,K)
for k=1:N
    for n=1:G
        array_sel(k)=n;
        H_FIX2=H(array_sel,:);
        ecap(n)=det(eye(K)+H_FIX2'*H_FIX2);
    end
    [~,indm]=max(ecap);
    array_sel(k)=indm;
end
H_SWAP=H(array_sel,:);
end

