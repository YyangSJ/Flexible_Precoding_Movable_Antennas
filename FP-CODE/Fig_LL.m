clear all;
close all;
rng(0)
f=3e9; % system frequency
lambda=3e8/f; % antenna wavelength
d=lambda/2; % inter-element spacing for virtual channel representation
d_min=lambda/2; % allowed minimal inter-element spacing in practice
K=4; % number of users
N=4; % number of movable antennas
G_16=16; % number of virtual antennas 16 
G_64=64; % number of virtual antennas 64 
LL=1:1:15; % number of channel paths
sigma_2=1; % noise power
P_max=1; % transmit power
alpha=1; % regularized factor
I=100; % iteration number for OffG procedure
realization=33;
for reali=1:realization
    for ll=1:length(LL)
        L=LL(ll);
        beta=normrnd(0,1,L,K)+1i*normrnd(0,1,L,K); % path gains
        phi_t=unifrnd(-1,1,L,K); % BS path DoA-azimuth
        theta_t=unifrnd(-1,1,L,K); % BS path DoA-elevation
        phi_r=zeros(L,K);  theta_r=zeros(L,K); % user antenna is fixed, the angles are set to 0.
        for k=1:K
            H_16(:,k)=Channel(beta(:,k),phi_r(:,k),theta_r(:,k),phi_t(:,k),theta_t(:,k),lambda,d,G_16,L,0,0);
            H_64(:,k)=Channel(beta(:,k),phi_r(:,k),theta_r(:,k),phi_t(:,k),theta_t(:,k),lambda,d,G_64,L,0,0);
        end
        
        %% Fixed Antenna Position for Precoding
        array_sel=1:N;
        H_FIX=H_64(array_sel,:);
        F_FIX= H_FIX/(H_FIX'*H_FIX+alpha*eye(K));
        pow_FIX=trace(F_FIX'*F_FIX);
        H_FIX=H_FIX'; 
        
        for k=1:K
            F_FIX(:,k) = F_FIX(:,k) / norm(F_FIX(:,k)); % normalization
        end
        %% Flexible Precoding
        [F_FLE_16,H_FLE_16,x_16,z_16]= Flexible_Precoding(eye(K),H_16',G_16,N,alpha,phi_t,theta_t,beta,K,L,lambda,I,d,d_min);
        [F_FLE_64,H_FLE_64,x_64,z_64]= Flexible_Precoding(eye(K),H_64',G_64,N,alpha,phi_t,theta_t,beta,K,L,lambda,I,d,d_min);
         
        
        for k=1:K
            F_FLE_16(:,k) = F_FLE_16(:,k) / norm(F_FLE_16(:,k)); % normalization 
            F_FLE_64(:,k) = F_FLE_64(:,k) / norm(F_FLE_64(:,k));  
        end
        %% Antenna Selection-  Fast greedy selection
        H_sel_16= AS_SEL(H_16,N,G_16,K,P_max,sigma_2); 
        H_sel_64= AS_SEL(H_64,N,G_64,K,P_max,sigma_2);
        
        F_SEL_16= H_sel_16/(H_sel_16'*H_sel_16+alpha*eye(K));
        H_SEL_16=H_sel_16';
    	F_SEL_64= H_sel_64/(H_sel_64'*H_sel_64+alpha*eye(K));
        H_SEL_64=H_sel_64';
        for k=1:K
            F_SEL_16(:,k) = F_SEL_16(:,k) / norm(F_SEL_16(:,k)); % normalization
            F_SEL_64(:,k) = F_SEL_64(:,k) / norm(F_SEL_64(:,k)); 
        end
         
        
        %% SINR-sum rate calculation
        for k=1:K
            IUI_FIX=0; 
            IUI_FLE_16=0; 
            IUI_FLE_64=0; 
            IUI_SEL_16=0;
            IUI_SEL_64=0;
            for kk=1:K
                if kk~=k
                    IUI_FIX=IUI_FIX+abs(H_FIX(k,:)*F_FIX(:,kk))^2; 
                    IUI_FLE_16=IUI_FLE_16+abs(H_FLE_16(k,:)*F_FLE_16(:,kk))^2; 
                    IUI_FLE_64=IUI_FLE_64+abs(H_FLE_64(k,:)*F_FLE_64(:,kk))^2; 
                    IUI_SEL_16=IUI_SEL_16+abs(H_SEL_16(k,:)*F_SEL_16(:,kk))^2;
                    IUI_SEL_64=IUI_SEL_64+abs(H_SEL_64(k,:)*F_SEL_64(:,kk))^2;
                end
            end
            SINR_FIX(k) =abs(H_FIX(k,:)*F_FIX(:,k))^2/(sigma_2+IUI_FIX); 
            SINR_FLE_16(k) =abs(H_FLE_16(k,:)*F_FLE_16(:,k))^2/(sigma_2+IUI_FLE_16); 
            SINR_FLE_64(k) =abs(H_FLE_64(k,:)*F_FLE_64(:,k))^2/(sigma_2+IUI_FLE_64); 
            SINR_SEL_16(k) =abs(H_SEL_16(k,:)*F_SEL_16(:,k))^2/(sigma_2+IUI_SEL_16);
            SINR_SEL_64(k) =abs(H_SEL_64(k,:)*F_SEL_64(:,k))^2/(sigma_2+IUI_SEL_64);
        end
        rate_FIX(reali,ll)=sum(log2(1+P_max/K*SINR_FIX)); % equal power allocation 
        rate_FLE_16(reali,ll)=sum(log2(1+P_max/K*SINR_FLE_16)); 
        rate_FLE_64(reali,ll)=sum(log2(1+P_max/K*SINR_FLE_64)); 
        rate_SEL_16(reali,ll)=sum(log2(1+P_max/K*SINR_SEL_16));
        rate_SEL_64(reali,ll)=sum(log2(1+P_max/K*SINR_SEL_64));
    end
end
rate_FIX_mean=mean(rate_FIX) 
rate_FLE_16_mean=mean(rate_FLE_16) 
rate_FLE_64_mean=mean(rate_FLE_64) 
rate_SEL_16_mean=mean(rate_SEL_16)
rate_SEL_64_mean=mean(rate_SEL_64)

co1= [0, 161, 241]/255;
co2=[29, 191, 151]/255
co3= [70, 158, 180]/255
co4=[253,185,106]/255
co5=[214,64,78]/255
figure
plot(LL, abs(rate_FIX_mean), '^k-', 'linewidth', 1.1, 'markerfacecolor', co2,'markersize', 7.5)
hold on  
plot(LL, abs(rate_SEL_16_mean), 'ok-', 'linewidth', 1.1, 'markerfacecolor', co1,'markersize', 7.5)
hold on
plot(LL, abs(rate_SEL_64_mean), 'ok--', 'linewidth', 1.1, 'markerfacecolor', co1,'markersize', 7.5)
hold on
plot(LL, abs(rate_FLE_16_mean), 'sk-', 'linewidth', 1.1, 'markerfacecolor', co5,'markersize', 7.5)
hold on 
plot(LL, abs(rate_FLE_64_mean), 'sk--', 'linewidth', 1.1, 'markerfacecolor', co5,'markersize', 7.5)
 
grid on
lgh=legend('Fixed Antenna Position','Fast AS, $G=16$','Fast AS, $G=64$','Flexible Precoding, $G=16$','Flexible Precoding, $G=64$');
set(lgh,'interpreter','latex', 'fontsize', 14);
xlabel('Number of Channel Paths $L$','interpreter','latex','fontsize',14)
ylabel('Sum Rate [bit/s/Hz]','interpreter','latex','fontsize',14)
axis([1,15,2,9])