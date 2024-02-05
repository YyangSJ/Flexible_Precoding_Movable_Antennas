clear all;
close all;
rng(0)
f=3e9; % system frequency
lambda=3e8/f; % antenna wavelength
d=lambda/2; % inter-element spacing for virtual channel representation
d_min=lambda/2; % allowed minimal inter-element spacing in practice
K=4; % number of users
N=4; % number of movable antennas
GG=[9,16,25,36,49,64,81,100]; % numbe of virtual antennas
L=3; % number of channel paths
sigma_2=1; % noise power
P_max=1; % transmit power
alpha=1; % regularized factor
I=100; % iteration number for OffG procedure
realization=666;
for reali=1:realization
    beta=normrnd(0,1,L,K)+1i*normrnd(0,1,L,K); % path gains
    phi_t=unifrnd(-1,1,L,K); % BS path DoA-azimuth
    theta_t=unifrnd(-1,1,L,K); % BS path DoA-elevation
    phi_r=zeros(L,K);  theta_r=zeros(L,K); % user antenna is fixed, the angles are set to 0.
    for gg=1:length(GG)
        H=[];
        G=GG(gg); 
        for k=1:K
            H(:,k)=Channel(beta(:,k),phi_r(:,k),theta_r(:,k),phi_t(:,k),theta_t(:,k),lambda,d,G,L,0,0);
        end
        %% Fixed Antenna Position for Precoding
        array_sel=1:N;
        H_FIX=H(array_sel,:);
        F_FIX= H_FIX/(H_FIX'*H_FIX+alpha*eye(K));
        pow_FIX=trace(F_FIX'*F_FIX)
        H_FIX=H_FIX';
        %% Flexible Precoding
        [F_FLE,H_FLE,x,z]= Flexible_Precoding(eye(K),H',G,N,alpha,phi_t,theta_t,beta,K,L,lambda,I,d,d_min);
        pow_FFZF=trace(F_FLE'*F_FLE)
        for k=1:K
            F_FIX(:,k) = F_FIX(:,k) / norm(F_FIX(:,k)); % normalization
        end
        
        for k=1:K
            F_FLE(:,k) = F_FLE(:,k) / norm(F_FLE(:,k)); % normalization
        end
        %% Antenna Selection-  Fast greedy selection
        H_sel= AS_SEL(H,N,G,K,P_max,sigma_2);
      
        F_SEL= H_sel/(H_sel'*H_sel+alpha*eye(K));
          H_SEL=H_sel';
        
        for k=1:K
            F_SEL(:,k) = F_SEL(:,k) / norm(F_SEL(:,k)); % normalization
        end
        
        
        %% SINR-sum rate calculation
        for k=1:K
            IUI_FIX=0;
            IUI_FLE=0;
            IUI_SEL=0;
            for kk=1:K
                if kk~=k
                    IUI_FIX=IUI_FIX+abs(H_FIX(k,:)*F_FIX(:,kk))^2;
                    IUI_FLE=IUI_FLE+abs(H_FLE(k,:)*F_FLE(:,kk))^2;
                    IUI_SEL=IUI_SEL+abs(H_SEL(k,:)*F_SEL(:,kk))^2;
                end
            end
            SINR_FIX(k) =abs(H_FIX(k,:)*F_FIX(:,k))^2/(sigma_2+IUI_FIX);
            SINR_FLE(k) =abs(H_FLE(k,:)*F_FLE(:,k))^2/(sigma_2+IUI_FLE);
            SINR_SEL(k) =abs(H_SEL(k,:)*F_SEL(:,k))^2/(sigma_2+IUI_SEL);
        end
        rate_FIX(reali,gg)=sum(log2(1+P_max/K*SINR_FIX)); % equal power allocation
        rate_FLE(reali,gg)=sum(log2(1+P_max/K*SINR_FLE));
        rate_SEL(reali,gg)=sum(log2(1+P_max/K*SINR_SEL));
    end
end
rate_FIX_mean=mean(rate_FIX)
rate_FLE_mean=mean(rate_FLE)
rate_SEL_mean=mean(rate_SEL)

co1= [0, 161, 241]/255;
co2=[29, 191, 151]/255
co3= [70, 158, 180]/255
co4=[253,185,106]/255
co5=[214,64,78]/255
figure
plot(GG, abs(rate_FIX_mean), '^k-', 'linewidth', 1.1, 'markerfacecolor', co2,'markersize', 7.5)
hold on
plot(GG, abs(rate_SEL_mean), 'sk-', 'linewidth', 1.1, 'markerfacecolor', co5,'markersize', 7.5)
hold on
plot(GG, abs(rate_FLE_mean), 'ok-', 'linewidth', 1.1, 'markerfacecolor', co1,'markersize', 7.5)

grid on
lgh=legend('Fixed Antenna Position','Fast AS','Flexible Precoding');
set(lgh,'interpreter','latex', 'fontsize', 14);
xlabel('Size of Movable Region $G$','interpreter','latex','fontsize',14)
ylabel('Sum Rate [bit/s/Hz]','interpreter','latex','fontsize',14)
