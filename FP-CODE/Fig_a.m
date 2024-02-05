clear all;
close all;
rng(0)
f=3e9; % system frequency
lambda=3e8/f; % antenna wavelength
d=lambda/2; % inter-element spacing for virtual channel representation
d_min=lambda/2; % allowed minimal inter-element spacing in practice
K=4; % number of users
N=4; % number of movable antennas
G=36; % numbe of virtual antennas
L=15; % number of channel paths
sigma_2=1; % noise power
P_max=1; % transmit power
alphaa=[1e-2,1,1e2]; % regularized factor
I=100; % iteration number for OffG procedure
realization=10000;
for reali=1:realization
    beta=normrnd(0,1,L,K)+1i*normrnd(0,1,L,K); % path gains
    phi_t=unifrnd(-1,1,L,K); % BS path DoA-azimuth
    theta_t=unifrnd(-1,1,L,K); % BS path DoA-elevation
    phi_r=zeros(L,K);  theta_r=zeros(L,K); % user antenna is fixed, the angles are set to 0.
    for aa=1:length(alphaa)
        H=[];
        alpha=alphaa(aa); 
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
        rate_FIX(reali,aa)=sum(log2(1+P_max/K*SINR_FIX)); % equal power allocation
        rate_FLE(reali,aa)=sum(log2(1+P_max/K*SINR_FLE));
        rate_SEL(reali,aa)=sum(log2(1+P_max/K*SINR_SEL));
    end
end
for aa=1:3
[f_FIX(:,aa), x_FIX(:,aa)] = ecdf(rate_FIX(:,aa));   
[f_SEL(:,aa), x_SEL(:,aa)] = ecdf(rate_SEL(:,aa));  
[f_FLE(:,aa), x_FLE(:,aa)] = ecdf(rate_FLE(:,aa));  
end

co1= [0, 161, 241]/255;
co2=[29, 191, 151]/255
co3= [70, 158, 180]/255
co4=[253,185,106]/255
co5=[214,64,78]/255
figure
plot(x_FIX(:,1), f_FIX(:,1), '--', 'linewidth', 1.1, 'color',co2)
hold on
plot(x_FIX(:,2), f_FIX(:,2), '-', 'linewidth', 1.1,'color',co2)
hold on
plot(x_FIX(:,3), f_FIX(:,3), '-.', 'linewidth', 1.1, 'color',co2)
hold on
plot(x_SEL(:,1), f_SEL(:,1), '--', 'linewidth', 1.1,'color',co1)
hold on
plot(x_SEL(:,2), f_SEL(:,2), '-', 'linewidth', 1.1, 'color', co1)
hold on
plot(x_SEL(:,3), f_SEL(:,3), '-.', 'linewidth', 1.1, 'color', co1)
hold on
plot(x_FLE(:,1), f_FLE(:,1), '--', 'linewidth', 1.1, 'color', co5)
hold on
plot(x_FLE(:,2), f_FLE(:,2), '-', 'linewidth', 1.1,  'color', co5)
hold on
plot(x_FLE(:,3), f_FLE(:,3), '-.', 'linewidth', 1.1,  'color', co5)
hold on


 
lgh=legend('Fixed Antenna Position, $\alpha=10^{-2}$','Fixed Antenna Position, $\alpha=1$',...,
    'Fixed Antenna Position, $\alpha=10^2$','Fast AS, $\alpha=10^{-2}$','Fast AS, $\alpha=1$',...,
    'Fast AS, $\alpha=10^2$','Flexible Precoding, $10^{-2}$','Flexible Precoding, 1','Flexible Precoding, $10^{2}$');
set(lgh,'interpreter','latex', 'fontsize', 14);
xlabel('Sum Rate [bit/s/Hz]','interpreter','latex','fontsize',14)
ylabel('CDF','interpreter','latex','fontsize',14)
