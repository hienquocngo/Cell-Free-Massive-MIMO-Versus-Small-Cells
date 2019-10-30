clc
clear all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Downlink
%Consider a square are of DxD m^2
%M distributed APs serves K terminals, they all randomly located in the area
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Inital parameters
M=100; %number of access points
K=40; %number of terminals

D=1; %in kilometer
tau=20;%training length
[U,S,V]=svd(randn(tau,tau));%U includes tau orthogonal sequences 

B=20; %Mhz
Hb = 15; % Base station height in m
Hm = 1.65; % Mobile height in m
f = 1900; % Frequency in MHz
aL = (1.1*log10(f)-0.7)*Hm-(1.56*log10(f)-0.8);
L = 46.3+33.9*log10(f)-13.82*log10(Hb)-aL;

power_f=0.2; %downlink power: 200 mW
noise_p = 10^((-203.975+10*log10(20*10^6)+9)/10); %noise power
Pd = power_f/noise_p;%nomalized receive SNR
Pp=Pd;%pilot power

d0=0.01;%km
d1=0.05;%km

N=200;
R_cf_min=zeros(1,N);%min rate, cell-free, without power allocation
R_cf_opt_min=zeros(1,N); %min rate, cell-free, with power allocation

R_sc_min=zeros(1,N);%small cell
R_sc_opt_min=zeros(1,N);

R_cf_user=zeros(N,K);
R_sc_user=zeros(N,K);

for n=1:N
n
%%%%%Randomly locations of M APs%%%%
AP=zeros(M,2,9);
AP(:,:,1)=unifrnd(-D/2,D/2,M,2);

%Wrapped around (8 neighbor cells)
D1=zeros(M,2);
D1(:,1)=D1(:,1)+ D*ones(M,1);
AP(:,:,2)=AP(:,:,1)+D1;

D2=zeros(M,2);
D2(:,2)=D2(:,2)+ D*ones(M,1);
AP(:,:,3)=AP(:,:,1)+D2;

D3=zeros(M,2);
D3(:,1)=D3(:,1)- D*ones(M,1);
AP(:,:,4)=AP(:,:,1)+D3;

D4=zeros(M,2);
D4(:,2)=D4(:,2)- D*ones(M,1);
AP(:,:,5)=AP(:,:,1)+D4;

D5=zeros(M,2);
D5(:,1)=D5(:,1)+ D*ones(M,1);
D5(:,2)=D5(:,2)- D*ones(M,1);
AP(:,:,6)=AP(:,:,1)+D5;

D6=zeros(M,2);
D6(:,1)=D6(:,1)- D*ones(M,1);
D6(:,2)=D6(:,2)+ D*ones(M,1);
AP(:,:,7)=AP(:,:,1)+D6;

D7=zeros(M,2);
D7=D7+ D*ones(M,2);
AP(:,:,8)=AP(:,:,1)+D7;

D8=zeros(M,2);
D8=D8- D*ones(M,2);
AP(:,:,9)=AP(:,:,1)+D8;

%Randomly locations of K terminals:
Ter=zeros(K,2,9);
Ter(:,:,1)=unifrnd(-D/2,D/2,K,2);

%Wrapped around (8 neighbor cells)
D1=zeros(K,2);
D1(:,1)=D1(:,1)+ D*ones(K,1);
Ter(:,:,2)=Ter(:,:,1)+D1;

D2=zeros(K,2);
D2(:,2)=D2(:,2)+ D*ones(K,1);
Ter(:,:,3)=Ter(:,:,1)+D2;

D3=zeros(K,2);
D3(:,1)=D3(:,1)- D*ones(K,1);
Ter(:,:,4)=Ter(:,:,1)+D3;

D4=zeros(K,2);
D4(:,2)=D4(:,2)- D*ones(K,1);
Ter(:,:,5)=Ter(:,:,1)+D4;

D5=zeros(K,2);
D5(:,1)=D5(:,1)+ D*ones(K,1);
D5(:,2)=D5(:,2)- D*ones(K,1);
Ter(:,:,6)=Ter(:,:,1)+D5;

D6=zeros(K,2);
D6(:,1)=D6(:,1)- D*ones(K,1);
D6(:,2)=D6(:,2)+ D*ones(K,1);
Ter(:,:,7)=Ter(:,:,1)+D6;

D7=zeros(K,2);
D7=D7+ D*ones(K,2);
Ter(:,:,8)=Ter(:,:,1)+D7;

D8=zeros(K,2);
D8=D8- D*ones(K,2);
Ter(:,:,9)=Ter(:,:,1)+D8;

sigma_shd=8; %in dB
D_cor=0.1;

%%%%%%Create the MxK correlated shadowing matrix %%%%%%%

    %%%%M correlated shadowing cofficients of M APs:
    Dist=zeros(M,M);%distance matrix
    Cor=zeros(M,M);%correlation matrix

    for m1=1:M
        for m2=1:M
            Dist(m1,m2) = min([norm(AP(m1,:,1)-AP(m2,:,1)), norm(AP(m1,:,1)-AP(m2,:,2)),norm(AP(m1,:,1)-AP(m2,:,3)),norm(AP(m1,:,1)-AP(m2,:,4)),norm(AP(m1,:,1)-AP(m2,:,5)),norm(AP(m1,:,1)-AP(m2,:,6)),norm(AP(m1,:,1)-AP(m2,:,7)),norm(AP(m1,:,1)-AP(m2,:,8)),norm(AP(m1,:,1)-AP(m2,:,9)) ]); %distance between AP m1 and AP m2
            Cor(m1,m2)=exp(-log(2)*Dist(m1,m2)/D_cor);
        end
    end
    A1 = chol(Cor,'lower');
    x1 = randn(M,1);
    sh_AP = A1*x1;
    for m=1:M
        sh_AP(m)=(1/sqrt(2))*sigma_shd*sh_AP(m)/norm(A1(m,:));
    end

    %%%%K correlated shadowing matrix of K terminal:
    Dist=zeros(K,K);%distance matrix
    Cor=zeros(K,K);%correlation matrix

    for k1=1:K
        for k2=1:K
            Dist(k1,k2)=min([norm(Ter(k1,:,1)-Ter(k2,:,1)), norm(Ter(k1,:,1)-Ter(k2,:,2)),norm(Ter(k1,:,1)-Ter(k2,:,3)),norm(Ter(k1,:,1)-Ter(k2,:,4)),norm(Ter(k1,:,1)-Ter(k2,:,5)),norm(Ter(k1,:,1)-Ter(k2,:,6)),norm(Ter(k1,:,1)-Ter(k2,:,7)),norm(Ter(k1,:,1)-Ter(k2,:,8)),norm(Ter(k1,:,1)-Ter(k2,:,9)) ]); %distance between Terminal k1 and Terminal k2
            Cor(k1,k2)=exp(-log(2)*Dist(k1,k2)/D_cor);
        end
    end
    A2 = chol(Cor,'lower');
    x2 = randn(K,1);
    sh_Ter = A2*x2;
    for k=1:K
        sh_Ter(k)=(1/sqrt(2))*sigma_shd*sh_Ter(k)/norm(A2(k,:));
    end

%%% The shadowing matrix:
Z_shd=zeros(M,K);
for m=1:M
    for k=1:K
        Z_shd(m,k)= sh_AP(m)+ sh_Ter(k);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Create an MxK large-scale coefficients beta_mk
BETAA = zeros(M,K);
dist=zeros(M,K);
for m=1:M  
    for k=1:K
    [dist(m,k),index] = min([norm(AP(m,:,1)-Ter(k,:,1)), norm(AP(m,:,2)-Ter(k,:,1)),norm(AP(m,:,3)-Ter(k,:,1)),norm(AP(m,:,4)-Ter(k,:,1)),norm(AP(m,:,5)-Ter(k,:,1)),norm(AP(m,:,6)-Ter(k,:,1)),norm(AP(m,:,7)-Ter(k,:,1)),norm(AP(m,:,8)-Ter(k,:,1)),norm(AP(m,:,9)-Ter(k,:,1)) ]); %distance between Terminal k and AP m
    if dist(m,k)<d0
         betadB=-L - 35*log10(d1) + 20*log10(d1) - 20*log10(d0);
    elseif ((dist(m,k)>=d0) && (dist(m,k)<=d1))
         betadB= -L - 35*log10(d1) + 20*log10(d1) - 20*log10(dist(m,k));
    else
    betadB = -L - 35*log10(dist(m,k)) + Z_shd(m,k); %large-scale in dB
    end
    
    BETAA(m,k)=10^(betadB/10); 
    end

end

%% Pilot Asignment: (random choice)
Phii=zeros(tau,K);
for k=1:K
    Point=randi([1,tau]);
    Phii(:,k)=U(:,Point);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%                       CELL FREE                     %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Phii_cf = Phii; % pilot set of cell-free systems
%% Create Gamma matrix (variances of the channel estimates)
Gammaa = zeros(M,K);
mau=zeros(M,K);
for m=1:M
    for k=1:K
        mau(m,k)=norm( (BETAA(m,:).^(1/2)).*(Phii_cf(:,k)'*Phii_cf))^2;
    end
end

for m=1:M
    for k=1:K
        Gammaa(m,k)=tau*Pp*BETAA(m,k)^2/(tau*Pp*mau(m,k) + 1);
    end
end

%% 1) Each AP has equal power allocations for K terminals
%%Compute etaa(m): (each AP transmits equal power to K terminals)
etaa=zeros(M,1);
for m=1:M
    etaa(m)=1/(sum(Gammaa(m,:)));
end

%%Compute Rate
SINR=zeros(1,K);
R_cf=zeros(1,K);

%Pilot contamination
PC = zeros(K,K);
for ii=1:K
    for k=1:K
        PC(ii,k)=sum((etaa.^(1/2)).*((Gammaa(:,ii)./BETAA(:,ii)).*BETAA(:,k)))*Phii_cf(:,ii)'*Phii_cf(:,k);
    end
end
PC1=(abs(PC)).^2;

for k=1:K
    num=0;
    for m=1:M
        num=num + (etaa(m)^(1/2))*Gammaa(m,k);
    end
    SINR(k) = Pd*num^2/(1 + Pd*sum(BETAA(:,k)) + Pd*sum(PC1(:,k)) - Pd*PC1(k,k));
    %Rate:
    R_cf(k) = log2(1+ SINR(k));
end

%%Greedy Pilot Assignment:
stepp=5;
Ratestep=zeros(stepp,K);
Ratestep(1,:)=R_cf;

%%Find the pilot sequences using greedy method
for st=2:stepp
    [minvalue minindex]=min(Ratestep(st-1,:));
    
        Mat=zeros(tau,tau)-Pd*sum(BETAA(:,minindex))*Phii_cf(:,minindex)*Phii_cf(:,minindex)';
        for mm=1:M
            for kk=1:K
                Mat = Mat + Pd*BETAA(mm,kk)*Phii_cf(:,kk)*Phii_cf(:,kk)';
            end
        end
        [U1,S1,V1] = svd(Mat);
        Phii_cf(:,minindex) = U1(:,tau);
    
 
  %% Create Gamma matrix with new pilot sequences
Gammaa = zeros(M,K);
mau=zeros(M,K);
for m=1:M
    for k=1:K
        mau(m,k)=norm( (BETAA(m,:).^(1/2)).*(Phii_cf(:,k)'*Phii_cf))^2;
    end
end

for m=1:M
    for k=1:K
        Gammaa(m,k)=tau*Pp*BETAA(m,k)^2/(tau*Pp*mau(m,k) + 1);
    end
end

%% Compute etaa(m): (each AP transmits equal power to K terminals)
etaa=zeros(M,1);
for m=1:M
    etaa(m)=1/(sum(Gammaa(m,:)));
end

%% Compute Rate
SINR=zeros(1,K);
R_cf=zeros(1,K);

%Pilot contamination
PC = zeros(K,K);
for ii=1:K
    for k=1:K
        PC(ii,k)=sum((etaa.^(1/2)).*((Gammaa(:,ii)./BETAA(:,ii)).*BETAA(:,k)))*Phii_cf(:,ii)'*Phii_cf(:,k);
    end
end
PC1=(abs(PC)).^2;

for k=1:K
    num=0;
    for m=1:M
        num=num + (etaa(m)^(1/2))*Gammaa(m,k);
    end
    SINR(k) = Pd*num^2/(1 + Pd*sum(BETAA(:,k)) + Pd*sum(PC1(:,k)) - Pd*PC1(k,k));
    %Rate:
    Ratestep(st,k) = log2(1+ SINR(k));
end      
    

end

R_cf_min(n)=min(Ratestep(stepp,:));
R_cf_user(n,:)=Ratestep(stepp,:);

%% 2) Max-Min power allocations%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tmin=2^R_cf_min(n)-1;
tmax=2^(2*R_cf_min(n)+1.2)-1;
epsi=max(tmin/5,0.01);
PhiPhi=zeros(K,K);
for ii=1:K
    for k=1:K
        PhiPhi(ii,k)=norm(Phii_cf(:,ii)'*Phii_cf(:,k));
    end
end
BETAAn=BETAA*Pd;
Gammaan=Gammaa*Pd;
%cvx_solver sedumi
cvx_quiet true
            while( tmax - tmin > epsi)

            tnext = (tmax+tmin)/2; 
           cvx_begin sdp
              variables X(M,K) y(M,1) Z(K,K)
              minimize(0)
              subject to
                for k=1:K
                     norm([sqrt(tnext)*[Z((1:(k-1)),k);Z(((k+1):K),k)].*[PhiPhi((1:(k-1)),k);PhiPhi(((k+1):K),k)]  ;  sqrt(tnext)*(BETAAn(:,k).^(1/2)).*y ; sqrt(tnext*Pd)]) <= (Gammaan(:,k))'*X(:,k) ;
                end
                for m=1:M
                    norm(((Gammaan(m,:)).^(1/2)).*X(m,:)) <= y(m); 
                    y(m)<=sqrt(Pd);
                end
                
                for k=1:K
                    for ii=1:K
                        sum( ((Gammaan(:,ii)./BETAA(:,ii)).*BETAA(:,k)).*X(:,ii)  ) <= Z(ii,k);
                    end
                end
                
                for m=1:M
                 for k=1:K
                    X(m,k)>=0;
                 end
                end
               
            cvx_end


            % bisection
            if strfind(cvx_status,'Solved') % feasible
            fprintf(1,'Problem is feasible ',tnext);
            tmin = tnext;
            else % not feasible
            fprintf(1,'Problem not feasible ',tnext);
            tmax = tnext;   
            end

            end

R_cf_opt_min(n) = log2(1+tmin);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%                       SMALL CELL                     %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Each terminal chooses an AP which has largest large-scale fading
m_index = zeros(1,K); %indexes of K chosen APs 
BETAA1=BETAA;
for k=1:K
  [Betamax m_index(k)]=max(BETAA1(:,k));
  BETAA1(m_index(k),:)=zeros(1,K);
end
Pd1=Pd;
Phii_sc=Phii; %pilot set of small-cell systems
%% 1) No power control (the selected APs use full power)
MUU=zeros(1,K);
mauu=zeros(1,K);

for k=1:K
    for kk=1:K
        mauu(k)=mauu(k) + BETAA(m_index(kk),k)*(norm(Phii_sc(:,k)'*Phii_sc(:,kk)))^2;
    end
end

for k=1:K
    MUU(k) = tau*Pp*BETAA(m_index(k),k)^2/(tau*Pp*mauu(k) + 1);
end

%% Compute rate
R_sc=zeros(1,K);
gam=zeros(1,K);
for k=1:K
   a=-Pd1*MUU(k);
   for m=1:K
   a=a + Pd1*BETAA(m_index(m),k);
   end
   gam(k) = Pd1*MUU(k)/(1+a);
end

for k=1:K
R_sc(k) = (1/log(2))*exp(1/gam(k))*expint(1/gam(k));
end

stepp=5;
Ratestep=zeros(stepp,K);
Ratestep(1,:)=R_sc;

%% Find the pilot sequences using greedy method
for st=2:stepp
    [minvalue minindex]=min(Ratestep(st-1,:));
    
        Mat=zeros(tau,tau)-Pd1*BETAA(m_index(minindex),minindex)*Phii_sc(:,minindex)*Phii_sc(:,minindex)';
        for kk=1:K
                Mat = Mat + Pd1*BETAA(m_index(kk),minindex)*Phii_sc(:,kk)*Phii_sc(:,kk)';
        end
        [U1,S1,V1] = svd(Mat);
        Phii_sc(:,minindex) = U1(:,tau);
    
 
MUU=zeros(1,K);
mauu=zeros(1,K);

for k=1:K
    for kk=1:K
        mauu(k)=mauu(k) + BETAA(m_index(kk),k)*(norm(Phii_sc(:,k)'*Phii_sc(:,kk)))^2;
    end
end

for k=1:K
    MUU(k) = tau*Pp*BETAA(m_index(k),k)^2/(tau*Pp*mauu(k) + 1);
end

  %% Computate rate
R_sc=zeros(1,K);
gam=zeros(1,K);
for k=1:K
   a=-Pd1*MUU(k);
   for m=1:K
   a=a + Pd1*BETAA(m_index(m),k);
   end
   gam(k) = Pd1*MUU(k)/(1+a);
end

for k=1:K
Ratestep(st,k) = (1/log(2))*exp(1/gam(k))*expint(1/gam(k));
end

       
end

R_sc_min(n)=min(Ratestep(stepp,:));
R_sc_user(n,:)=Ratestep(stepp,:);

%% 2) Max-Min power allocations%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gam_min=min(gam);
tmin=gam_min;
tmax=2*tmin + 1;
epsi=max(tmin/5,0.01);

BETAAn=BETAA*Pd1;
MUUn=MUU*Pd1;

%cvx_solver sedumi
cvx_quiet true
            while( tmax - tmin > epsi)

            tnext = (tmax+tmin)/2; 
           cvx_begin sdp
              variables x(K,1) 
              minimize(0)
              subject to
                for k=1:K
                   BETAAn(m_index,k)'*x  - MUUn(k)*x(k) +  1 <= (1/tnext)*MUUn(k)*x(k) ;
                end              
                for k=1:K
                    x(k)>=0;
                    x(k)<=1;
                end
                            
            cvx_end


            % bisection
            if strfind(cvx_status,'Solved') % feasible
            fprintf(1,'Problem is feasible ',tnext);
            tmin = tnext;
            else % not feasible
            fprintf(1,'Problem not feasible ',tnext);
            tmax = tnext;   
            end

            end

R_sc_opt_min(n) = (1/log(2))*exp(1/tmin)*expint(1/tmin);

end

Y=linspace(0,1,N);
hold on 

plot(sort(R_cf_opt_min),Y(:),'r')
plot(sort(R_sc_min),Y(:),'--b')







