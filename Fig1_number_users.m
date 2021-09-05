clear all;
close all;
clc;


%---------- coordinates to create a cell ----------%


%the minimum x and y coordinates of BD (backscatter device) from BS
b1=0; 
%the maximum x and y coordinates of BD from BS
b2=10; 

%the minimum x and y coordinates of users from BS
a1=10; 
%the maximum x and y coordinates of users from BS
a2=20; 

%---------- input parameters ----------%


%pathloss exponent
alpha_path=3;

%power budget in dbm
Pmax_dbm=80;

%circuit power in dbm
Pc_dbm=30; 

%noise power in dbm
sigma_dbm=-100;

%call the function 'dbm_to_Watt' to convert from dbm to Watt
Pmax=dbm_to_Watt(Pmax_dbm);
Pc=dbm_to_Watt(Pc_dbm);
sigma=dbm_to_Watt(sigma_dbm);

%minimum  rate QoS constraint 
Rmin=1;

%number of realizations
N=1e5;


%-------- compute the GEE as a function of the number of users K ----------%

%number of users 
K=2:10;

%for each channel realization 
for n=1:N
    %for each number of users K
    for i=1:length(K)
        
        A=(2.^(2*Rmin)).*ones(K(i),1);
        
        %---- generating channel gains ----%
        
        %randomly generate x and y coordinates for users
        users_location=(a2-a1).*rand(2,K(i))+a1; 
        %randomly generate x and y coordinates for BD
        BD_location=(b2-b1).*rand(2,1)+b1; 
        
        %compute the distance BS-BD
        d_BS_BD=sqrt(sum(BD_location.^2)); 
        %compute the distance BS-users
        d_BS_users=sqrt(sum(users_location.^2))'; 
        %compute the distance BD-users
        d_BD_users=sqrt(sum((ones(1,K(i)).*BD_location-users_location).^2))'; 
        
        
        %generate BS-BD channel gain
        G_BS_BD=channel_gain_pathloss(d_BS_BD,alpha_path,sigma);
        %generate BS-users channel gains (in descending order -> SIC)
        [G_BS_users,I]=channel_gain_pathloss(d_BS_users,alpha_path,sigma);
        %generate BD-users channel gains before SIC order
        G_BD_users_unordered=channel_gain_pathloss_BD(d_BD_users,alpha_path);
        %ordering channel gains BD-users
        G_BD_users=G_BD_users_unordered(I);
                             
              
        
        %---- calling the function that computes R according to equation (6) in the paper ----%
        R = rho_plus(G_BS_users,G_BS_BD,G_BD_users);
        
        %---- compute the minimum power budget (Pmin) required for meeting QoS constraints in conventional NOMA ----%
        Pmin_NOMA_conv=(A(end)-1)/G_BS_users(end);
        for j=1:(K(i)-1)
            Pmin_NOMA_conv=Pmin_NOMA_conv+(A(j)-1)/G_BS_users(j)*prod(A((j+1):K(i)));
        end
        
        %---- compute Pmin required for meeting QoS constraints in conventional OMA ----%
        Pmin_OMA_conv=sum((A.^K(i)-1)./G_BS_users);
        
        %compute Pmin in NOMA with BD and OMA with BD + computing optimal rho 
        [rho_NOMA,G_NOMA_BD,Pmin_NOMA_BD,G_OMA_BD,Pmin_OMA_BD] = optimal_rho(G_BS_users,G_BS_BD,G_BD_users,A,R);
        
        
        %---- checking the feasability condition ----%
        while (Pmin_NOMA_conv>Pmax || Pmin_NOMA_BD>Pmax || Pmin_OMA_BD>Pmax || Pmin_OMA_conv>Pmax)
            
            %repeat generating channel gains until the feasibility condition is met
            users_location=(a2-a1).*rand(2,K(i))+a1; 
            BD_location=(b2-b1).*rand(2,1)+b1;
            
            d_BS_BD=sqrt(sum(BD_location.^2));
            d_BS_users=sqrt(sum(users_location.^2))';
            d_BD_users=sqrt(sum((ones(1,K(i)).*BD_location-users_location).^2))'; 
            
            
            G_BS_BD=channel_gain_pathloss(d_BS_BD,alpha_path,sigma);
            [G_BS_users,I]=channel_gain_pathloss(d_BS_users,alpha_path,sigma);
            G_BD_users_unordered=channel_gain_pathloss_BD(d_BD_users,alpha_path);
            G_BD_users=G_BD_users_unordered(I);
            
            R = rho_plus(G_BS_users,G_BS_BD,G_BD_users);
            
            Pmin_NOMA_conv=(A(end)-1)/G_BS_users(end);
            for j=1:(K(i)-1)
                Pmin_NOMA_conv=Pmin_NOMA_conv+(A(j)-1)/G_BS_users(j)*prod(A((j+1):K(i)));
            end
            Pmin_OMA_conv=sum((A.^K(i)-1)./G_BS_users);
            
            %compute the optimal reflection coefficient(rho)
            
            if (isempty(R))
                rho_NOMA=1;
            else
                rho_NOMA=min(1,min(R));
            end
            
            %compute Gamma according to the notations in the paper
            G_OMA_BD=(G_BS_users+G_BS_BD*G_BD_users).^2;
            G_NOMA_BD=(G_BS_users+sqrt(rho_NOMA)*G_BS_BD*G_BD_users).^2;
           
            %compute Pmin in OMA with BD and NOMA with BD
            Pmin_OMA=sum((A.^length(G_BS_users)-1)./G_OMA_BD);
            Pmin_NOMA=0;
            for i=1:length(G_BS_users)
                Pmin_NOMA=Pmin_NOMA+(A(i)-1)/G_NOMA_BD(i)*prod(A(i+1:length(G_BS_users)));
            end             
            
            %[rho_NOMA,G_NOMA_BD,Pmin_NOMA_BD,G_OMA_BD,Pmin_OMA_BD] = optimal_rho(G_BS_users,G_BS_BD,G_BD_users,A,R);
            
        end
        
        %---- compute the optimal energy efficiency ----% 
        EE_NOMA_BD = optimal_solution_NOMA(G_NOMA_BD,A,Pmax,Pmin_NOMA_BD,Pc);
        EE_NOMA_conv = optimal_solution_NOMA(G_BS_users,A,Pmax,Pmin_NOMA_conv,Pc);
        EE_OMA_BD = optimal_solution_OMA(G_OMA_BD,A,Pmax,Pc);
        EE_OMA_conv = optimal_solution_OMA(G_BS_users,A,Pmax,Pc);
        
        
        %---- stock results for each number of users ----%
        EE_op_NOMA_BD(i)=EE_NOMA_BD;
        EE_op_NOMA_conv(i)=EE_NOMA_conv;
        EE_op_OMA_BD(i)=EE_OMA_BD;
        EE_op_OMA_conv(i)=EE_OMA_conv;
        
    end
    %-------- stock results for each channel realization --------%
    EE_op_NOMA_BD_n(n,:)=EE_op_NOMA_BD;
    EE_op_NOMA_conv_n(n,:)=EE_op_NOMA_conv;
    EE_op_OMA_BD_n(n,:)=EE_op_OMA_BD;
    EE_op_OMA_conv_n(n,:)=EE_op_OMA_conv;
    
end

%------------- averaging over channel realizations -----------%
EE_op_NOMA_BD_mean=mean(EE_op_NOMA_BD_n);
EE_op_NOMA_conv_mean=mean(EE_op_NOMA_conv_n);
EE_op_OMA_BD_mean=mean(EE_op_OMA_BD_n);
EE_op_OMA_conv_mean=mean(EE_op_OMA_conv_n);

%------------- plot figures -----------%
figure(1)
plot(K,EE_op_NOMA_BD_mean,'-o','MarkerSize',7,'LineWidth',1.3,'MarkerIndices', 1:length(K));
hold on;
plot(K,EE_op_NOMA_conv_mean,'-d','MarkerSize',7,'LineWidth',1.3,'MarkerIndices', 1:length(K));
hold on;
plot(K,EE_op_OMA_BD_mean,'-s','MarkerSize',7,'LineWidth',1.3,'MarkerIndices', 1:length(K));
hold on;
plot(K,EE_op_OMA_conv_mean,'-*','MarkerSize',7,'LineWidth',1.3,'MarkerIndices', 1:length(K));
ylabel('GEE (bits/J)');
xlabel('Number of receivers');
legend('NOMA+backscatter','NOMA','OMA+backscatter','OMA','Location=Best');
grid on;
