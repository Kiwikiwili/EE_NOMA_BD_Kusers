function [alpha,Rsum,Ptot] = optimal_solution_OMA(G,A,Pmax,Pc)

epsilon=1e-7;
alpha=1e-5;
R_P=1;

P1=(A.^length(G)-1)./G;

P3=Pmax*ones(length(G),1);

P1_1=(A(1)^length(G)-1)/G(1);
P2_1=(A(2)^length(G)-1)/G(2);

P1_3=Pmax;
P2_3=Pmax;

while(R_P>epsilon)

    P2=1/(2*log(2)*alpha)-1./G;
    P1_2=1/(2*log(2)*alpha)-1/G(1);
    P2_2=1/(2*log(2)*alpha)-1/G(2);
    
    
    p1=min(P1_3,max(P1_1,P1_2));
    p2=min(P2_3,max(P2_1,P2_2));

    for k=1:length(G)
        p(k)=min(P3(k),max(P1(k),P2(k)));
    end
    
    R_P=(1/(2*length(G)))*sum(log2(1+G.*p'))-alpha*(sum(p)/length(G)+Pc);
    alpha=(1/(2*length(G)))*sum(log2(1+G.*p'))/(sum(p)/length(G)+Pc);
    
    Rsum=(1/(2*length(G)))*sum(log2(1+G.*p'));
    Ptot=(sum(p)/length(G))+Pc;
end
end