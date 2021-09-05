%------ Dinkelbach's algorithm -------% 

%% input: channel gain (G), A=2^(2*Rmin), maximum power budget (Pmax), circuit power (Pc)
%% ouput: global energy efficiency (alpha)


function alpha = optimal_solution_NOMA(G,A,Pmax,Pmin,Pc)

%---- initialization ----%
epsilon=1e-7;
alpha=1e-5;
R_P=1;


b=(A(1)-1)/G(1);
c=1/prod(A(2:length(G)))*(Pmax-Pmin+b*prod(A(2:length(G))));

%----- Dinkelbach's algorithm ----%
while(R_P>epsilon)
    a=1/(2*log(2)*alpha*prod(A(2:length(G))))-1/G(1);
    
    %compute the optimal power allocation
    p(1)=min(c,max(a,b));
    p(2)=(A(2)-1)*(1/G(2)+p(1));
    
    if length(G)>2
        for k=3:length(G)
            B=0;
            for i=2:k-1
                B=B+(A(i)-1)/G(i)*prod(A(i+1:k-1));
            end
            p(k)=(A(k)-1)*(1/G(k)+p(1)*prod(A(2:k-1))+B);
        end
    end
    
    %compute the sum rate
    for k=1:length(G)
        r(k)=1/2*log2(1+G(k)*p(k)/(1+G(k)*(sum(p(1:k-1)))));
    end
    
    R_P=sum(r)-alpha*(sum(p)+Pc);  
    alpha=sum(r)/(sum(p)+Pc);

end


end
