function [Wopt,W,error] = RLS(u, d, N, lambda, W_init, P_init)

    arguments
        u
        d
        N
        lambda
        W_init = zeros(N,1) ;  % default value
        P_init =  (0.00000000001)^-1 * eye(N) ;  % default value
    end
    P = P_init ;
    W0 = zeros(N,1) + W_init ;
    W_n = W0 ;
    M = length(W0) ;
    W = zeros(N,length(u) - M + 1) ;
    error = zeros(length(u),1) ;
    
    for n = M:length(u)

        x = u(n:-1:n-M+1); 
        y_n = W_n' * x ;
        k = ( (lambda^-1) * P * x ) / (1 + (lambda^-1) * x' * P * x) ;
        zeta = d(n) - y_n ;
        
        W_n1 = W_n + k * conj(zeta);
        P = (lambda^-1) * P - (lambda^-1)* k * x' * P ;
        W_n = W_n1 ;
        
        W(:,n-M+1) = W_n ;
        error(n,1) = zeta ;
    end

    Wopt = W_n ;

end