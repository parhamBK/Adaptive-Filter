function [Wopt,W,error] = problem2_LMS(u, d, N, mu, W_init)

    arguments
        u
        d
        N
        mu
        W_init = zeros(N,1) ;  % default value
    end

    W0 = zeros(N,1) + W_init ;
    W_n = W0 ;
    M = length(W0) ;
    W = zeros(N,length(u) - M + 1) ;
    error = zeros(length(u),1) ;
    k = 0 ;
    for n = M:length(u)
        k = k + 1 ;
        x = u(n:-1:n-M+1); 
        y_n = W_n' * x ;
        e_n = d(n) - y_n ;
        W_n1 = W_n + mu * x * conj(e_n);
        W_n = W_n1 ;
        W(:,n-M+1) = W_n ;
        error(k,1) = e_n ;
    end

    Wopt = W_n ;

end