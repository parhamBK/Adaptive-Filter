function [Wopt,W,error] = VLMS(u, d, N, mu, mu_min, mu_max, p, method, W_init)

    arguments
        u
        d
        N % taps
        mu
        mu_min
        mu_max
        p = 2 ;
        method = 1 ;
        W_init = zeros(N,1) ;  % default value
    end

    W0 = zeros(N,1) + W_init ;
    W_n = W0 ;
    M = length(W0) ;
    W = zeros(N,length(u) - M + 1) ;
    W(:,1) = W0 ;
    error = zeros(length(u),1) ;
    g = zeros(M,2) ;
    mu_i = mu * ones(N,1) ;
    k = 2 ;
    for n = M:length(u)

        x = u(n:-1:n-M+1); 
        y_n = W_n' * x ;
        e_n = d(n) - y_n ;

        for i = 0:M-1

            g(i+1,2) = g(i+1,1) ;             % g(i+1 , 2 ) --> gi(n-1)
            g(i+1,1) = conj(e_n) * x(i+1) ;    % g(i+1 , 1 ) --> gi(n)


            switch method
            case 1
                mu_i(i+1) = mu_i(i+1) + p * sign(g(i+1,1)) * sign(g(i+1,2)) ;
            case 2
                mu_i(i+1) = mu_i(i+1) + p * g(i+1,1) * g(i+1,2) ;
            otherwise

                if  sign(g(i+1,1)) == sign(g(i+1,2))
                     mu_i(i+1) = p * mu_i(i+1);
                else
                     mu_i(i+1) = mu_i(i+1) / p ;
                end

            end


            if mu_i(i+1) > mu_max
                mu_i(i+1) = mu_max ;
            elseif mu_i(i+1) < mu_min
                mu_i(i+1) = mu_min ;
            end


            W(i+1,n-M + 2) = W(i+1,n-M + 1) + mu_i(i+1) * g(i+1,1) ;
        end

        W_n = W(:, k);
        k = k + 1 ;
        error(n,1) = e_n ;
    end

    Wopt =  W_n ;

end
