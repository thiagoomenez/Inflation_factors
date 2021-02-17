function [alpha, N, gama, v] = inflation_factors(a, Na, Ne, Nd, d, dobs, Cd)

dD = (1/sqrt(Ne))*(d - mean(d, 2));
 
Gd = sqrt(Cd)\dD;
 
[~, S, ~] = svd(Gd);
v = diag(S);
ng = rank(Gd);
v = v(1:ng);
 
y = sqrt(Cd)\(dobs - mean(d, 2));
xa = @(x)( (Gd'/(Gd*Gd' + x*eye(Nd)))*y );
rho = 0.5;
tau = sqrt(Nd);
eta = 1;
 
f1 = @(x)( norm(Gd*xa(x) - y) - tau*eta );
f2 = @(x)( x*x*norm( (Gd*Gd' + x*eye(Nd))\y )^2 - rho*rho*norm(y)^2 );
fsum = @(x, a, Na)((1 - 1/(x^Na))/(1 - 1/x) - a);

switch a
    
    case 0 % ES-MDA-EQL
        alpha = ones(Na, 1)*Na;
        N = Na;
        gama = 1;
        
    case 1 % ES-MDA-GEO1
        
        vm = mean(v);
        vmsq = (rho/(1 - rho))*vm^2;
        f = @(x)(fsum(x,vmsq, Na));
        
        auxc1 = 0.001;
        auxc2 = 1;
        
        c1 = f(auxc1);
        c2 = f(auxc2);
        den = 2;
        
        while (c1*c2 > 0)
            auxc2 = auxc1;
            auxc1 = auxc1/den;
            c1 = f(auxc1);
        end
       
        gama = bisection(f, auxc1, auxc2);
        
        alpha = compute_alpha(gama, vmsq, Na);
        N = Na;
        
    case 2 % ES-MDA-GEO2
        anam = 1.5; amax = 10e5; Namax = 100;
        
        g = @(x)(Gd'*( (Gd*Gd' + x*eye(Nd))\y ));
        f = @(x)(norm(Gd*g(x) - y)^2);
        em_f = @(x)(f(x) - tau*tau);
        
        alphaprime = -1;
        xteste = 10;
        while( alphaprime < 0)
            alphaprime = fzero(em_f, xteste);
            xteste = xteste*1.2;
        end
        
        if(alphaprime > amax)
            alphaprime = amax;
        end
        
        a1 = alphaprime - 1;
        
        while (a1 < alphaprime)
            gam = @(x)(fsum_em_f(x, anam, Na));
            gama = bisection(gam, 1e-5, 1);
            a1 = anam*gama^(1 - Na);
            if (a1 < alphaprime)
                Na = Na + 1;
            end
            if (Na > Namax)
                Na = Namax;
                break;
            end
        end
        alpha = compute_alpha(gama, a1, Na);
        N = Na;
        
    case 3 % ES-MDA-GEO3
        anamax = Na;
        mu_alpha = 1.1;
        vm = mean(v);
        vmsq = (rho/(1 - rho))*vm^2;
        
        if (v(end) > 1)
            anam = mean(v(end))^2;
        else
            anam = mean(v(end))^2;
            p = 1;
            while (anam <= mu_alpha)
                anam = mean(v(end-p:end))^2;
                p = p+1;
            end
        end
        
        if(anam > anamax)
            anam = mu_alpha;
        end
        
        f = @(x)(descobrir_Na(x, anam) - Na);
        
        while(f(vmsq) > 0)
            Na = Na + 1;
            f = @(x)(descobrir_Na(x, anam) - Na);
        end
        
        a1 = bisection(f, vmsq, 10e50);
        gama = ( anam/a1 )^( 1/(Na - 1));
        alpha = compute_alpha(gama, a1, Na);
        N = Na;
        
        
    case 4 % ES-MDA-GEO2 (mu_alpha)
        Namax = 100; 
        anamax = Na;
        mu_alpha = 1.1;
        
        anam = mean(v(end))^2;
        p = 1;
        while (anam <= mu_alpha)
            anam = mean(v(end-p:end))^2;
            p = p+1;
        end
        
        if(anam > anamax)
            anam = mu_alpha;
        end
        
        g = @(x)(Gd'*( (Gd*Gd' + x*eye(Nd))\y ));
        f = @(x)(norm(Gd*g(x) - y)^2);
        em_f = @(x)(f(x) - tau*tau);
        
        alphaprime = -1;
        xteste = 10;
        while( alphaprime < 0)
            alphaprime = fzero(em_f, xteste);
            xteste = xteste*1.2;
        end
        
        a1 = alphaprime - 1;
        
        while (a1 < alphaprime)
            gam = @(x)(fsum_em_f(x, anam, Na));
            gama = bisection(gam, 1e-5, 1);
            a1 = anam*gama^(1 - Na);
            if (a1 < alphaprime)
                Na = Na + 1;
            end
            if (Na > Namax)
                Na = Namax;
                break;
            end
        end
        alpha = compute_alpha(gama, a1, Na);
        N = Na;
        
    case 5 % alpha_1 = Hanke
        
        a1 = 1;
        a2 = 100;
        
        while (f2(a1)*f2(a2) > 0)
            a2 = a2*2;
        end
        
        zerof2 = bisection(f2, a1, a2);
        
        gsum = @(x)(fsum(x, zerof2, Na));
        xteste = 0.9;
        gama = fzero(gsum, xteste);
        
        while (or(gama < 0, isnan(gama)))
            xteste = xteste*0.9;
            gama = fzero(gsum, xteste);
        end
        
        alpha = ones(Na, 1);
        alpha(1) = zerof2;
        
        for i = 2:Na
            alpha(i) = alpha(i-1)*gama;
        end
        
        N = Na;
        
    case 6 % alpha_1 = Discrepancy Principle
        
        a1 = 1;
        a2 = 100;
        
        while (f1(a1)*f1(a2) > 0)
            a2 = a2*2;
        end
        
        zerof1 = bisection(f1, a1, a2);
        
        gsum = @(x)(fsum(x, zerof1, Na));
        xteste = 0.9;
        gama = fzero(gsum, xteste);
        
        while (or(gama < 0, isnan(gama)))
            xteste = xteste*0.9;
            gama = fzero(gsum, xteste);
        end
        
        alpha = ones(Na, 1);
        alpha(1) = zerof1;
        
        for i = 2:Na
            alpha(i) = alpha(i-1)*gama;
        end
        
        N = Na;
        
end


end

function p = bisection(f,a,b)
if f(a)*f(b)>0
    disp('Wrong choice')
else
    p = (a + b)/2;
    err = abs(f(p));
    while err > 1e-10
        if f(a)*f(p)<0
            b = p;
        else
            a = p;
        end
        p = (a + b)/2;
        err = abs(f(p));
    end
end
end

function [alp] = compute_alpha(a, a1, Na)

alp = zeros(Na, 1);
for i = 1:Na
    alp(i) = fteste(a, a1, i);
end

end

function [fsum] = fsum(a, a1, Na)

alpha = compute_alpha(a, a1, Na);

aux = sum(1./alpha);

fsum = aux - 1;

end

function [fsum_m] = fsum_em_f(gama, ana, Na)

fsum_m = gama.^Na - ana.*gama + ana - 1;

end

function [teste] = fteste(a, a1, k)

if(k == 1)
    teste = a1;
else
    teste = a.*(fteste(a, a1, k-1));
    
end
end

function [Na] = descobrir_Na(a1, ana)

a = (a1 - a1.*ana)./(ana - a1.*ana);
b = ana./a1;

Na = log(b)./log(a) + 1;

end
