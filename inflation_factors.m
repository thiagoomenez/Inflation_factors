function [alpha, N, gama, v] = inflation_factors(a, Na, Ne, Nd, d, dobs, LCd)

dD = (1/sqrt(Ne - 1))*(d - mean(d, 2));

Gd = LCd*dD;

[U, S, V] = svd(Gd);
v = diag(S);
ng = rank(Gd);
v = v(1:ng);

y = LCd*(dobs - mean(d, 2));

switch a
    
    case 1 % ES-MDA-EQL
        alpha = ones(Na, 1)*Na;
        N = Na;
        gama = 1;
        
    case 2 % ES-MDA-GEO3
        anamax = Na;
        mu_alpha = 1.2;
        rho = 0.5;
        vm = mean(v);
        vmsq = (rho/(1 - rho))*vm^2;
        
        anam = mean(v(end))^2;
        p = 1;
        while (anam <= mu_alpha)
            anam = mean(v(end-p:end))^2;
            p = p+1;
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
        
    case 3 % ES-MDA-GEO1
        rho = 0.5;
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
        
        
    case 4 % ES-MDA-GEO2
        anam = 1.5; amax = 10e5; Namax = 100;
        
        g = @(x)(V*(S'*S + x*eye(size(S'*S)))*S'*U'*y);
        f = @(x)(norm(S*V'*g(x) - U'*y)^2);
        em_f = @(x)(f(x) - sqrt(Nd));
        
        alphaprime = bisection(em_f, 1e-5, 10e50);
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
        
    case 6 % ES-MDA-GEO2 (mu_alpha)
        amax = 10e5; Namax = 100; anamax = Na;
        
        anam = mean(v(end))^2;
        p = 1;
        while (anam <= 1.55)
            anam = mean(v(end-p:end))^2;
            p = p+1;
        end
        
        if(anam > anamax)
            anam = mu_alpha;
        end
        
        g = @(x)(V*(S'*S + x*eye(size(S'*S)))*S'*U'*y);
        f = @(x)(norm(S*V'*g(x) - U'*y)^2);
        em_f = @(x)(f(x) - sqrt(ng));
        
        alphaprime = bisection(em_f, 1e-5, 10e50);
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
