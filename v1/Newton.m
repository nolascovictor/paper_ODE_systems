function [convergiu, solaprox, num_itr] = Newton(f, d_uf, a, errtol, maxitr)
[m,n] = size(a);
%Newton for F_cosm
convergiu = 0;
num_itr = 0;

F_Bnm = @(a) integral(@(t) ([1 sqrt(2).*cos((1:(m-1)).*pi.*t)]')*([1 sqrt(2).*cos((1:(m-1)).*pi.*t)]*a - f(t, compute_u(a, t))), 0, 1,'ArrayValued',true);

erro = sqrt(sum(integral(@(t) ([1 sqrt(2).*cos((1:(m-1)).*pi.*t)]*F_Bnm(a)).^2, 0, 1,'ArrayValued',true)));

if erro < errtol
   convergiu = 1;
   solaprox = a;
   return
end

xold = a;
while (num_itr < maxitr)
   
    xnew = xold-reshape(compute_nu0T1(d_uf, xold)\reshape(F_Bnm(xold), [m*n,1]), [m,n]);
    
    erro = sqrt(sum(integral(@(t) ([1 sqrt(2).*cos((1:(m-1)).*pi.*t)]*F_Bnm(xnew)).^2, 0, 1,'ArrayValued',true)));
    num_itr = num_itr+1;
    if (erro < errtol)
        convergiu = 1; 
        solaprox = xnew;
        return
    end 
    xold = xnew;
end

solaprox = xnew;

end

%%%%%%%%%%%%%%%%%%%% Support functions%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% Support functions%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function YY = compute_nu0T1(d_uf, a)
% Está função calcula DF(a)^(-1) para uma sequecia 'a' de coeficientes de uma
% solução numérica. Note que isso representa calcular o lambda.
% A função F_Bnm abaixo representa a F (operador) para uma sequecia 'a' de coeficientes de uma
% solução numérica.
% F_Bnm = @(a) integral(@(t) ([1 sqrt(2).*cos((1:(m-1)).*pi.*t)]')*([1 sqrt(2).*cos((1:(m-1)).*pi.*t)]*a - new_f(x, compute_u(a, t))), 0, 1,'ArrayValued',true);

[m,n] = size(a);
% Compute A = DF_cosm(w)
w = @(x) compute_u(a, x);
YY = zeros(n*m);
Id_n = eye(n);
for k = 1:n
    for j = 1:m
        if j>=2
            uj = @(x) sqrt(2)*sin((j-1)*pi*x) / ((j-1)*pi);
            duj = @(x) sqrt(2)*cos((j-1)*pi*x);
        else
            uj = @(x) x;
            duj = @(x) 1;
        end
        YY(:,(k-1)*m+j) = reshape(integral(@(x)([1 sqrt(2).*cos((1:(m-1)).*pi.*x)]')*(duj(x).*Id_n(k,:) - (d_uf(x,w(x))*(uj(x).*Id_n(k,:))')'), 0, 1, 'ArrayValued',true), [1,m*n]);
    end
end       
end