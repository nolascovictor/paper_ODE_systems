function u = compute_u(a, x)
[m, n] = size(a); % size of truncation, a is mxn
k = size(x,2); % size of input, x is 1xk
if ~isa(x,'taylor')
    if m>1 % size of truncation greater than 1
        if k==1 % x is a number
            baseH1 = [x sqrt(2).*sin((1:(m-1)).*pi.*x)./((1:(m-1)).*pi)];  %H^{1,0}_B base
            u = baseH1*a; % 1xn vector function
        else % x is a vector
            A = zeros(k,m);% each row i of A is the base on x(i)
            for i=1:k
                A(i,:) = [x(i) sqrt(2).*sin((1:(m-1)).*pi.*x(i))./((1:(m-1)).*pi)];
            end
            u = reshape(A*a, [1,k*n]); % 1x(kxn), [A(x(1))*a(:,1),...,A(x(k))*a(:,1),...,A(x(1))*a(:,n),......,A(x(k))*a(:,n)]
        end
    else  % size of truncation equal 1, that is, a is row vector
        u = zeros(1,k*n); %1x(kxn) [a(1,1),...,a(1,1),...,a(1,n),......,a(1,n)]
        for i =1:n
            u(((i-1)*k+1):i*k) = x.*a(1,i);
        end
    end
else
    v = x(1).t;
    taylorIndex = size(v,2)-1;
    f_aux = @(y) y;
    u = f_aux(taylorinit(1:n, taylorIndex));
    for j = 1:n
        u(1,j) = x*a(1,j);
        for i = 2:m
            u(1,j) = u(1,j) + (sqrt(2)*sin((i-1)*pi*x)/((i-1)*pi))*a(i,j);
        end
    end
end
end