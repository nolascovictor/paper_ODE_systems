function [eta, nu, K, t_star, t_2star] = verify_solution(a0, f, t0, L, x0, d_uf, R)

%% Define new f
new_f = @(t,x) f((L - t0)*t + t0, x + x0).*(L - t0);
new_d_uf = @(t,x) (L-t0).*d_uf((L - t0)*t + t0, x + x0);

%% Compute eta 
eta = compute_eta(new_f, a0);

%% Compute nu
nu = compute_nu(new_d_uf, a0); 

%% Compute r and K
K = compute_K(new_f, a0);

%% Compute t*
Delta = nu^2 - 2 * eta * K;
if Delta > 0
    t_star = (nu - sqrt(Delta)) / K;
    if (inf(t_star) > 0) && (sup(t_star) < R)
        disp('Proof was successful!')
        disp(['eta = ', num2str(sup(eta))])
        disp(['nu = ', num2str(inf(nu))])
        disp(['K = ', num2str(sup(K))])
        disp(['t_star = ', num2str(sup(t_star))])
        t_2star = sup((nu + sqrt(Delta)) / K);
        disp(['t_2star = ', num2str(t_2star)])
    elseif (inf(t_star) < 0)
        disp('Proof Failed! Root out of bounds. Please, try to increase m length.')
        disp(['eta = ', num2str(sup(eta))])
        disp(['nu = ', num2str(inf(nu))])
        disp(['K = ', num2str(sup(K))])
    else
        disp('Proof Failed! Root out of bounds. Please, try to increase R length.')
        disp(['eta = ', num2str(sup(eta))])
        disp(['nu = ', num2str(inf(nu))])
        disp(['K = ', num2str(sup(K))])
    end
else
    disp('Proof Failed! No real roots.')
    disp(['eta = ', num2str(sup(eta))])
    disp(['nu = ', num2str(inf(nu))])
    disp(['K = ', num2str(sup(K))])
end
end
%%%%%%%%%%%%%%%%%%%% Support functions%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% Support functions%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function iu = compute_iu(a, x)
[m, n] = size(a); % size of truncation, a is mxn
ipi = intval('pi');
if ~isa(x,'taylor')
    k = size(x,2); % size of input, x is 1xk
    if m>1 % size of truncation greater than 1
        if k==1 % x is a number
            baseH1 = [intval(x) sqrt(intval(2)).*sin((1:(m-1)).*ipi.*intval(x))./((1:(m-1)).*ipi)];  %H^{1,0}_B base
            iu = baseH1*intval(a); % 1xn vector function
        else % x is a vector
            A = intval(0).*zeros(k,m);% each row i of A is the base on x(i)
            for i=1:k
                A(i,:) = [intval(x(i)) sqrt(intval(2)).*sin((1:(m-1)).*ipi.*intval(x(i)))./((1:(m-1)).*ipi)];
            end
            iu = reshape(A*intval(a), [1,k*n]); % 1x(kxn), [A(x(1))*a(:,1),...,A(x(k))*a(:,1),...,A(x(1))*a(:,n),......,A(x(k))*a(:,n)]
        end
    else  % size of truncation equal 1, that is, a is row vector
        iu = intval(0).*zeros(1,k*n); %1x(kxn) [a(1,1),...,a(1,1),...,a(1,n),......,a(1,n)]
        for i =1:n
            iu(((i-1)*k+1):i*k) = intval(x).*intval(a(1,i));
        end
    end
else
    v = x(1).t;
    taylorIndex = size(v,2)-1;
    f_aux = @(y) y;
    iu = f_aux(taylorinit(intval(1)*(1:n), taylorIndex));
    for j = 1:n
        iu(1,j) = intval(x)*intval(a(1,j));
        for i = 2:m
            iu(1,j) = iu(1,j) + (sqrt(intval(2))*sin((i-1)*ipi*intval(x))/((i-1)*ipi))*intval(a(i,j));
        end
    end
end
end

function idu = compute_idu(a, x)
[m, n] = size(a); % size of truncation, a is mxn
ipi = intval('pi');
if ~isa(x,'taylor')
    k = size(x,2); % size of input, x is 1xk
    if m>1 % size of truncation greater than 1
        if k==1 % x is a number
            idbaseH1 = [intval(1) sqrt(intval(2)).*cos((1:(m-1)).*ipi.*intval(x))];  % derivative of H^{1,0}_B base
            idu = idbaseH1*a; % 1xn vector function
        else % x is a vector
            A = intval(0).*zeros(k,m);% each row i of A is the base on x(i)
            for i=1:k
                A(i,:) = [intval(1) sqrt(intval(2)).*cos((1:(m-1)).*ipi.*intval(x(i)))];
            end
            idu = reshape(A*a, [1,k*n]); % 1x(kxn), [A(x(1))*a(:,1),...,A(x(k))*a(:,1),...,A(x(1))*a(:,n),......,A(x(k))*a(:,n)]
        end
    else % size of truncation equal 1, that is, a is row vector
        idu = intval(0).*zeros(1,k*n); %1x(kxn) [a(1,1),...,a(1,1),...,(1,n),......,(1,n)]
        for i =1:n
            idu(((i-1)*k+1):i*k) = intval(a(1,i));
        end
    end
else
    v = x(1).t;
    taylorIndex = size(v,2)-1;
    f_aux = @(y) y;
    idu = f_aux(taylorinit(intval(1)*(1:n), taylorIndex));
    for j = 1:n
        idu(1,j) = intval(a(1,j));
        for i = 2:m
            idu(1,j) = idu(1,j) + sqrt(intval(2))*cos((i-1)*ipi*intval(x))*intval(a(i,j));
        end
    end
end
end
%%%%%%%%%%%%% vect_int_sympson function %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function X = vect_int_sympson(f, a, b, n)
m = size(f(a),2);
D = intval(b) - intval(a); H = D / n;
x = a + (intval(0:n)/n)*D;
w = 2 * ones(1,n+1); w(2:2:n) = 4; w(1) = 1; w(n+1) = 1;
vect_imgf = intval(zeros(n+1,m));
for i=1:(n+1)
    vect_imgf(i,:) = intval(f(x(i)));
end
V = H/3 .* (w * vect_imgf);
E = (H^4*D / intval('180'));
X = V- E; % Integral inclusion
end
%%%%%%%%%%%%% eta function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function eta = compute_eta(f, a) 
w = @(x) compute_iu(a, x); 
dw = @(x) compute_idu(a, x);
F_u_2 = @(x) (dw(x) - f(x,w(x)));
eta = sqrt(sum(abs(vect_int_sympson(F_u_2, 0, 1, 1000))));%(L2)^n norm
eta = intval(sup(eta));
end
%%%%%%%%%%%%% nu function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function nu = compute_nu(d_uf, a)
m = size(a,1);
ipi = intval('pi');
N = compute_N(d_uf, a);
[~,nu0] =  compute_nu0(d_uf, a);  
L = inf(min(nu0, 1));
nu = intval(inf(L - N/(ipi*intval(m))));
end

function [XX,nu0] = compute_nu0(d_uf, a)
[m,n] = size(a);
ipi = intval('pi');
% Compute A = DF_cosm(w)
w = @(x) compute_iu(a, x);
XX = intval(0).*zeros(n*m);
Id_n = eye(n);
for j = 1:m
    if j>=2
        uj = @(x) sqrt(intval(2))*sin((j-1)*ipi*intval(x)) / ((j-1)*ipi);
        duj = @(x) sqrt(intval(2))*cos((j-1)*ipi*intval(x));
    else
        uj = @(x) intval(x);
        duj = @(x) intval(1);
    end
    for k = 1:n
        A = intval(0).*zeros(m,n);
        for i = 1:m
            if i > 1
                Aij_base = @(x) (duj(x).*Id_n(k,:) - (d_uf(x,w(x))*(uj(x).*Id_n(k,:))')').*(sqrt(intval(2))*cos((i-1)*ipi*intval(x)));
            else
                Aij_base = @(x) duj(x).*Id_n(k,:) - (d_uf(x,w(x))*(uj(x).*Id_n(k,:))')';
            end
            %A(i,:) = integral(@(x) inf(Aij_base(x)), 0, 1, 'ArrayValued',true);
            A(i,:) = vect_int_sympson(Aij_base, 0, 1, 30);
        end
        XX(:,(k-1)*m+j) = (reshape(A, [1,n*m]))';  %%%%%%%%% XX(:,1:m) 
    end
end
% Compute nu
YY = XX^(-1);
% Frobenius norm of YY
nu0 = intval(norm(YY,'fro')^(-1));
nu0 = intval(inf(nu0));
end

function N = compute_N(d_uf, a)
w = @(x) compute_iu(a, x);
d_uf_w = @(x) d_uf(x, w(x)); 
d_uf_wInit = d_uf_w(taylorinit(infsup(0,1), 1));
if ~isa(d_uf_wInit,'taylor')
    N_matrix = (intval(1) + sqrt(intval(2))).*sup(abs(d_uf_wInit));
else
    XX = d_uf_wInit.t(:,:,1);
    YY = d_uf_wInit.t(:,:,2);    
    N_matrix = sup(abs(XX)) + sqrt(intval(2)).*sup(sqrt((abs(XX).^2 + abs(YY).^2))); 
end
N = intval(sup(norm(N_matrix,'fro')));
end
%%%%%%%%%%%%% K function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function K = compute_K(f, a)
n = size (a,2);
w = @(x) compute_iu(a, x);
y = @(x) f(x, hessianinit(w(x)));
YY = y(infsup(0,1));
K = intval(0);
for i = 1:n
    K_fi = sup(norm(YY(i).hx,'fro')); % sup of Frobenius norm of Hessian_fi
    K = K + K_fi^2;
end
K = sup(sqrt(K));
K = intval(max(K, 0.05)); 
end