function [b, t0, L, x0] = compute_solution(f, d_uf, m, n, intvSol, initCond, finalCond)
t0 = intvSol(1); L = intvSol(2);
x0 = initCond;
if isempty(finalCond) %BVP or IVP?
    x0 = cell2mat(x0);
    new_f = @(t,x) (f((L - t0).*t + t0, x + x0').*(L - t0))';
    opts = odeset('RelTol',1e-12,'AbsTol',1e-20.*ones(1,n));
    xmesh = linspace(0,1);
    [x, y] = ode45(new_f, xmesh, zeros(n,1), opts);
else
    x1 = finalCond;
    index1 = zeros(1,n);
    index0 = zeros(1,n);
    for i=1:n % to know which index is filled with numbers
        if ~isempty(x1{i})
            index1(i) = i;
        end
        if ~isempty(x0{i})
            index0(i) = i;
        end
    end
    index1 = index1(index1>0);
    index0 = index0(index0>0);
    if isempty(index0)
        disp('Please fulfil initial condition!')
        return
    elseif (length(index1)+length(index0)>n)
        if (length(index1)==n)&&(length(index0)==n)
            x0 = cell2mat(x0); x1 = cell2mat(x1);
            new_f = @(t,x) (f((L - t0).*t + t0, x + x0').*(L - t0))';
            init_guess = new_f(0,zeros(n,1));
            opt = optimoptions('fsolve','Algorithm','levenberg-marquardt','Display','off');
            rng default
            bvp_guess = fsolve(@(s) residual(new_f, s, n, x0, x1, 0, 1),init_guess,opt);%guess to zero u2_0
            options = odeset('RelTol',1e-12,'AbsTol',1e-20.*ones(1,n));
            xmesh = linspace(0,1);
            [t,y] = ode45(new_f,  xmesh, bvp_guess,options);
            C = zeros(length(t),m);
            for i=1:length(t)
                C(i,:) = [t(i) sqrt(2).*sin((1:(m-1)).*pi.*t(i))./((1:(m-1)).*pi)];
            end
            b = ((C'*C)\C')*y;
            guess_w = @(x) compute_u(b,x)';
            solinit = bvpinit(xmesh,guess_w);
            Mod_bc_eq = @(yt0,yL) bc_eq(yt0,yL,x0,x1);
            sol = bvp4c(new_f, Mod_bc_eq, solinit);
            x = sol.x; y = sol.y';
        else
            disp('Length of initial condition plus length of final condition must be 2n(each one=n) or equal to n!')
            return
        end
    else
        v0 = zeros(1,n); v1 = zeros(1,n);
        v0(index0(1:end)) = x0{index0(1:end)};
        v1(index1(1:end)) = x1{index1(1:end)};
        new_f = @(t,x) f(t,x)';
        init_guess = v0';
        opt = optimoptions('fsolve','Algorithm','levenberg-marquardt','Display','off');
        rng default
        xmesh = linspace(t0,L);
        bvp_guess = fsolve(@(s) residual(new_f, s, n, zeros(1,n), v1, t0, L),init_guess,opt);%guess to zero u2_0
        options = odeset('RelTol',1e-12,'AbsTol',1e-20.*ones(1,n));
        [~,y] = ode45(new_f,  xmesh, bvp_guess,options);
        Mod_bc_eq = @(yt0,yL) bc_eq2(yt0,yL,x0,x1,index0,index1);
        solinit = bvpinit(xmesh,y(1,:));
        sol = bvp4c(new_f, Mod_bc_eq, solinit);
        y = sol.y';
        X0 = y(1,:); X1 = y(end,:);
        X0(index0(1:end)) = x0{index0(1:end)}; X1(index1(1:end)) = x1{index1(1:end)};
        new_f = @(t,x) (f((L - t0).*t + t0, x + X0').*(L - t0))';
        init_guess = new_f(0,zeros(n,1));
        opt = optimoptions('fsolve','Algorithm','levenberg-marquardt','Display','off');
        rng default
        bvp_guess = fsolve(@(s) residual(new_f, s, n, X0, X1, 0, 1),init_guess,opt);%guess to zero u2_0
        options = odeset('RelTol',1e-12,'AbsTol',1e-20.*ones(1,n));
        xmesh = linspace(0,1);
        [t,y] = ode45(new_f,  xmesh, bvp_guess,options);
        C = zeros(length(t),m);
        for i=1:length(t)
            C(i,:) = [t(i) sqrt(2).*sin((1:(m-1)).*pi.*t(i))./((1:(m-1)).*pi)];
        end
        b = ((C'*C)\C')*y;
        guess_w = @(x) compute_u(b,x)';
        solinit = bvpinit(xmesh,guess_w);
        Mod_bc_eq = @(yt0,yL) bc_eq3(yt0,yL,X0,X1);
        sol = bvp4c(new_f, Mod_bc_eq, solinit);
        x = sol.x; y = sol.y';
        x0 = X0;
    end
end
C = zeros(length(x),m);
for i=1:length(x)
    C(i,:) = [x(i) sqrt(2).*sin((1:(m-1)).*pi.*x(i))./((1:(m-1)).*pi)];
end
b = ((C'*C)\C')*y;
b(abs(b) < 1e-14) = 0;
new_f = @(t,x) (f((L - t0).*t + t0, x + x0).*(L - t0));
new_d_uf = @(t,x) (L-t0).*d_uf((L - t0)*t + t0, x + x0);
Eta = @(x) sqrt(sum(x.^2) - 2*(x(1:m:end)*integral(@(t) new_f(t, compute_u(reshape(x, [m,n]), t)), 0, 1,'ArrayValued',true)') + sum(integral(@(t) new_f(t, compute_u(reshape(x, [m,n]), t)).^2, 0, 1,'ArrayValued',true)) - 2*sqrt(2)*sum(integral(@(t) (((reshape(x, [m,n])')*([0 cos((1:(m-1)).*(pi*t))]'))')*(new_f(t, compute_u(reshape(x, [m,n]), t))'), 0, 1,'ArrayValued',true)));
options = optimoptions('fminunc', 'Algorithm', 'quasi-newton', 'HessUpdate', 'bfgs', 'StepTolerance', 1e-16, 'OptimalityTolerance', 1e-16); %  'CheckGradients', 'on'
options = optimoptions(options,'UseParallel',true);
theta = fminunc(Eta, reshape(b, [1, m*n]), options);
b = reshape(theta, [m, n]);
[convergiu, solaprox, num_itr] = Newton(new_f, new_d_uf, b, 1e-6, 100);
b = solaprox;
end

% F_Bnm = @(x) integral(@(t) [1 sqrt(2).*cos((1:(m-1)).*pi.*x)]*(([1 sqrt(2).*cos((1:(m-1)).*pi.*t)]')*(compute_du(b, t) - f(x, compute_u(b, t)))), 0, 1,'ArrayValued',true);

%% Support functions

function res = bc_eq(yt0,yL,initCond,finalCond)%
res = [yt0(1);yL(1)-(finalCond(1)-initCond(1))];
end

function res = bc_eq2(yt0,yL,x0,x1,index0,index1)
v0(index0(1:end)) = x0{index0(1:end)};
v1(index1(1:end)) = x1{index1(1:end)};
res = [yt0(index0(1:end))'-v0(index0(1:end)) yL(index1(1:end))'-v1(index1(1:end))]';
end

function res = bc_eq3(yt0,yL,initCond,finalCond)%
v1 = finalCond(1:(end-1));
v0 = initCond(1:(end-1));
res = [yt0(1) yL(1:(end-1))'-(v1-v0)]';
end

function res = residual(f, s, n, initCond, finalCond, t0, L)
options = odeset('RelTol',1e-12,'AbsTol',1e-20.*ones(1,n));
x = linspace(t0,L);
sol = ode45(f,x,s, options);
u = sol.y(:,end)';
res = u-(finalCond-initCond);
end

function F_a = F_Bnm(f, a) % coefficients of F(w) in L^2(I)^n cosine base
[m,n] = size(a); % truncations size m, and systems size n
F_u = @(x) compute_du(a, x) - f(x, compute_u(a, x));
F_a = zeros(m,n);
F_a(1,:) = integral(F_u, 0, 1,'ArrayValued',true);
for i = 2:m
    F_a(i,:) = integral(@(x) F_u(x).*sqrt(2).*cos((i-1).*pi.*x), 0, 1, 'ArrayValued',true);
end
end

function du = compute_du(a, x)
[m, n] = size(a); % size of truncation, a is mxn
k = size(x,2); % size of input, x is 1xk
if ~isa(x,'taylor')
    if m>1 % size of truncation greater than 1
        if k==1 % x is a number
            dbaseH1 = [1 sqrt(2).*cos((1:(m-1)).*pi.*x)];  % derivative of H^{1,0}_B base
            du = dbaseH1*a; % 1xn vector function
        else % x is a vector
            A = zeros(k,m);% each row i of A is the base on x(i)
            for i=1:k
                A(i,:) = [1 sqrt(2).*cos((1:(m-1)).*pi.*x(i))];
            end
            du = reshape(A*a, [1,k*size(a,2)]); % 1x(kxn), [A(x(1))*a(:,1),...,A(x(k))*a(:,1),...,A(x(1))*a(:,n),......,A(x(k))*a(:,n)]
        end
    else % size of truncation equal 1, that is, a is row vector
        du = zeros(1,k*size(a,2)); %1x(kxn) [a(1,1),...,a(1,1),...,(1,n),......,(1,n)]
        for i =1:size(a,2)
            du(((i-1)*k+1):i*k) = a(1,i);
        end
    end
else
    v = x(1).t;
    taylorIndex = size(v,2)-1;
    f_aux = @(y) y;
    du = f_aux(taylorinit(1:n, taylorIndex));
    for j = 1:n
        du(1,j) = a(1,j);
        for i = 2:m
            du(1,j) = du(1,j) + sqrt(2).*cos((i-1).*pi.*x)*a(i,j);
        end
    end
end
end