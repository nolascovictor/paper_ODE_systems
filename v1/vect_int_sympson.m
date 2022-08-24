function X = vect_int_sympson(f, a, b, n)
  m = size(f(a),2);  
  D = intval(b) - intval(a); H = D / n;
  x = a + intval(0:n) / n*D;
  w = 2 * ones(1,n+1); w(2:2:n) = 4; w(1) = 1; w(n+1) = 1;
  V = H/3 .* (w * reshape(f(intval(x)),[(n+1),m]));
  Y = f(taylorinit(infsup(a,b), 4)); % Approximation inclusion mx1
  E = (H^4*D / intval('7.5')) .* (Y{4}'); % Error term
  X = V- E; % Integral inclusion
end