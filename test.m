clear all;
more off;
format long;

global L = 100;
global sigma = 0.2;
K = 0.01;

alpha = 0;
beta = 1;

M = 2000;
global h = 2*L/M;


function xi = xi(i)
  global L;
  global h;
  xi = -L + (i-1)*h;    
end

function v = v(alpha, x)
  global L;
  v = alpha*(x + L);
end

function w = w(beta, T)
  global sigma;
  w = (T - sigma) * (1 - T) * exp(-beta/T);
end

  t0 = clock();

  b = zeros(M+1, 1);

  Told = zeros(M+1, 1);
  for i = 1:M+1
    Told(i) = 1 + (sigma - 1) * (xi(i) + L) / (2*L);
  end

  l = u = d = zeros(1, M+1)';

  d(1) = 1; u(2) = 0;
  d(M+1) = 1; l(M) = 0;
  b(1) = 1;
  b(M+1) = sigma;
  for i = 2:M
    l(i-1) = 1/h^2 - v(alpha, xi(i))/(2*h);
    u(i+1) = 1/h^2 + v(alpha, xi(i))/(2*h);
    d(i) = -2/h^2;
    b(i) = -K*w(beta, Told(i));
  end
  A = spdiags([l d u], -1:1, M+1, M+1);

  elapsed_time = etime(clock(), t0);
  disp("matrix construction time");
  elapsed_time

  t0 = clock();
  T = A\b;
  elapsed_time = etime(clock(), t0);
  disp("system solving time");
  elapsed_time
