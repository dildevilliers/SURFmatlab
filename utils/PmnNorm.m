function P = PmnNorm(x,m,n)

% Computes the normalised Legendre function of the first kind 
% P^m_n(x) of integer order and degree m and n
% Not tested for x dimensions higher than 1...
% Returns a column vector of length(x)

Pmat = legendre(n,x,'norm');
P = Pmat(m+1,:);
P = P(:);