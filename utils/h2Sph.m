function h = h2Sph(x,n)

% Function to compute the second kind spherical Hankel function used in the
% SWE: h^(2)_n(x).  Order n, and argument x

% Direct using ordinary Hankel functions
h = sqrt(pi./(2.*x)).*besselh(n+0.5,2,x);

% % Recursion
% h0 = 1i.*exp(-1i.*x)./x;
% h1 = (-1 + 1i./x).*h0./1i;
% if n == 0
%     hn = h0;
% elseif n == 1
%     hn = h1;
% else
%     hn = h1;
%     hnm1 = h0;
%     for nn = 1:n-1
%         hnp1 = (2*nn+1).*hn./x - hnm1;
%         hnm1 = hn;
%         hn = hnp1;
%     end
%     hn = hnp1;    
% end
% h = hn;

% % Recursion with attempt to sort out small argument issues
% xh0 = 1i.*exp(-1i.*x);
% xh1 = (1i.*x + 1).*xh0./x;
% if n == 0
%     xhn = xh0;
% elseif n == 1
%     xhn = xh1;
% else
%     xhn = xh1;
%     xhnm1 = xh0;
%     for nn = 1:n-1
%         xhnp1 = (2*nn+1).*xhn./x - xhnm1;
%         xhnm1 = xhn;
%         xhn = xhnp1;
%     end
%     xhn = xhnp1;    
% end
% h = xhn./x;

%% Farfield asymptote
h(isinf(x)) = 1i^(n+1); % Suppress the e^(-jkr)/kr factor

