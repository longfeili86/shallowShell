% This is the default script that defines given functions:
% such as the external forcing terms of the W and phi equations: f.phi(x,y), f.w(x,y)
% and the exact solution (if known): exact.phi(x,y), exact.w(x,y)
% and the initial shell shape: w0(x,y).
% 
% we can write other .m script to define these function. To use a different
% funcDefFile, simply pass the file name to the runShell -funcDefFile=<filename>



% forcing term
fx=2*pi/(xb-xa);
fy=2*pi/(yb-ya);
lambda=1;
if(parameters.isLinear)
    lambda=0;
end
f.w=@(x,y) -lambda.*((pi.^2.*sin((pi.*(x-xa).*2.0)./(xa-xb)).^3.*sin((pi.*(y-ya).*2.0)./(ya-yb)).^3.*1.0./(xa-xb).^2.*1.2e1-pi.^2.*cos((pi.*(x-xa).*2.0)./(xa-xb)).^2.*sin((pi.*(x-xa).*2.0)./(xa-xb)).*sin((pi.*(y-ya).*2.0)./(ya-yb)).^3.*1.0./(xa-xb).^2.*2.4e1).*(pi.^2.*sin((pi.*(x-xa).*2.0)./(xa-xb)).^5.*sin((pi.*(y-ya).*2.0)./(ya-yb)).^5.*1.0./(ya-yb).^2.*2.0e1-pi.^2.*cos((pi.*(y-ya).*2.0)./(ya-yb)).^2.*sin((pi.*(x-xa).*2.0)./(xa-xb)).^5.*sin((pi.*(y-ya).*2.0)./(ya-yb)).^3.*1.0./(ya-yb).^2.*8.0e1)+(pi.^2.*sin((pi.*(x-xa).*2.0)./(xa-xb)).^3.*sin((pi.*(y-ya).*2.0)./(ya-yb)).^3.*1.0./(ya-yb).^2.*1.2e1-pi.^2.*cos((pi.*(y-ya).*2.0)./(ya-yb)).^2.*sin((pi.*(x-xa).*2.0)./(xa-xb)).^3.*sin((pi.*(y-ya).*2.0)./(ya-yb)).*1.0./(ya-yb).^2.*2.4e1).*(pi.^2.*sin((pi.*(x-xa).*2.0)./(xa-xb)).^5.*sin((pi.*(y-ya).*2.0)./(ya-yb)).^5.*1.0./(xa-xb).^2.*2.0e1-pi.^2.*cos((pi.*(x-xa).*2.0)./(xa-xb)).^2.*sin((pi.*(x-xa).*2.0)./(xa-xb)).^3.*sin((pi.*(y-ya).*2.0)./(ya-yb)).^5.*1.0./(xa-xb).^2.*8.0e1)-pi.^4.*cos((pi.*(x-xa).*2.0)./(xa-xb)).^2.*cos((pi.*(y-ya).*2.0)./(ya-yb)).^2.*sin((pi.*(x-xa).*2.0)./(xa-xb)).^6.*sin((pi.*(y-ya).*2.0)./(ya-yb)).^6.*1.0./(xa-xb).^2.*1.0./(ya-yb).^2.*7.2e3)+pi.^4.*sin((pi.*(x-xa).*2.0)./(xa-xb)).^5.*sin((pi.*(y-ya).*2.0)./(ya-yb)).^5.*1.0./(xa-xb).^4.*1.04e3+pi.^4.*sin((pi.*(x-xa).*2.0)./(xa-xb)).^5.*sin((pi.*(y-ya).*2.0)./(ya-yb)).^5.*1.0./(ya-yb).^4.*1.04e3+pi.^4.*cos((pi.*(x-xa).*2.0)./(xa-xb)).^4.*sin((pi.*(x-xa).*2.0)./(xa-xb)).*sin((pi.*(y-ya).*2.0)./(ya-yb)).^5.*1.0./(xa-xb).^4.*1.92e3+pi.^4.*cos((pi.*(y-ya).*2.0)./(ya-yb)).^4.*sin((pi.*(x-xa).*2.0)./(xa-xb)).^5.*sin((pi.*(y-ya).*2.0)./(ya-yb)).*1.0./(ya-yb).^4.*1.92e3+pi.^4.*sin((pi.*(x-xa).*2.0)./(xa-xb)).^5.*sin((pi.*(y-ya).*2.0)./(ya-yb)).^5.*1.0./(xa-xb).^2.*1.0./(ya-yb).^2.*8.0e2-pi.^4.*cos((pi.*(x-xa).*2.0)./(xa-xb)).^2.*sin((pi.*(x-xa).*2.0)./(xa-xb)).^3.*sin((pi.*(y-ya).*2.0)./(ya-yb)).^5.*1.0./(xa-xb).^4.*7.04e3-pi.^4.*cos((pi.*(y-ya).*2.0)./(ya-yb)).^2.*sin((pi.*(x-xa).*2.0)./(xa-xb)).^5.*sin((pi.*(y-ya).*2.0)./(ya-yb)).^3.*1.0./(ya-yb).^4.*7.04e3-pi.^4.*cos((pi.*(x-xa).*2.0)./(xa-xb)).^2.*sin((pi.*(x-xa).*2.0)./(xa-xb)).^3.*sin((pi.*(y-ya).*2.0)./(ya-yb)).^5.*1.0./(xa-xb).^2.*1.0./(ya-yb).^2.*3.2e3-pi.^4.*cos((pi.*(y-ya).*2.0)./(ya-yb)).^2.*sin((pi.*(x-xa).*2.0)./(xa-xb)).^5.*sin((pi.*(y-ya).*2.0)./(ya-yb)).^3.*1.0./(xa-xb).^2.*1.0./(ya-yb).^2.*3.2e3+pi.^4.*cos((pi.*(x-xa).*2.0)./(xa-xb)).^2.*cos((pi.*(y-ya).*2.0)./(ya-yb)).^2.*sin((pi.*(x-xa).*2.0)./(xa-xb)).^3.*sin((pi.*(y-ya).*2.0)./(ya-yb)).^3.*1.0./(xa-xb).^2.*1.0./(ya-yb).^2.*1.28e4;
f.phi=@(x,y) -lambda.*((pi.^2.*sin((pi.*(x-xa).*2.0)./(xa-xb)).^5.*sin((pi.*(y-ya).*2.0)./(ya-yb)).^5.*1.0./(xa-xb).^2.*2.0e1-pi.^2.*cos((pi.*(x-xa).*2.0)./(xa-xb)).^2.*sin((pi.*(x-xa).*2.0)./(xa-xb)).^3.*sin((pi.*(y-ya).*2.0)./(ya-yb)).^5.*1.0./(xa-xb).^2.*8.0e1).*(pi.^2.*sin((pi.*(x-xa).*2.0)./(xa-xb)).^5.*sin((pi.*(y-ya).*2.0)./(ya-yb)).^5.*1.0./(ya-yb).^2.*2.0e1-pi.^2.*cos((pi.*(y-ya).*2.0)./(ya-yb)).^2.*sin((pi.*(x-xa).*2.0)./(xa-xb)).^5.*sin((pi.*(y-ya).*2.0)./(ya-yb)).^3.*1.0./(ya-yb).^2.*8.0e1)-pi.^4.*cos((pi.*(x-xa).*2.0)./(xa-xb)).^2.*cos((pi.*(y-ya).*2.0)./(ya-yb)).^2.*sin((pi.*(x-xa).*2.0)./(xa-xb)).^8.*sin((pi.*(y-ya).*2.0)./(ya-yb)).^8.*1.0./(xa-xb).^2.*1.0./(ya-yb).^2.*1.0e4)-pi.^4.*sin((pi.*(x-xa).*2.0)./(xa-xb)).^3.*sin((pi.*(y-ya).*2.0)./(ya-yb)).^3.*1.0./(xa-xb).^4.*3.36e2-pi.^4.*sin((pi.*(x-xa).*2.0)./(xa-xb)).^3.*sin((pi.*(y-ya).*2.0)./(ya-yb)).^3.*1.0./(ya-yb).^4.*3.36e2+pi.^4.*cos((pi.*(x-xa).*2.0)./(xa-xb)).^2.*sin((pi.*(x-xa).*2.0)./(xa-xb)).*sin((pi.*(y-ya).*2.0)./(ya-yb)).^3.*1.0./(xa-xb).^4.*9.6e2+pi.^4.*cos((pi.*(y-ya).*2.0)./(ya-yb)).^2.*sin((pi.*(x-xa).*2.0)./(xa-xb)).^3.*sin((pi.*(y-ya).*2.0)./(ya-yb)).*1.0./(ya-yb).^4.*9.6e2-pi.^4.*sin((pi.*(x-xa).*2.0)./(xa-xb)).^3.*sin((pi.*(y-ya).*2.0)./(ya-yb)).^3.*1.0./(xa-xb).^2.*1.0./(ya-yb).^2.*2.88e2+pi.^4.*cos((pi.*(x-xa).*2.0)./(xa-xb)).^2.*sin((pi.*(x-xa).*2.0)./(xa-xb)).*sin((pi.*(y-ya).*2.0)./(ya-yb)).^3.*1.0./(xa-xb).^2.*1.0./(ya-yb).^2.*5.76e2+pi.^4.*cos((pi.*(y-ya).*2.0)./(ya-yb)).^2.*sin((pi.*(x-xa).*2.0)./(xa-xb)).^3.*sin((pi.*(y-ya).*2.0)./(ya-yb)).*1.0./(xa-xb).^2.*1.0./(ya-yb).^2.*5.76e2-pi.^4.*cos((pi.*(x-xa).*2.0)./(xa-xb)).^2.*cos((pi.*(y-ya).*2.0)./(ya-yb)).^2.*sin((pi.*(x-xa).*2.0)./(xa-xb)).*sin((pi.*(y-ya).*2.0)./(ya-yb)).*1.0./(xa-xb).^2.*1.0./(ya-yb).^2.*1.152e3;

%exact solution is known for this test
knownExactSolution=true;
exact.w=@(x,y) (sin(fx*(x-xa)).*sin(fy*(y-ya))).^5;
exact.phi=@(x,y) (sin(fx*(x-xa)).*sin(fy*(y-ya))).^3;

% initial shell shape: w0(x,y).
w0=@(x,y) 0.*x;