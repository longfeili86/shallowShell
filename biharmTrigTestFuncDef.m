% This is the script that defines given functions for biharmTrigTest:
% 
% use this file: runShell -funcDefFile=biharmTrigTestFuncDef


fx=2*pi/(xb-xa);
fy=2*pi/(yb-ya);

f.w=@(x,y) pi.^4.*cos((pi.*(x-xa).*2.0)./(xa-xb)).^4.*sin((pi.*(y-ya).*2.0)./(ya-yb)).^4.*1.0./(xa-xb).^4.*3.84e2+pi.^4.*cos((pi.*(y-ya).*2.0)./(ya-yb)).^4.*sin((pi.*(x-xa).*2.0)./(xa-xb)).^4.*1.0./(ya-yb).^4.*3.84e2+pi.^4.*sin((pi.*(x-xa).*2.0)./(xa-xb)).^4.*sin((pi.*(y-ya).*2.0)./(ya-yb)).^4.*1.0./(xa-xb).^4.*6.4e2+pi.^4.*sin((pi.*(x-xa).*2.0)./(xa-xb)).^4.*sin((pi.*(y-ya).*2.0)./(ya-yb)).^4.*1.0./(ya-yb).^4.*6.4e2+pi.^4.*sin((pi.*(x-xa).*2.0)./(xa-xb)).^4.*sin((pi.*(y-ya).*2.0)./(ya-yb)).^4.*1.0./(xa-xb).^2.*1.0./(ya-yb).^2.*5.12e2-pi.^4.*cos((pi.*(x-xa).*2.0)./(xa-xb)).^2.*sin((pi.*(x-xa).*2.0)./(xa-xb)).^2.*sin((pi.*(y-ya).*2.0)./(ya-yb)).^4.*1.0./(xa-xb).^4.*3.072e3-pi.^4.*cos((pi.*(y-ya).*2.0)./(ya-yb)).^2.*sin((pi.*(x-xa).*2.0)./(xa-xb)).^4.*sin((pi.*(y-ya).*2.0)./(ya-yb)).^2.*1.0./(ya-yb).^4.*3.072e3-pi.^4.*cos((pi.*(x-xa).*2.0)./(xa-xb)).^2.*sin((pi.*(x-xa).*2.0)./(xa-xb)).^2.*sin((pi.*(y-ya).*2.0)./(ya-yb)).^4.*1.0./(xa-xb).^2.*1.0./(ya-yb).^2.*1.536e3-pi.^4.*cos((pi.*(y-ya).*2.0)./(ya-yb)).^2.*sin((pi.*(x-xa).*2.0)./(xa-xb)).^4.*sin((pi.*(y-ya).*2.0)./(ya-yb)).^2.*1.0./(xa-xb).^2.*1.0./(ya-yb).^2.*1.536e3+pi.^4.*cos((pi.*(x-xa).*2.0)./(xa-xb)).^2.*cos((pi.*(y-ya).*2.0)./(ya-yb)).^2.*sin((pi.*(x-xa).*2.0)./(xa-xb)).^2.*sin((pi.*(y-ya).*2.0)./(ya-yb)).^2.*1.0./(xa-xb).^2.*1.0./(ya-yb).^2.*4.608e3;
f.phi=@(x,y) 0.*x+0.*y+0.;


%exact solution is known for this test
knownExactSolution=true;
exact.w=@(x,y) (sin(fx*(x-xa)).*sin(fy*(y-ya))).^4;
