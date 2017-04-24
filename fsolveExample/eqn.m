function F = eqn(x,a)

%F(1) = exp(-exp(-(x(1)+x(2)))) - x(2)*(1+x(1)^2);
%F(2) = x(1)*cos(x(2)) + x(2)*sin(x(1)) - 0.5;


F = x.*x-exp(x)*a;
F(end)=1-x(end);

end