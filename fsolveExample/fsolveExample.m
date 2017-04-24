% we can specify Algorithm to use by setting
%('Algorithm','trust-region-dogleg');
% available algorithms are: 'trust-region-dogleg' (default), 'trust-region-reflective', or 'levenberg-marquardt'.

problem.options = optimoptions('fsolve','Display','iter','Algorithm', 'trust-region-dogleg');
problem.objective = @(x) [eqn(x(1:3),10),eqn2(x(4:5),2)];
problem.x0 = [0,0,0,0,0];
problem.solver = 'fsolve';


[x,fval,exitflag,output] = fsolve(problem);
