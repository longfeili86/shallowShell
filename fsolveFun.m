function [F,J]=fsolveFun(x,n,W0,Aphi,Aw,mtx,Index,RHSphi,RHSw,R,parameters,xOld,dx)
% the syntax of fsolve requires one function returns fval and J.
%
% -- Longfei Li


F=FEvaluation(x,n,Aphi,Aw,RHSphi,RHSw,R,Index,parameters,xOld,dx);
J=JEvaluation(x,n,W0,Aphi,Aw,mtx,Index,parameters,dx);



end