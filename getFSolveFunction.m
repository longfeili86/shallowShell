function F=getFSolveFunction(x,Aphi,Aw,RHSphi,RHSw,bcType,R)
% build F for fsolve
%
% -- Longfei Li

% num of phi eqn. It is never singular, so n is also number of nodes
n=length(Aphi(:,1)); 

phi=x(1:n);
w=x(n+1:2*n);

F(1:n)= Aphi*phi-RHSphi(w,phi); % phi eqn
if(bcType==3) % w eqn is singular, we need extra equations for lambda
    lambda=x(2*n+1:2*n+3); % 3 extra unknowns
    F(n+1:2*n+3) = Aw*[w;lambda]-[RHSw(w,phi);R]; % w + lambda eqn
else
    F(n+1:2*n) = Aw*w-RHSw(w,phi); % w eqn
end

end