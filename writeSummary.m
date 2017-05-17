function writeSummary(fval,exitflag,output)


algorithm=output.algorithm;
iterations=output.iterations;
message=output.message;

fprintf('----------------------Iteration Summary----------------------\n'); 
fprintf('algorithm: %s\n',algorithm);

if(exitflag==0)
    fprintf('Iteration does not converge after %i steps (maxIter reached)\n',iterations);
else
    fprintf('Iteration converges after %i steps\n',iterations); 
end
fprintf('L2 norm of fval=%e, Lmax norm of fval=%s\n',sqrt(sum(fval.^2)/length(fval)),max(abs(fval)));
fprintf('%s\n',message);
fprintf('-------------------------------------------------------------\n');

end
