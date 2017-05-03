function Q=getKernalOfSingularMatrix(myGrid,Index)
% get the kernal/null of the singular matrix associated with bcType==3(free)
% its kernal is [x,y,1]
%
% -- Longfei Li

Xvec = myGrid.XX(:);%column vector
Yvec = myGrid.YY(:);%column vector
a = zeros(length(Xvec),1); b=a; c=a;
a(Index.UsedPoints) = Xvec(Index.UsedPoints);
b(Index.UsedPoints) = Yvec(Index.UsedPoints);
c(Index.UsedPoints) = ones(length(Index.UsedPoints),1);  

Q=([a,b,c]); % right kernal 

end