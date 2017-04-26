function [A,Q]=removeMatrixSingularity(A,myGrid,Index,quiet)
% remove singularity for Free BC
% --Longfei Li
if(nargin==3)
    quiet=false;
end

infoPrefix = '--removeMatrixSingularity--: '; % all info displayed by this function includes this prefix
Xvec = myGrid.XX(:);%column vector
Yvec = myGrid.YY(:);%column vector
a = zeros(length(Xvec),1); b=a; c=a;
a(Index.UsedPoints) = Xvec(Index.UsedPoints);
b(Index.UsedPoints) = Yvec(Index.UsedPoints);
c(Index.UsedPoints) = ones(length(Index.UsedPoints),1);  
Q=([a,b,c]); % right kernal of Aused
A= [A,Q;
    Q',zeros(3)]; % sum of [(xW),(yW),(W)]

if(~quiet)
    fprintf('%sMatrix is augmented to remove sigularity for free bc.\n',infoPrefix);
end

end