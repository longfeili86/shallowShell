function A=removeMatrixSingularity(A,myGrid,Index,quiet)
% remove singularity for Free BC
% --Longfei Li
if(nargin==3)
    quiet=false;
end

infoPrefix = '--removeMatrixSingularity--: '; % all info displayed by this function includes this prefix

Q=getKernalOfSingularMatrix(myGrid,Index);
A= [A,Q;
    Q',zeros(3)]; 

if(~quiet)
    fprintf('%sMatrix is augmented to remove sigularity for free bc.\n',infoPrefix);
end

end