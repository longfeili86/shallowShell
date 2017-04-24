function [Aused,RHSused]=removeSingularity(Aused,RHSused,myGrid,Index,addRHS)
% remove singularity for Free BC
% --Longfei Li

    infoPrefix = '--removeSingularity--: '; % all info displayed by this function includes this prefix
    Xvec = myGrid.XX(:);%column vector
    Yvec = myGrid.YY(:);%column vector
    a = Xvec(Index.UsedPoints);
    b = Yvec(Index.UsedPoints);
    c = ones(length(Index.UsedPoints),1);  
    Q=([a,b,c]); % right kernal of Aused
    P=([a,b,c]); % sum of [(xW),(yW),(W)]
    R=P'*addRHS;
    fprintf('%sFree BC additional rhs: r1=%f;r2=%f;r3=%f\n',infoPrefix,R(1),R(2),R(3));
    RHSused=[RHSused;R];
    Aused= [Aused,Q;
               P',zeros(3)];


end