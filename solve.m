function solve(parameters)
% This function do the actual solve for various case with given parameters
% --Longfei Li

infoPrefix = '--solve--: '; % all info displayed by this function includes this prefix

% simplify names
caseName=parameters.caseName;


switch (caseName)
    case 'biharmonic'
        fprintf('%sSolve the biharmonic equation\n',infoPrefix);
        biharmonic(parameters);
    case 'coupledSystem'
        fprintf('%sSolve the the coupled system\n',infoPrefix);
        coupledSystem(parameters);
    otherwise
        fprintf('%sSupported cases: biharmonic,exPicard,imPicard,newton.\n',infoPrefix);
        fprintf('%sTerminated.\n',infoPrefix);
        return
end







end