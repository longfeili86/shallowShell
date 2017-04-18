function solve(parameters)
% This function do the actual solve for various case with given parameters
% --Longfei Li

infoPrefix = '--solve--: '; % all info displayed by this function includes this prefix

% simplify names
caseName=parameters.caseName;


switch (caseName)
    case 'biharmonic'
        fprintf('%sTest the biharmonic solver\n',infoPrefix);
        biharmonic(parameters);
    case 'exPicard'
        fprintf('%sTest the explicit Picard iterations\n',infoPrefix);
        %explicitPicardIterations(parameters);
        fprintf('%sFinish me...\n',infoPrefix);
        return
    case 'imPicard'
        fprintf('%sTest the semi-implicit Picard iterations\n',infoPrefix);
        %implicitPicardIterations(parameters);
        fprintf('%sFinish me...\n',infoPrefix);
        return
    case 'newton'
        fprintf('%sTest the newton iterations (direct solve)\n',infoPrefix);
        %newtonsIterations(parameters);
        fprintf('%sFinish me...\n',infoPrefix);
        return
    otherwise
        fprintf('%sSupported cases: biharmonic,exPicard,imPicard,newton.\n',infoPrefix);
        fprintf('%sTerminated.\n',infoPrefix);
        return
end







end