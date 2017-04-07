function RHS=assignBoundaryConditionsRHS(RHS,Index,parameter)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this function assigns the rhs for corresponding boundary conditions
% Input:
%   RHS: the RHS vector of the system
%   Index: index object
%   parameter: parameter object
% Output:
% Output:
%   RHS: the RHS vector with BC implemented
% 
% by L. Li 07/01/2015
% modified by H. Ji 12/25/2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%simplify names
bcType = parameter.bcType;
mixedbcType = parameter.mixedbcType;
nx = Index.nx;
ny = Index.ny;

    if(bcType==1) % simply supported
        RHS(Index.Boundary) = 0;
        RHS(Index.Corners) = 0;

    elseif(bcType==2) % clamped edge 
        RHS(Index.Boundary) = 0;
        RHS(Index.Corners) = 0;
        
    elseif(bcType==3) % free edge 
         % RHS should be the forcing itself.
    end
    

     if(mixedbcType==1) % simply supported partial edge
        RHS(Index.mixedBoundary) = 0;
        switch bcType % address singularity where BCs changes by analytical solution              
            case 1 % simply supported bc, do nothing

            case 2 % clamped bc, correct
                  for i = 1:8
                    x0y0 = Index.mixedBoundaryEnd(i,1);

                    if x0y0>0
                        x0yp = x0y0-1;           
                        xmy0 = x0y0-ny;
                        xpy0 = x0y0+ny;           
                        x0ym = x0y0+1;

                        switch Index.mixedBoundaryEnd(i,2)
                        case 1 % mixedB1
                            RHS(x0yp)=0;
                        case 2 % mixedB2
                            RHS(x0yp)=0;
                        case 3 % mixedT1
                            RHS(x0ym)=0;
                        case 4 % mixedT2
                            RHS(x0ym)=0;
                        case 5 % mixedR1
                            RHS(xmy0)=0;
                        case 6 % mixedR2
                            RHS(xmy0)=0;
                        case 7 % mixedL1
                            RHS(xpy0)=0;
                        case 8 % mixedL2
                            RHS(xpy0)=0;
                        otherwise
                            disp('error in constructing the matrix!')
                            return;
                        end
                    end
                 end
                

            case 3 % free bc

                
            otherwise

        end
        
        
        elseif(mixedbcType==2) % clamped partial edge
             RHS(Index.mixedBoundary) = 0;
            switch bcType % address singularity where BCs changes by analytical solution              
                case 1 % simply supported bc, correct
                   for i = 1:8
                        x0y0 = Index.mixedBoundaryEnd(i,1);

                        if x0y0>0
                            x0yp = x0y0-1;           
                            xmy0 = x0y0-ny;
                            xpy0 = x0y0+ny;           
                            x0ym = x0y0+1;

                            switch Index.mixedBoundaryEnd(i,2)
                            case 1 % mixedB1
                                RHS(x0yp)=0;
                            case 2 % mixedB2
                                RHS(x0yp)=0;
                            case 3 % mixedT1
                                RHS(x0ym)=0;
                            case 4 % mixedT2
                                RHS(x0ym)=0;
                            case 5 % mixedR1
                                RHS(xmy0)=0;
                            case 6 % mixedR2
                                RHS(xmy0)=0;
                            case 7 % mixedL1
                                RHS(xpy0)=0;
                            case 8 % mixedL2
                                RHS(xpy0)=0;
                            otherwise
                                disp('error in constructing the matrix!')
                                return;
                            end
                        end
                    end

                case 2 % clamped bc, do nothing

                    return;
                case 3 % free bc, correct
                    for i = 1:8                       
                        x0y0 = Index.mixedBoundaryEnd(i,1);
                        if x0y0>0
                            x0yp = x0y0-1;
                            xmy0 = x0y0-ny;
                            xpy0 = x0y0+ny;           
                            x0ym = x0y0+1;
                            switch Index.mixedBoundaryEnd(i,2)
                                case 1 % mixedB1
                                    RHS(x0yp)=0;
                                case 2 % mixedB2
                                    RHS(x0yp)=0;
                                case 3 % mixedT1
                                    RHS(x0ym)=0;                                    
                                case 4 % mixedT2
                                    RHS(x0ym)=0;
                                case 5 % mixedR1
                                    RHS(xmy0)=0;
                                case 6 % mixedR2
                                    RHS(xmy0)=0;
                                case 7 % mixedL1
                                    RHS(xpy0)=0;
                                case 8 % mixedL2
                                    RHS(xpy0)=0;
                                otherwise
                                    disp('error in constructing the matrix!')
                            end
                        end
                    end
                    
                otherwise

            end
    end
    
end