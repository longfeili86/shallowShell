function RHS=assignBoundaryConditionsRHSHalfdomain(RHS,Index,parameter)
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
% modified by H. Ji 4/24/2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% For x>= a/2, clamped BC w = w_n = 0
%
% For x< a/2, Simply Supported BC w = w_nn = 0
%               or free boudary condition
% 


%simplify names
LeftbcType = parameter.LeftbcType;
RightbcType = parameter.RightbcType;

    if(LeftbcType==1) % simply supported
        RHS(Index.Boundary) = 0;
        RHS(Index.Corners) = 0;

    elseif(LeftbcType==2) % clamped edge 
        RHS(Index.Boundary) = 0;
        RHS(Index.Corners) = 0;
        
    elseif(LeftbcType==3) % free edge 
        % RHS should be the forcing itself.
    end
    

     if(RightbcType==1) % simply supported partial edge
        RHS(Index.RboundaryAll) = 0;
        switch LeftbcType % address singularity where BCs changes by analytical solution              
            case 1 % simply supported bc, do nothing

            case 2 % clamped bc, correct
                x0y0 = Index.Tsingular;
                x0ym = x0y0+1;
                RHS(x0ym)=0;
                
                x0y0 = Index.Bsingular;
                x0yp = x0y0-1;              
                RHS(x0yp)=0;
  
            case 3 % free bc, do nothing
                
            otherwise

        end
        
        
        elseif(RightbcType==2) % clamped partial edge
             RHS(Index.RboundaryAll) = 0;
            switch LeftbcType % address singularity where BCs changes by analytical solution              
                case 1 % simply supported bc, correct
                    x0y0 = Index.Tsingular;
                    x0ym = x0y0+1;

                    RHS(x0ym)=0;
                    

                    x0y0 = Index.Bsingular;
                    x0yp = x0y0-1;

                    RHS(x0yp)=0;

                case 2 % clamped bc, do nothing

                case 3 % free bc
                    x0y0 = Index.Tsingular;
                    x0ym = x0y0+1;
                    RHS(x0ym)=0;

                    x0y0 = Index.Bsingular;
                    x0yp = x0y0-1;                
                    RHS(x0yp)=0;
% 
%                     x0y0 = Index.TsingularAdjacent;          
%                     RHS(x0y0) = 0;
% 
%                     x0y0 = Index.BsingularAdjacent;          
%                     RHS(x0y0) = 0;

                otherwise

            end

    end
    
    


end