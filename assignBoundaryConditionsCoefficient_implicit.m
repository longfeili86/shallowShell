function A=assignBoundaryConditionsCoefficient_implicit(A,Index,parameter)
%assign boundary conditions for the coefficient matrix
%Input:
% A: differentiation matrix for biharmonic operator
% Index: 
% diffMtx: differentiation matrices
% bcType: boundary condition types: Supported=1, clamped=2, free=3


    %simplify names
    bcType = parameter.bcType;
    mixedbcType = parameter.mixedbcType;
    nu = parameter.nu;  
    I = speye(size(A));
    
    nx = Index.nx;
    ny = Index.ny;
    
    if (bcType==1) % simply supported w = w_nn=0
        
        %disp('simply supported bc'); 
        
        A(Index.Boundary,:) = I(Index.Boundary,:);   % w = 0 on boundary 
        A(Index.Corners,:) = I(Index.Corners,:);   % w = 0 on boundary 

         for i = 1:length(Index.InBoundaryCorner) 
            x0y0 = Index.InBoundaryCorner(i,1);

            switch Index.InBoundaryCorner(i,2) 
                case 1 % left boundary
                    A(x0y0,x0y0) = A(x0y0,x0y0)-1;               
                case 2 % right boundary
                    A(x0y0,x0y0) = A(x0y0,x0y0)-1;
                case 3 % bottom boundary
                    A(x0y0,x0y0) = A(x0y0,x0y0)-1;
                case 4 % top boundary
                    A(x0y0,x0y0) = A(x0y0,x0y0)-1; 
                case 5 % InCornerBL
                    A(x0y0,x0y0) = A(x0y0,x0y0)-2; 
                case 6 % InCornerBR
                    A(x0y0,x0y0) = A(x0y0,x0y0)-2; 
                case 7 % InCornerTL
                    A(x0y0,x0y0) = A(x0y0,x0y0)-2; 
                case 8 % InCornerTR
                    A(x0y0,x0y0) = A(x0y0,x0y0)-2; 
                otherwise
                    disp('error in constructing the matrix!')
                    return;
            end          
        end


    elseif (bcType==2) % clamped edge
        A(Index.Boundary,:) = I(Index.Boundary,:);   % w = 0 on boundary 
        A(Index.Corners,:) = I(Index.Corners,:);   % w = 0 on boundary 

         for i = 1:length(Index.InBoundaryCorner) 
            x0y0 = Index.InBoundaryCorner(i,1);

            switch Index.InBoundaryCorner(i,2) 
                case 1 % left boundary
                    A(x0y0,x0y0) = A(x0y0,x0y0)+1;               
                case 2 % right boundary
                    A(x0y0,x0y0) = A(x0y0,x0y0)+1;
                case 3 % bottom boundary
                    A(x0y0,x0y0) = A(x0y0,x0y0)+1;
                case 4 % top boundary
                    A(x0y0,x0y0) = A(x0y0,x0y0)+1; 
                case 5 % InCornerBL
                    A(x0y0,x0y0) = A(x0y0,x0y0)+2; 
                case 6 % InCornerBR
                    A(x0y0,x0y0) = A(x0y0,x0y0)+2; 
                case 7 % InCornerTL
                    A(x0y0,x0y0) = A(x0y0,x0y0)+2; 
                case 8 % InCornerTR
                    A(x0y0,x0y0) = A(x0y0,x0y0)+2; 
                otherwise
                    disp('error in constructing the matrix!')
                    return;
            end          
        end
  
     end
        
    if (mixedbcType==1 && bcType==2) % simply supported for center portions  (cannot include corners) with clamped bcs in the remainder
        
        A(Index.mixedBoundary,:) = I(Index.mixedBoundary,:);% w = 0 on boundary 
        
       for i = 1:length(Index.mixedInBoundaryIdx) 
            x0y0 = Index.mixedInBoundaryIdx(i,1);
            
            if(x0y0>0)

                switch Index.mixedInBoundaryIdx(i,2) 
                    case 1 % bottom boundary
                         A(x0y0,x0y0) = A(x0y0,x0y0)-2;               
                    case 2 % top boundary
                        A(x0y0,x0y0) = A(x0y0,x0y0)-2;    
                    case 3 % right boundary
                        A(x0y0,x0y0) = A(x0y0,x0y0)-2;    
                    case 4 % left boundary
                         A(x0y0,x0y0) = A(x0y0,x0y0)-2;    
                    otherwise
                        disp('error in constructing the matrix!')
                        return;
                end
            end
           
       end
        
        singular_stencil = ClampedSimplysupportedCoefficients(parameter.hx);

             for i = 1:8
                x0y0 = Index.mixedBoundaryEnd(i,1);
                x0ypp = x0y0-2;
                xmyp = x0y0-ny-1;
                x0yp = x0y0-1;
                xpyp = x0y0+ny-1;
                xmmy0 = x0y0-2*ny;             
                xmy0 = x0y0-ny;
                xpy0 = x0y0+ny;           
                xppy0 = x0y0+2*ny;
                xmym = x0y0-ny+1;
                x0ym = x0y0+1;
                xpym = x0y0+ny+1;
                x0ymm = x0y0+2;
                xpymm = x0y0+ny+2;
                xmymm = x0y0-ny+2;
                xppym = x0y0+2*ny+1;
                xmmym = x0y0-2*ny+1;
                xpypp = x0y0+ny-2;
                xmypp = x0y0-ny-2;
                xmmyp = x0y0-2*ny-1;
                xppyp = x0y0+2*ny-1;

                if x0y0>0 % address singularity where BCs changes by analytical solution                 
                    switch Index.mixedBoundaryEnd(i,2)
                    case 1 % mixedB1
                        A(x0yp,:)=0;
                        idx_vec = [x0yp,xpyp,xmyp,x0ypp,xppyp,xmmyp,xpypp,xmypp];
                        A(x0yp,idx_vec) = singular_stencil;
                    case 2 % mixedB2
                        A(x0yp,:)=0;
                        idx_vec = [x0yp,xmyp,xpyp,x0ypp,xmmyp,xppyp,xmypp,xpypp];
                        A(x0yp,idx_vec) = singular_stencil;
                    case 3 % mixedT1
                        A(x0ym,:)=0;
                        idx_vec = [x0ym,xpym,xmym,x0ymm,xppym,xmmym,xpymm,xmymm];
                        A(x0ym,idx_vec) = singular_stencil;
                    case 4 % mixedT2
                        A(x0ym,:)=0;
                        idx_vec = [x0ym,xmym,xpym,x0ymm,xmmym,xppym,xmymm,xpymm];
                        A(x0ym,idx_vec) = singular_stencil;
                    case 5 % mixedR1
                        A(xmy0,:)=0;
                        idx_vec = [xmy0,xmym,xmyp,xmmy0,xmymm,xmypp,xmmym,xmmyp];
                        A(xmy0,idx_vec) = singular_stencil;
                    case 6 % mixedR2
                        A(xmy0,:)=0;
                        idx_vec = [xmy0,xmyp,xmym,xmmy0,xmypp,xmymm,xmmyp,xmmym];
                        A(xmy0,idx_vec) = singular_stencil;
                    case 7 % mixedL1
                        A(xpy0,:)=0;
                        idx_vec = [xpy0,xpym,xpyp,xppy0,xpymm,xpypp,xppym,xppyp];
                        A(xpy0,idx_vec) = singular_stencil;
                    case 8 % mixedL2
                        A(xpy0,:)=0;
                        idx_vec = [xpy0,xpyp,xpym,xppy0,xpypp,xpymm,xppyp,xppym];
                        A(xpy0,idx_vec) = singular_stencil;
                    otherwise
                        disp('error in constructing the matrix!')
                        return;
                    end
                end
            end

    elseif (mixedbcType==2 && bcType==1) % simply supported for center portions  (cannot include corners) with clamped bcs in the remainder

        A(Index.mixedBoundary,:) = I(Index.mixedBoundary,:);% w = 0 on boundary 

       for i = 1:length(Index.mixedInBoundaryIdx) 
            x0y0 = Index.mixedInBoundaryIdx(i,1);
            
            if(x0y0>0)

                switch Index.mixedInBoundaryIdx(i,2) 
                    case 1 % bottom boundary
    %                     x0ymm = x0y0;
                        A(x0y0,x0y0) = A(x0y0,x0y0)+2;               
                    case 2 % top boundary
                        A(x0y0,x0y0) = A(x0y0,x0y0)+2;               
                    case 3 % right boundary
                        A(x0y0,x0y0) = A(x0y0,x0y0)+2;  
                    case 4 % left boundary
                        A(x0y0,x0y0) = A(x0y0,x0y0)+2; 
                    otherwise
                        disp('error in constructing the matrix!')
                        return;
                end
            end
           
       end  
       
        singular_stencil = ClampedSimplysupportedCoefficients(parameter.hx);

         for i = 1:8
            x0y0 = Index.mixedBoundaryEnd(i,1);
            x0ypp = x0y0-2;
            xmyp = x0y0-ny-1;
            x0yp = x0y0-1;
            xpyp = x0y0+ny-1;
            xmmy0 = x0y0-2*ny;             
            xmy0 = x0y0-ny;
            xpy0 = x0y0+ny;           
            xppy0 = x0y0+2*ny;
            xmym = x0y0-ny+1;
            x0ym = x0y0+1;
            xpym = x0y0+ny+1;
            x0ymm = x0y0+2;
            xpymm = x0y0+ny+2;
            xmymm = x0y0-ny+2;
            xppym = x0y0+2*ny+1;
            xmmym = x0y0-2*ny+1;
            xpypp = x0y0+ny-2;
            xmypp = x0y0-ny-2;
            xmmyp = x0y0-2*ny-1;
            xppyp = x0y0+2*ny-1;

            if x0y0>0 % address singularity where BCs changes by analytical solution                 
                switch Index.mixedBoundaryEnd(i,2)
                case 1 % mixedB1
                    A(x0yp,:)=0;
                    idx_vec = [x0yp,xmyp,xpyp,x0ypp,xmmyp,xppyp,xmypp,xpypp];
                    A(x0yp,idx_vec) = singular_stencil;
                case 2 % mixedB2
                    A(x0yp,:)=0;

                    idx_vec = [x0yp,xpyp,xmyp,x0ypp,xppyp,xmmyp,xpypp,xmypp];
                    A(x0yp,idx_vec) = singular_stencil;
                case 3 % mixedT1
                    A(x0ym,:)=0;
                    idx_vec = [x0ym,xmym,xpym,x0ymm,xmmym,xppym,xmymm,xpymm];
                    A(x0ym,idx_vec) = singular_stencil;
                case 4 % mixedT2
                    A(x0ym,:)=0;

                    idx_vec = [x0ym,xpym,xmym,x0ymm,xppym,xmmym,xpymm,xmymm];
                    A(x0ym,idx_vec) = singular_stencil;
                case 5 % mixedR1
                    A(xmy0,:)=0;
                    idx_vec = [xmy0,xmyp,xmym,xmmy0,xmypp,xmymm,xmmyp,xmmym];
                    A(xmy0,idx_vec) = singular_stencil;
                case 6 % mixedR2
                    A(xmy0,:)=0;

                    idx_vec = [xmy0,xmym,xmyp,xmmy0,xmymm,xmypp,xmmym,xmmyp];
                    A(xmy0,idx_vec) = singular_stencil;
                case 7 % mixedL1
                    A(xpy0,:)=0;
                    idx_vec = [xpy0,xpyp,xpym,xppy0,xpypp,xpymm,xppyp,xppym];
                    A(xpy0,idx_vec) = singular_stencil;
                case 8 % mixedL2
                    A(xpy0,:)=0;

                    idx_vec = [xpy0,xpym,xpyp,xppy0,xpymm,xpypp,xppym,xppyp];
                    A(xpy0,idx_vec) = singular_stencil;
                otherwise
                    disp('error in constructing the matrix!')
                    return;
                end
            end
         end
        
 
        
    end
    
end

