function A=assignBoundaryConditionsCoefficient(A,Index,parameter)
%assign boundary conditions for the coefficient matrix
%Input:
% A: differentiation matrix for biharmonic operator
% Index: 
% diffMtx: differentiation matrices
% bcType: boundary condition types: simply supported=1, clamped=2, free=3

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
        diagonalIdx = (Index.InBoundary-1)*(nx*ny+1)+1;
        A(diagonalIdx) = 19;
        diagonalIdx = (Index.InCorners-1)*(nx*ny+1)+1;
        A(diagonalIdx) = 18;
        
        diagonalIdx = (Index.InB1-1)*(nx*ny+1)+1+nx*ny;
        A(diagonalIdx) = -8;
        
        diagonalIdx = (Index.InR1-1)*(nx*ny+1)+1+nx*ny*ny;
        A(diagonalIdx) = -8;
        
        A(Index.InCornerBL,Index.InCornerBL+1) = -8;
        A(Index.InCornerBR,Index.InCornerBR+1) = -8;
        A(Index.InCornerBR,Index.InCornerBR+ny) = -8;
        A(Index.InCornerTR,Index.InCornerTR+ny) = -8;


    elseif (bcType==2) % clamped edge
       % disp('clamped bc'); 
        A(Index.Boundary,:) = I(Index.Boundary,:);   % w = 0 on boundary 
        A(Index.Corners,:) = I(Index.Corners,:);   % w = 0 on boundary 
        diagonalIdx = (Index.InBoundary-1)*(nx*ny+1)+1;
        A(diagonalIdx) = 21;

        diagonalIdx = (Index.InCorners-1)*(nx*ny+1)+1;
        A(diagonalIdx) = 22;
        
        diagonalIdx = (Index.InB1-1)*(nx*ny+1)+1+nx*ny;
        A(diagonalIdx) = -8;
        
        diagonalIdx = (Index.InR1-1)*(nx*ny+1)+1+nx*ny*ny;
        A(diagonalIdx) = -8;
        
        A(Index.InCornerBL,Index.InCornerBL+1) = -8;
        A(Index.InCornerBR,Index.InCornerBR+1) = -8;
        A(Index.InCornerBR,Index.InCornerBR+ny) = -8;
        A(Index.InCornerTR,Index.InCornerTR+ny) = -8;
        

    elseif (bcType==3) % free edge
        %disp('free bc');
        A(Index.Boundary,:) = 0;   
        A(Index.Corners,:) = 0;
        A(Index.InBoundaryCorner(:,1),:) = 0;

        tmp0 = 16-8*nu-6*nu^2;
        tmp1 = -8+4*nu+4*nu^2;
        tmp2 = 1-nu^2;
        tmp3 = 4-2*nu;
        tmp4 = -12+4*nu;

        tmp5 = 12-4*nu^2-8*nu;
        tmp6 = -12+8*nu+4*nu^2;
        tmp7 = 8-8*nu;
        tmp8 = 2-2*nu^2;
        
        tmp9 = 15-8*nu-5*nu^2;
        tmp10 = 2*nu^2+4*nu-6;
        tmp11 = 4*nu^2+4*nu-8;
        tmp12 = 4-2*nu;
        tmp13 = 4*nu-12;


        for i = 1:length(Index.BoundaryCorners)
            x0y0 = Index.BoundaryCorners(i,1);
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
            
             
            switch Index.BoundaryCorners(i,2)
                case 1 % left boundary
                    A(x0y0,[x0y0,x0yp,x0ym,x0ypp,x0ymm,xpyp,xpym,xpy0,xppy0]) = [tmp0,tmp1,tmp1,tmp2,tmp2,tmp3,tmp3,tmp4,2];
                case 2 % right boundary
                    A(x0y0,[x0y0,x0yp,x0ym,x0ypp,x0ymm,xmyp,xmym,xmy0,xmmy0]) = [tmp0,tmp1,tmp1,tmp2,tmp2,tmp3,tmp3,tmp4,2];
                case 3 % bottom boundary
                    A(x0y0,[x0y0,xpy0,xmy0,xppy0,xmmy0,xpyp,xmyp,x0yp,x0ypp]) = [tmp0,tmp1,tmp1,tmp2,tmp2,tmp3,tmp3,tmp4,2];
                case 4 % top boundary
                    A(x0y0,[x0y0,xpy0,xmy0,xppy0,xmmy0,xmym,xpym,x0ym,x0ymm]) = [tmp0,tmp1,tmp1,tmp2,tmp2,tmp3,tmp3,tmp4,2];
                case 5 % CornerBL
                    A(x0y0,[x0y0,x0yp,xpy0,xpyp,x0ypp,xppy0]) = [tmp5,tmp6,tmp6,tmp7,tmp8,tmp8];
                case 6 % CornerBR
                    A(x0y0,[x0y0,x0yp,xmy0,xmyp,x0ypp,xmmy0]) = [tmp5,tmp6,tmp6,tmp7,tmp8,tmp8];
                case 7 % CornerTL
                    A(x0y0,[x0y0,x0ym,xpy0,xpym,x0ymm,xppy0]) = [tmp5,tmp6,tmp6,tmp7,tmp8,tmp8];
                case 8 % CornerTR
                    A(x0y0,[x0y0,x0ym,xmy0,xmym,x0ymm,xmmy0]) = [tmp5,tmp6,tmp6,tmp7,tmp8,tmp8];
                case 9 % CornerBLL
                    A(x0y0,[x0y0,x0ym,x0yp,xpyp,xpym,xpy0,x0ypp,xppy0]) = [tmp9,tmp10,tmp11,tmp12,tmp12,tmp13,tmp2,2];
                case 10 %CornerBLB
                    A(x0y0,[x0y0,xmy0,xpy0,xmyp,xpyp,x0yp,xppy0,x0ypp]) = [tmp9,tmp10,tmp11,tmp12,tmp12,tmp13,tmp2,2];
                case 11 %CornerBRR
                    A(x0y0,[x0y0,x0ym,x0yp,xmym,xmyp,xmy0,x0ypp,xmmy0]) = [tmp9,tmp10,tmp11,tmp12,tmp12,tmp13,tmp2,2];
                case 12 %CornerBRB
                    A(x0y0,[x0y0,xpy0,xmy0,xmyp,xpyp,x0yp,xmmy0,x0ypp]) = [tmp9,tmp10,tmp11,tmp12,tmp12,tmp13,tmp2,2];
                case 13 %CornerTLL
                    A(x0y0,[x0y0,x0yp,x0ym,xpym,xpyp,xpy0,x0ymm,xppy0]) = [tmp9,tmp10,tmp11,tmp12,tmp12,tmp13,tmp2,2];
                case 14 %CornerTLT
                    A(x0y0,[x0y0,xmy0,xpy0,xpym,xmym,x0ym,xppy0,x0ymm]) = [tmp9,tmp10,tmp11,tmp12,tmp12,tmp13,tmp2,2];
                case 15 %CornerTRR
                    A(x0y0,[x0y0,x0yp,x0ym,xmym,xmyp,xmy0,x0ymm,xmmy0]) = [tmp9,tmp10,tmp11,tmp12,tmp12,tmp13,tmp2,2];
                case 16 %CornerTRT
                    A(x0y0,[x0y0,xpy0,xmy0,xmym,xpym,x0ym,xmmy0,x0ymm]) = [tmp9,tmp10,tmp11,tmp12,tmp12,tmp13,tmp2,2];
                    
                otherwise
                    disp('error in constructing the matrix!')
                    return;
            end
           
        end
        

        for i = 1:length(Index.InBoundaryCorner) 
            x0y0 = Index.InBoundaryCorner(i,1);
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

            switch Index.InBoundaryCorner(i,2)
                case 1 % left boundary
                    A(x0y0,[x0y0,xpy0,x0yp,x0ym,xmy0,xpyp,xpym,xmyp,xmym,xppy0,x0ypp,x0ymm]) = [19,-8,-8,-8,-6+2*nu,2,2,2-nu,2-nu,1,1,1];
                case 2 % right boundary
                    A(x0y0,[x0y0,xmy0,x0yp,x0ym,xpy0,xmyp,xmym,xpyp,xpym,xmmy0,x0ypp,x0ymm]) = [19,-8,-8,-8,-6+2*nu,2,2,2-nu,2-nu,1,1,1]; 
                case 3 % bottom boundary
                    A(x0y0,[x0y0,x0yp,xpy0,xmy0,x0ym,xmyp,xpyp,xmym,xpym,xmmy0,xppy0,x0ypp]) = [19,-8,-8,-8,-6+2*nu,2,2,2-nu,2-nu,1,1,1]; 
                case 4 % top boundary
                    A(x0y0,[x0y0,x0ym,xpy0,xmy0,x0yp,xmym,xpym,xmyp,xpyp,xmmy0,xppy0,x0ymm]) = [19,-8,-8,-8,-6+2*nu,2,2,2-nu,2-nu,1,1,1]; 
                case 5 % InCornerBL
                    A(x0y0,[x0y0,xmy0,x0ym,x0yp,xpy0,xmym,xmyp,xpym,xpyp,x0ypp,xppy0]) = [18,-6+2*nu,-6+2*nu,-8,-8,2-2*nu,2-nu,2-nu,2,1,1];
                case 6 % InCornerBR
                    A(x0y0,[x0y0,xpy0,x0ym,x0yp,xmy0,xpym,xmym,xpyp,xmyp,x0ypp,xmmy0]) = [18,-6+2*nu,-6+2*nu,-8,-8,2-2*nu,2-nu,2-nu,2,1,1];
                case 7 % InCornerTL
                    A(x0y0,[x0y0,xmy0,x0yp,x0ym,xpy0,xmyp,xmym,xpyp,xpym,x0ymm,xppy0]) = [18,-6+2*nu,-6+2*nu,-8,-8,2-2*nu,2-nu,2-nu,2,1,1];
                case 8 % InCornerTR
                    A(x0y0,[x0y0,xpy0,x0yp,x0ym,xmy0,xpyp,xmyp,xpym,xmym,x0ymm,xmmy0]) = [18,-6+2*nu,-6+2*nu,-8,-8,2-2*nu,2-nu,2-nu,2,1,1];
                otherwise
                    disp('error in constructing the matrix!')
                    return;
            end
           
        end
        
     end

    if (mixedbcType==1) % simply supported  (cannot include corners)     

        %disp('mixed simply supported bc'); 
        
        A(Index.mixedBoundary,:) = I(Index.mixedBoundary,:);% w = 0 on boundary 
        A(Index.mixedInBoundary,:) = 0; 
        
       for i = 1:length(Index.mixedInBoundaryIdx) 
            x0y0 = Index.mixedInBoundaryIdx(i,1);
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
            
            if(x0y0>0)

                switch Index.mixedInBoundaryIdx(i,2) 
                    case 1 % bottom boundary
    %                     x0ymm = x0y0;
                        A(x0y0,[x0y0,xpy0,xmy0,x0yp,x0ym,xpyp,xmyp,xmym,xpym,xppy0,xmmy0,x0ypp]) = [19,-8,-8,-8,-8,2,2,2,2,1,1,1];
                    case 2 % top boundary
    %                     x0ypp = x0y0;
                        A(x0y0,[x0y0,xpy0,xmy0,x0yp,x0ym,xpyp,xmyp,xmym,xpym,xppy0,xmmy0,x0ymm]) = [19,-8,-8,-8,-8,2,2,2,2,1,1,1];
                    case 3 % right boundary
    %                     xppy0 = x0y0;
                         A(x0y0,[x0y0,xpy0,xmy0,x0yp,x0ym,xpyp,xmyp,xmym,xpym,xmmy0,x0ypp,x0ymm]) = [19,-8,-8,-8,-8,2,2,2,2,1,1,1];
                    case 4 % left boundary
    %                     xmmy0 = x0y0;
                        A(x0y0,[x0y0,xpy0,xmy0,x0yp,x0ym,xpyp,xmyp,xmym,xpym,xppy0,x0ypp,x0ymm]) = [19,-8,-8,-8,-8,2,2,2,2,1,1,1];
                    otherwise
                        disp('error in constructing the matrix!')
                        return;
                end
            end
           
        end

        
        switch bcType % address singularity where BCs changes by analytical solution              
            case 1 % simply supported bc,do nothing

            case 2 % clamped bc
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
                
            case 3 % free bc
                
 
                    
            otherwise
            return;

         end


    elseif (mixedbcType==2) % clamped edge (cannot include corners)     
        
      %disp('mixed simply supported bc'); 
        
        A(Index.mixedBoundary,:) = I(Index.mixedBoundary,:);% w = 0 on boundary 
        A(Index.mixedInBoundary,:) = 0; 
        
       for i = 1:length(Index.mixedInBoundaryIdx) 
            x0y0 = Index.mixedInBoundaryIdx(i,1);
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
            
            if(x0y0>0)

                switch Index.mixedInBoundaryIdx(i,2) 
                    case 1 % bottom boundary
    %                     x0ymm = x0y0;
                        A(x0y0,[x0y0,xpy0,xmy0,x0yp,x0ym,xpyp,xmyp,xmym,xpym,xppy0,xmmy0,x0ypp]) = [21,-8,-8,-8,-8,2,2,2,2,1,1,1];
                    case 2 % top boundary
    %                     x0ypp = x0y0;
                        A(x0y0,[x0y0,xpy0,xmy0,x0yp,x0ym,xpyp,xmyp,xmym,xpym,xppy0,xmmy0,x0ymm]) = [21,-8,-8,-8,-8,2,2,2,2,1,1,1];
                    case 3 % right boundary
    %                     xppy0 = x0y0;
                         A(x0y0,[x0y0,xpy0,xmy0,x0yp,x0ym,xpyp,xmyp,xmym,xpym,xmmy0,x0ypp,x0ymm]) = [21,-8,-8,-8,-8,2,2,2,2,1,1,1];
 
                    case 4 % left boundary
    %                     xmmy0 = x0y0;
                        A(x0y0,[x0y0,xpy0,xmy0,x0yp,x0ym,xpyp,xmyp,xmym,xpym,xppy0,x0ypp,x0ymm]) = [21,-8,-8,-8,-8,2,2,2,2,1,1,1];
                    otherwise
                        disp('error in constructing the matrix!')
                        return;
                end
            end
           
        end        
        switch bcType % address singularity where BCs changes by analytical solution              
            case 1 % simply supported bc
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
                
            case 2 % clamped bc, do nothing

            case 3 % free bc
                
                singular_stencil = ClampedFreeCoefficient(parameter.hx,parameter.nu);
   
                 for i = 1:8
                    x0y0 = Index.mixedBoundaryEnd(i,1);
                    xmmy0 = x0y0-2*ny;             
                    xmy0 = x0y0-ny;
                    xpy0 = x0y0+ny;           
                    xppy0 = x0y0+2*ny;
                    x0ypp = x0y0-2;
                    xmyp = x0y0-ny-1;
                    x0yp = x0y0-1;
                    xpyp = x0y0+ny-1;
                    xmmyp = x0y0-2*ny-1;
                    xppyp = x0y0+2*ny-1;
                    xmypp = x0y0-ny-2;
                    xpypp = x0y0+ny-2;
                    xmmypp = x0y0-2*ny-2;
                    xppypp = x0y0+2*ny-2;
                    xmym = x0y0-ny+1;
                    x0ym = x0y0+1;
                    xpym = x0y0+ny+1;
                    x0ymm = x0y0+2;
                    xmmym = x0y0-2*ny+1;
                    xppym = x0y0+2*ny+1;
                    xmymm = x0y0-ny+2;
                    xpymm = x0y0+ny+2;
                    xmmymm = x0y0-2*ny+2;
                    xppymm = x0y0+2*ny+2;

                    if x0y0>0 % address singularity where BCs changes by analytical solution   
                        switch Index.mixedBoundaryEnd(i,2)
                            case 1 % mixedB1
                                A(x0yp,:)=0;
                                idx_vec = [x0yp,xmyp,xpyp,x0ypp,xmmyp,xppyp];                           
                                A(x0yp,idx_vec) = singular_stencil;

                            case 2 % mixedB2
                                A(x0yp,:)=0;
                                idx_vec = [x0yp,xpyp,xmyp,x0ypp,xppyp,xmmyp];                
                                A(x0yp,idx_vec) = singular_stencil;
                               
                            case 3 % mixedT1
                                A(x0ym,:)=0;
                                idx_vec = [x0ym,xmym,xpym,x0ymm,xmmym,xppym];
                                A(x0ym,idx_vec) = singular_stencil;
                                     
                            case 4 % mixedT2
                                A(x0ym,:)=0;
                                idx_vec = [x0ym,xpym,xmym,x0ymm,xppym,xmmym];
                                A(x0ym,idx_vec) = singular_stencil;
                                                                
                            case 5 % mixedR1
                                A(xmy0,:)=0;
                                idx_vec = [xmy0,xmyp,xmym,xmmy0,xmypp,xmymm];                
                                A(xmy0,idx_vec) = singular_stencil;
                                                                
                            case 6 % mixedR2
                                A(xmy0,:)=0;
                                idx_vec = [xmy0,xmym,xmyp,xmmy0,xmymm,xmypp];                
                                A(xmy0,idx_vec) = singular_stencil;
                                
                             case 7 % mixedL1
                                A(xpy0,:)=0;
                                idx_vec = [xpy0,xpyp,xpym,xppy0,xpypp,xpymm];                
                                A(xpy0,idx_vec) = singular_stencil;
                                                                
                            case 8 % mixedL2
                                A(xpy0,:)=0;
                                idx_vec = [xpy0,xpym,xpyp,xppy0,xpymm,xpypp];                
                                A(xpy0,idx_vec) = singular_stencil;
                                
                            otherwise
                                disp('error in constructing the matrix!')
                        end
                    end
                    
                 end
            otherwise
                disp('BC Type not supported!')
            return;

         end
        
    end
    
end

