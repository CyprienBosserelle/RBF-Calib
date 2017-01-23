
function yy = evaluarmin( x, y, limmin, limmax)

 
[yy,fval] = fminbnd(@(ep)CostEps(ep,x,y),limmin,limmax,optimset('TolX',0.0001));
    
    function yy = CostEps(ep,x,y)
    
        [m,n]=size(x);
        coeff = rbfcreate_modificado_matrizA(x,y,'RBFFunction','gaussian','RBFConstant',ep);
        A1=coeff.matrizA;
        A=A1(1:n,1:n);
        invA=pinv(A);
        [m1,n1]=size(coeff.rbfcoeff);
        KK = y - coeff.rbfcoeff(m1-m);
        for i=1:m
            KK = KK - coeff.rbfcoeff(m1-m+i).*x(i,:);
        end
        ceps = (invA*KK')./diag(invA);
        yy=norm(ceps);
        
    end

end
        
    