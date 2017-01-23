function [qerr bmu] = dist_normalizada(M, D, PosEscalares, PosDireccionales)  

     dMD(:,PosEscalares) = abs(D(:,PosEscalares)-M(:,PosEscalares));
     dMD(:,PosDireccionales) = min(abs(D(:,PosDireccionales)-M(:,PosDireccionales)),2*pi-abs(D(:,PosDireccionales)-M(:,PosDireccionales)))./pi; 
        
     Dist = abs(sum(dMD.^2,2));
     [qerr bmu] = min(Dist);
    
     