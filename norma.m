function [Dist] = norma (M,D)  
    
    %Posiciones de las variables escalares de los puntos en indefinidas
    PosEscalares=[1 2 4 5];
    %Posiciones de las variables direccionales de los puntos en indefinidas
    PosDireccionales=[3 6];

    dx(PosEscalares,:)=M(PosEscalares,:) - D(PosEscalares,:); 
    dx(PosDireccionales,:)=min(abs(D(PosDireccionales,:)-M(PosDireccionales,:)),2*pi-abs(D(PosDireccionales,:)-M(PosDireccionales,:)))./pi;
                   
    Dist=sqrt(sum(abs(dx.^2),1));
