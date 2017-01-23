%addpath C:\Chisco\Proyectos\MARUCA\Soft\Ec_Dispersion\

function resultados=LecturaResultadosPropagaciones_Interpolacion_4parametros(path,caso,nodos,distancia)


load([path caso '_SF0131.mat']);

if sum(distancia)==0
    
    resultados(1)=Hsig(nodos(2,1),nodos(1,1));
    resultados(2)=Tm02(nodos(2,1),nodos(1,1));
    resultados(3)=RTpeak(nodos(2,1),nodos(1,1));
    resultados(4)=Dir(nodos(2,1),nodos(1,1));
    
else
    
    
    H=[Hsig(nodos(2,1),nodos(1,1)) Hsig(nodos(2,1),nodos(1,2)) Hsig(nodos(2,2),nodos(1,1)) Hsig(nodos(2,2),nodos(1,2))];
    
    for i=1:4
        if isnan(H(i))
            valido(i)=0;
        else
            valido(i)=1;
        end
    end
    
    T=[Tm02(nodos(2,1),nodos(1,1)) Tm02(nodos(2,1),nodos(1,2)) Tm02(nodos(2,2),nodos(1,1)) Tm02(nodos(2,2),nodos(1,2))];
    Tp=[RTpeak(nodos(2,1),nodos(1,1)) RTpeak(nodos(2,1),nodos(1,2)) RTpeak(nodos(2,2),nodos(1,1)) RTpeak(nodos(2,2),nodos(1,2))];
    D=[Dir(nodos(2,1),nodos(1,1)) Dir(nodos(2,1),nodos(1,2)) Dir(nodos(2,2),nodos(1,1)) Dir(nodos(2,2),nodos(1,2))];
         
    dum=find(valido==0);
    distancia(dum)=NaN;
    %Ponderacion por la inversa de distancia
    distT=nansum((1./distancia)');
    dist=(1./distancia)./distT;
    
    resultados(1)=nansum(H'.*dist'.*valido'); %Hs
    resultados(2)=nansum(T'.*dist'.*valido'); %T
    resultados(3)=nansum(Tp'.*dist'.*valido'); %Tp
    %Direccion media del oleaje y de la energia
    resultados(4)=180/pi*atan2(nansum(sin(D*pi/180)'.*dist'.*valido'),nansum(cos(D*pi/180)'.*dist'.*valido')); %Dir
    if resultados(4)<0 
        resultados(4)=resultados(4)+360; 
    end
    
end





