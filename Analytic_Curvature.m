clear all
clc
format long

linedata = importdata('data1.txt');
x = linedata(:,2);
y = linedata(:,3);

%% contadores
    nlines = 0; % contador de streamlines
    
for i=1:length(x)-1
    if x(i+1) < x(i)
        nlines = nlines +1 ;
    else
        nlines = nlines + 0;
    end    
end

 
contador = 1;
j = 1;
%lineas = {nlines};

for i=2:length(x)
    
    if x(i) == 0;
        j = j+1;
        contador = 1;
    else
        j = j+0;
        contador = contador +1;
    end
    
    lineas{j}(contador,1) = x(i);
    lineas{j}(contador,2) = y(i);
    
    if j == nlines
        break;
    else
        
    end
      
end

for i=1:nlines 
   sizearr = length(lineas{i});
   plot(lineas{i}(:,1),lineas{i}(:,2))
   hold on 
end

%% CALCULO DE CURVATURA %%

%primera derivada 
Dx = 0.01;

for i=1:nlines % recorro todas las lineas de corriente en la direccion vertical.
    
    for j=2:(length(lineas{i})-1)
        if abs(lineas{i}(j+1,1) - lineas{i}(j,1)) < Dx
            
            Ub{i}(j-1,1) = lineas{i}(j,1);
            Ub{i}(j-1,2) = lineas{i}(j,2);
            
            Ud = lineas{i}(j+1,1) - lineas{i}(j,1);
            Ub{i}(j-1,3) = Ud; 
            Ub{i}(j-1,4) = lineas{i}(j+1,2);

        else
            
            Ud = Dx;
            
            Ub{i}(j-1,1) = lineas{i}(j,1);
            Ub{i}(j-1,2) = lineas{i}(j,2);  
            
            Ub{i}(j-1,3) = Ud;
            m = (lineas{i}(j+1,2) - lineas{i}(j,2))/(lineas{i}(j+1,1) - lineas{i}(j,1));
            Ub{i}(j-1,4) = lineas{i}(j,2) + m*Dx;
            
            
        end
        
        if abs(lineas{i}(j-1,1) - lineas{i}(j,1)) < Dx
            
            Lb{i}(j-1,1) = lineas{i}(j,1);
            Lb{i}(j-1,2) = lineas{i}(j,2);
            
            Ld = lineas{i}(j-1,1) - lineas{i}(j,1);
            Lb{i}(j-1,4) = lineas{i}(j-1,1);
            Lb{i}(j-1,3) = Ld;
        else
            
            Lb{i}(j-1,1) = lineas{i}(j,1);
            Lb{i}(j-1,2) = lineas{i}(j,2);  
            
            Ld = Dx;
            Lb{i}(j-1,3) = Ud;
            m = (lineas{i}(j-1,2) - lineas{i}(j,2))/(lineas{i}(j-1,1) - lineas{i}(j,1));
            Lb{i}(j-1,4) = lineas{i}(j,2) + m*Dx;            
            
            
        end
        
    end   
end
j = i
%% S E G U N D A    D E R I V A D A %%

for i=1:nlines-1 % recorro todas las lineas de corriente en la direccion vertical.
    
    for j=2:(length(Ub{i})-1) % recorro todos los elementos de la linea de corriente
        if abs(Ub{i}(j+1,1) - Ub{i}(j,1)) < Dx
           
            Ubtwo{i}(j-1,1) = lineas{i}(j+1,1);
            Ubtwo{i}(j-1,2) = lineas{i}(j+1,2);
            
            Udtwo = Ub{i}(j+1,1) - Ub{i}(j,1);
            Ubtwo{i}(j-1,3) = Ud; 
            Ubtwo{i}(j-1,4) = Ub{i}(j+1,2);
            
        else
           
            Ubtwo{i}(j-1,1) = lineas{i}(j+1,1);
            Ubtwo{i}(j-1,2) = lineas{i}(j+1,2);
            Udtwo = Dx;
            Ubtwo{i}(j-1,3) = Ud;
            m = (Ub{i}(j+1,2) - Ub{i}(j,2))/(Ub{i}(j+1,1) - Ub{i}(j,1));
            Ubtwo{i}(j-1,4) = Ub{i}(j,2) + m*Dx;
            
            
        end
        
        if abs(Lb{i}(j-1,1) - Lb{i}(j,1)) < Dx
            
            Lbtwo{i}(j-1,1) = lineas{i}(j+1,1);
            Lbtwo{i}(j-1,2) = lineas{i}(j+1,2);
            
            Ldtwo = Lb{i}(j-1,1) - Lb{i}(j,1);
            Lbtwo{i}(j-1,4) = Lb{i}(j-1,1);
            Lbtwo{i}(j-1,3) = Ld;
        else
          
            Lbtwo{i}(j-1,1) = lineas{i}(j,1);
            Lbtwo{i}(j-1,2) = lineas{i}(j,2);            
            Ldtwo = Dx;
            Lbtwo{i}(j-1,3) = Ud;
            m = (Lb{i}(j-1,2) - Lb{i}(j,2))/(Lb{i}(j-1,1) - Lb{i}(j,1));
            Lbtwo{i}(j-1,4) = Lb{i}(j,2) + m*Dx;            
                     
        end
        
    end   
end

%% ACOPLAMIENTO DERIVADAS Y CALCULO DE RADIO DE CURVATURA %%

for i=1:length(Ub)
    for j=1:length(Ub{i})
        ddxsl{i}(j,1) = (Ub{i}(j,4) - Lb{i}(j,4))/(Ub{i}(j,3) + Lb{i}(j,3));
    end
end

for i=1:length(Ubtwo)
    for j=1:length(Ubtwo{i})
        ddx2sl{i}(j,1) = (Ubtwo{i}(j,4) - Lbtwo{i}(j,4))/(Ubtwo{i}(j,3) + Lbtwo{i}(j,3));
    end
end


for i=1:nlines -1
    for j=1:length(Ubtwo{i})  
        
        k{i}(j,1) = (abs(ddx2sl{i}(j,1)))/(((1 + ((ddxsl{i}(j+1,1))^(2))))^(3/2));
        
    end   
end

for i=1:length(Lbtwo)
    for j=1:(length(Lbtwo{i})-1)
    
      if Lbtwo{i}(j+1,1)>5 
          
         dx = (Lbtwo{i}(j+1,1) - Lbtwo{i}(j,1)) ;
         Dx = 5 - k{i}(j,1);
         m = (k{i}(j+1,1) - k{i}(j,1))/dx;
         
             curvature(i,1) = Lbtwo{i}(j,2);
             curvature(i,2) = k{i}(j,1) + m*Dx;
         
         break
      else
          continue
      end
      
    end
end
curvature(1,1) = 0.1;
curvature(1,2) = 0;

%%curl = log10(1./curvature(:,2));
%%curled = log10(curvature(:,2));
verticalline = [0 0.1 ; 0 2];

figure(2)
scatter(abs(curvature(:,2)),curvature(:,1),'filled');

figure(3)
plot(verticalline(:,1),verticalline(:,2),'b','LineWidth',2.5);
hold on
plot((curvature(:,2)),curvature(:,1),'r','LineWidth',2.5);
xlabel('Curvatura (k)')
ylabel('Direccion vertical [m]')



%%figure(4)
%%plot(abs(curl),curvature(:,1),'r');

%%figure(5)
%%plot(abs(curled),curvature(:,1),'r');

    
    %%plot(x,y);
%hold on 
% scatter(x,y,'r')