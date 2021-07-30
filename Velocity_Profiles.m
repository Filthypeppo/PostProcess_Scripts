clear all 
clc
numfiles = 11;
x = 4.5:0.1:5.5;
y = zeros(100,22,2); % primera dimension en z es para la posicion 
                     % segunda posicion en z es para la velocidad
Uref = 1; %velocidad aproximada
a = 0.1; % lado de la figura [m]
refx = 5; % posicion del centro del objeto en el eje "X";
for z=1:length(x)
    pos(z) = (x(z)-refx)/a;
    
end

for k=1:numfiles
    myfilename = sprintf('lineX%d_U.txt', k);
    mydata{k} = importdata(myfilename);
end

for j=1:(numfiles)
    for i=1:length(mydata{j})
    y(i,j,1) = mydata{j}(i,1);
    y(i,j,2) = mydata{j}(i,2);%/mydata{2}(i,2);
    %y(i,j,2) = sqrt(((mydata{j}(i,2))^2)+((mydata{j}(i,3))^2)+((mydata{j}(i,4))^2))/Uref;
    end
end

for i=1:(numfiles)
    plot(pos(i)+y(:,i,2),y(:,i,1));
    hold on   
end
ylabel({'$y$'},'Interpreter','latex','FontSize',20)
title({'${\bar{Ux}/U_{\infty}}$'},'Interpreter','latex','FontSize',20)
xlabel('$x$','Interpreter','latex','FontSize',35)
xlim([-10 45]);

figure(2)
for i=1:(numfiles)

    subplot(1,numfiles,i),plot(y(:,i,2),y(:,i,1));
    %%stdlim = std(y(:,i,1));
    %%minlim = min(y(:,i,1))-stdlim;
    %%maxlim = max(y(:,i,1)) + stdlim;
    %%xlim([minlim maxlim]);
    ylabel({'${y}$'},'Interpreter','latex','FontSize',10)
    xlabel({'${\bar{Ux}/U_{\infty}}$'},'Interpreter','latex','FontSize',10)
    title(['$x/a = $' num2str(pos(i))],'Interpreter','latex','FontSize',15)
    hold on 

end

