clearvars
close all

%clear all; %esborrem tot
x1=2.0; 
x2=5.0;
f = -3.0;
u1 = 3.0;
u2 = 8.0;

% x1=0.0;
% x2=2.0;
% u1 = 1.7;
% dudx2 = 2.0;
% f = 0.0;

L = x2-x1;
a1=3; a0=2; numElem=200; h=L/numElem;

%clear K; Why?
nodes=(x1:h:x2)'; %colum vector
%     for i=1:numElem
%            elem(i,:)=[i, i+1]; %creo elements i numeració global
%            %exemple e1=[1,2]
%     end   
for e = 1:numElem
    elem(e,:) = [e,e+1];
end

numNod = size(nodes,1);
numElem = size(elem,1);

K=zeros(numNod);
Q = zeros(numNod,1);
F = zeros(numNod,1);
u = zeros(numNod,1);

Fe = f*h*[1;1]/2.0;

for i=1:numElem
    nod1=elem(i,1);
    nod2=elem(i,2);
    K1=(a1/(3*h^2))*((nodes(nod2))^3-(nodes(nod1))^3)*[1,-1;-1,1];
    K0=a0*h/6*[2,1;1,2];
    Ke = K1+K0;
    rows=[elem(i,1);elem(i,2)];
    colums=rows;
    K(rows,colums)=K(rows,colums)+Ke;
    F(rows) = F(rows)+Fe;
end
%solucio= K(50,50)   %Sol. correcte 2.9921e+0
%hint=K(20,21)       %Sol. correcte -1.0511e+03
    
fixedNodes = [1,numNod];
freeNodes = setdiff(1:numNod, fixedNodes);
   
%Natural B.C.
u(fixedNodes) = [u1;u2];
    
%Essential B.C.
Q(freeNodes) = 0.0; %**In this case** not necessary, but...   
    
%Reduced system
Fm = Q(freeNodes) - K(freeNodes,fixedNodes)*u(fixedNodes);
Fm = Fm + F(freeNodes);
    
Km = K(freeNodes,freeNodes);
    
%solve the reduced system
um = Km\Fm;  
u(freeNodes)= um;
    
%Post process:
Q = K*u - F;
     
%Compute the mean value 
meanValueOfU = sum(u)/numNod;
idxs = find(abs(u - meanValueOfU) < 0.1);
numNodsIdxs =length(idxs);
   
fprintf('K50,50) = %.5e\n',K(50,50));
fprintf('Hint. K(20,21) = %.5e\n',K(20,21));
fprintf('<u> = %.5e\n',meanValueOfU);
fprintf('  N = %d\n',numNodsIdxs);