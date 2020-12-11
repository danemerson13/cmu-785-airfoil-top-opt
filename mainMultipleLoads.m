function [fval,exitflag,loop,maxStress,x,volfrac,lambda,output] = mainMultipleLoads(volfrac,penal,rmin)

% clc; 
close all; clearvars -except volfrac penal rmin;
% volfrac = 0.6;
% penal = 4.5;
% rmin = 0.006;

load('Airfoil_Uniform_1100_MeshData_v3.mat','ElementList','ElementList_Active_Bool',...
    'NodeList','NodeList_Active','NodeList_Active','NodeList_FixedDOF','NodeList_FixedDOF_Bool',...
    'NodeList_Spar1','NodeList_Spar1_Bool','NodeList_Spar2','NodeList_Spar2_Bool','NodeList_Surface');
load('LoadingZones.mat','zoneLeadBot','zoneLeadTop','zoneTrailBot','zoneTrailTop','zoneMountFront','zoneMountRear');

% INITIALIZE
% Find which elements are active and which nodes are fixed
ElementList_Active_Bool = [ElementList(:,1), any(ismember(ElementList(:,2:end), NodeList_Active), 2)];
NodeList_FixedDOF_Bool = [NodeList(:,1), ismember(NodeList(:,1), NodeList_FixedDOF)];

% Material properties
E = 68.9e9;
nu = 0.33;
oy = 276e6;

% Number of Elements and Nodes
nE = size(ElementList,1);
nN = size(NodeList, 1);

% Initialize x0 to all have density of volfrac
x = volfrac*ones(nE, 1);

% Adds fixed active areas to the airfoil, e.g. boundaries of airfoil
% and spars, adds to initial x matrix, added 25 Nov 2020 AJB
for a = 1:nE
    if ElementList_Active_Bool(a,2) == 1
        x(a) = 1;
    end
end

%Matrix for storing and plotting von Mises stresses
stressMat = zeros(nE, 1);

%Pre-compile element stiffness matrices
k_elements = zeros(8,8,nE); 
for a = 1:nE
   nodes_i = ElementList(a,[2 3 4 5 2]);
   nodePos_i = NodeList(nodes_i, 2:3)';
   k_elements(:,:, a) = calcStiffness(nodePos_i, E, nu);
end

%Get areas and centroids of elements
[area_e, centroid_e] = calcArea(ElementList,NodeList);

% Setup figures for plotting volume fraction x and vonMises stress
figure('Units', 'normalized', 'Position', [0.5, 0, 0.5, 1]);
ax(1) = subplot(2,1,1);
plot1 =  meshPlot(ElementList,NodeList, ElementList_Active_Bool, NodeList_FixedDOF_Bool, -x);
colorbar(ax(1))
title('Volume Fraction')
axis equal; axis tight; axis off;

ax(2) = subplot(2,1,2);
plot2 = meshPlot(ElementList,NodeList, ElementList_Active_Bool, NodeList_FixedDOF_Bool, stressMat);
axis equal; axis tight; axis off;
colorbar(ax(2));
title('von Mises Stress')
caxis(ax(2),[0, oy]);

% Initialize while loop counting and termination variables
loop = 0; 
change = 1.;
% START ITERATION
while change > 0.01  % Loop until the design does not change significantly. This is typical termination criteria in top-opt
    loop = loop + 1;
    xold = x;
    
    %FE-ANALYSIS
    [U, F]=FE(nE, nN, ElementList, x,penal, k_elements, NodeList_FixedDOF, NodeList_Spar1, NodeList_Spar2, NodeList_Surface);

    %Plot boundary conditions once
    if loop == 1
        figure('Units', 'normalized', 'Position', [0,0,0.5,1])
        title("Active Elements and Boundary Conditions")
        subplot(3,1,1); title('Steady-State Load Case');
        subplot(3,1,2); title('Acceleration Load Case');
        subplot(3,1,3); title('Braking Load Case');
        for n = 1:3
            U_plotting = reshape(U(:,n), 2, [])';
            NodeListU = NodeList;
            NodeListU(:,2:3) = NodeListU(:,2:3) + 100*U_plotting;
            
            subplot(3,1,n); hold on;
            meshPlot(ElementList,NodeList, ElementList_Active_Bool, NodeList_FixedDOF_Bool);
            hold off

            NodePosQuiver1 = NodeList(NodeList_Spar1, 2:3);
            FQuiver1 = [F(NodeList_Spar1*2-1,n), F(NodeList_Spar1*2,n)];

            NodePosQuiver2 = NodeList(NodeList_Spar2, 2:3);
            FQuiver2 = [F(NodeList_Spar2*2-1,n), F(NodeList_Spar2*2,n)];

            NodePosQuiver3 = NodeList(NodeList_Surface, 2:3);
            FQuiver3 = [F(NodeList_Surface*2-1,n), F(NodeList_Surface*2,n)];

            NodePosQuiverTotal = [NodePosQuiver1; NodePosQuiver2; NodePosQuiver3];
            FQuiverTotal = [FQuiver1; FQuiver2; FQuiver3];

            subplot(3,1,n); hold on;
            if n == 2
                quiver(NodePosQuiverTotal(:,1), NodePosQuiverTotal(:,2), FQuiverTotal(:,1), FQuiverTotal(:,2),'b', 'LineWidth',2, 'AutoScaleFactor', 0.2)
            else
                quiver(NodePosQuiverTotal(:,1), NodePosQuiverTotal(:,2), FQuiverTotal(:,1), FQuiverTotal(:,2),'b', 'LineWidth',2, 'AutoScaleFactor', 0.5)
            end
            hold off
        end
        subplot(3,1,1); axis off; axis equal; axis tight;
        subplot(3,1,2); axis off; axis equal; axis tight;
        subplot(3,1,3); axis off; axis equal; axis tight;
    end

    % DESIGN UPDATE USING FMINCON
    [x,fval,stressMat,exitflag,output,lambda] = solver(x);

    % PRINT RESULTS
    change = max(max(abs(x-xold)));
    disp([' It.: ' sprintf('%4i',loop) ' Obj.: ' sprintf('%10.4f',fval) ...
    ' Vol.: ' sprintf('%6.3f',sum(x.*area_e)/sum(area_e)) ...
    ' maxVMStress: ' sprintf('%6.3e',max(stressMat)) ...
    ' ch.: ' sprintf('%6.3f',change )])

    % PLOT DENSITIES
    set(plot1, 'CData', -x);
    set(plot2, 'CData', stressMat)
    colormap(ax(1), gray);
    colormap(ax(2), parula);
    
    pause(1e-6);
end

% Post loop thresholding filter to generate a more binary design
%added to show more easily where the final design components should be 6
%Dec 2020
%  x_post_filter = x;
%   criteria = false;
%   while criteria == false
%       for threshold = 0:0.01:1;
%           for i = 1:nE
%             if x_post_filter(i) == 1
%                 x_post_filter(i) = 1;
%             elseif x_post_filter(i) > threshold
%                 x_post_filter(i) = 0.8;
%             else
%                 x_post_filter(i) = 0.05;
%             end
%           end
%             if (x_post_filter'*area_e)/sum(area_e) < volfrac
%                 criteria = true;
%                 x = x_post_filter;
%                 break
%             else
%                 x_post_filter = x;
%             end
%       end
%   end

    %Showing the max VM stress on the figure
    figure(1); subplot(2,1,2);
    str1 = "Yield Stress: " + num2str(oy,'%.3e');
    str2 = "Max VM Stress: " + num2str(max(stressMat),'%.3e');
    annotation('textbox',[0.2 0.05 0.2 0.2],'String',{str1,str2},'FitBoxToText','on','HorizontalAlignment','left','VerticalAlignment','middle');

    set(plot1, 'CData', -x);
    colormap(ax(1), gray);
    
    disp(['Final Obj Func Value: ' sprintf('%f',fval) ' exitflag: ' sprintf('%d',exitflag)])
    output
    maxStress = max(stressMat)
    volfrac

%%%%%%%%%% SOLVER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [xnew,fval,stressMat,exitflag,output,lambda] = solver(x)
    % Intialize fmincon
    obj = @objective;
    x0 = x;
    lb = 0.001*ones(size(x));
    ub = 1*ones(size(x));
    nonlcon = @nonlinear;
    A = stress(x);
    b = ones(nE,1)*0;
    % Setting fmincon options
    options = optimoptions(@fmincon,...
    'SpecifyObjectiveGradient',true,...
    'StepTolerance', 1.0000e-10,...
    'MaxFunctionEvaluations',Inf,...
    'MaxIterations',1000,...
    'Algorithm','sqp');
%     % Using fmincon with finite differences
%     options = optimoptions(@fmincon,...
%         'Algorithm','sqp');
    
    % Calling fmincon
    [xnew,fval,exitflag,output,lambda] = fmincon(obj,x0,A,b,[],[],lb,ub,nonlcon,options);
    exitflag
    fval
    output

    % Calculate Stresses for new design xnew
    stressMat = zeros(nE, 1); %Matrix for storing and plotting vonMises stress
    for i = 1:nE
        %Get Ue
        nodes_i = ElementList(i,[2 3 4 5]);
        edof = reshape([nodes_i*2 - 1; nodes_i*2], [], 1);
        Ue = U(edof, 1);
        %Get vonMises stress
        nodePos_i = NodeList(nodes_i, 2:3)';
        [~, ~, stressVMe] = evalElement(nodePos_i, Ue, E, nu); %Calculate stress
        stressMat(i) = stressVMe*x(i)^(1/2);   
    end
end
%%%%%%%%%% OBJECTIVE FUNCTION%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [fx,dx] = objective(x)
    % First need to calculate U(x). K(x) is constant
    % Compute U(x)
    [U, F]=FE(nE,nN,ElementList,x,penal,k_elements,NodeList_FixedDOF,NodeList_Spar1,NodeList_Spar2,NodeList_Surface);
    % Now Compute c(x) and dc(x)
    c = 0.0;
    dc = zeros(nE,1);
    for m = 1:3
        for i =1:nE
            % First pull correct values for Ue(x)
            nodes_i = ElementList(i,[2 3 4 5]);
            edof = reshape([nodes_i*2 - 1; nodes_i*2], [], 1);
            Ue = U(edof, m);
            % Now compute c(x) and dc(x)
            c = c + x(i)^penal*Ue'*k_elements(:,:,i)*Ue;
            % This isnt actually the gradient. More correct to call it
            % sensitivity?
            dc(i) = dc(i)-penal*x(i)^(penal-1)*Ue'*k_elements(:,:,i)*Ue;
        end
    end
    % Now apply a filter to this sensitivity
    [dc] = check(nE,centroid_e,rmin,x,dc);
    % Assign Outputs
    fx = c;
    dx = dc;
end

%%%%%%%%%% LINEAR CONSTRAINT FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [A] = stress(x)
            % Calculate Stresses for new design xnew
            stressMat = zeros(nE, 1); %Matrix for storing and plotting vonMises stress
            for i = 1:nE
                %Get Ue
                nodes_i = ElementList(i,[2 3 4 5]);
                edof = reshape([nodes_i*2 - 1; nodes_i*2], [], 1);
                Ue = U(edof, 1);
                %Get vonMises stress
                nodePos_i = NodeList(nodes_i, 2:3)';
                [~, ~, stressVMe] = evalElement(nodePos_i, Ue, E, nu); %Calculate stress
                stressMat(i) = stressVMe*x(i)^(1/2);   
            end
        A = diag(stressMat-oy);
    end

%%%%%%%%%% NONLINEAR CONSTRAINT FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [c,ceq] = nonlinear(x)
    c = [];
    ceq(1) = sum(x.*area_e)/sum(area_e)-volfrac; % Volume Fraction Constraint
    ceq(2) = sum(x.*ElementList_Active_Bool(:,2))/sum(ElementList_Active_Bool(:,2)) - 1; % Active Elements Constraint
end

%%%%%%%%%% MESH-INDEPENDENCY FILTER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dcn]=check(nE,centroid_e, rmin,x,dc)
    dcn=zeros(nE,1);

    for i = 1:nE
        %Find elements within a range rmin
        %Use centroids of elements
        centroid_i = centroid_e(i,:);
        dCentroid = centroid_e - centroid_i;
        distCentroid = sqrt(dCentroid(:,1).^2 + dCentroid(:,2).^2);
        
        inRange = (distCentroid <= rmin);
     
        %Do the rest of filtering
        fac = rmin - distCentroid(inRange);
        total = sum(fac);
        
        dcn(i) = sum(fac.*x(inRange).*dc(inRange));
        dcn(i) = dcn(i)/(x(i)*total);        
    end
end

%%%%%%%%%% FE-ANALYSIS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [U, F]=FE(nE, nN, ElementList, x,penal, k_elements, NodeList_FixedDOF, NodeList_Spar1, NodeList_Spar2, NodeList_Surface)

    K = sparse(2*nN, 2*nN);
    F = sparse(2*nN,3);
    U = zeros(2*nN,3);

    %Construct stiffness matrix
    for j = 1:nE
        nodes_i = ElementList(j,[2 3 4 5]); %Get nodes for each element
        edof = reshape([nodes_i*2 - 1; nodes_i*2], [], 1); %Get indices of nodes (X and Y) in K/F/U matrices 
        K(edof, edof) = K(edof,edof) + x(j)^penal*k_elements(:,:,j); %Add to total stiffness matrix
    end

    % DEFINE LOADS CASES AND SUPPORTS
    %%%%%%%%%%%%%%%%%%%%%%%%% STEADY STATE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    F(NodeList_Spar1*2,1) = -400; % -Y load on spar 1
    F(NodeList_Spar1*2-1,1) = -80; % -X load on spar 1
    F(NodeList_Spar2*2, 1) = -400; % -Y load on spar 2
    F(NodeList_Spar2*2-1, 1) = -80; % -X load on spar 2
    
    F(zoneLeadTop*2,1) = -1000; %-Y load on leading edge top
    F(zoneTrailTop*2,1) = -1000; %-Y load on trailing edge top
    F(zoneLeadBot*2,1) = -200; %-Y load on leading edge bottom
    F(zoneTrailBot*2,1) = -200; %-Y load on trailing edge bottom
    
    F(zoneLeadTop*2-1,1) = -200; %-X load on leading edge top
    F(zoneTrailTop*2-1,1) = -200; %-X load on trailing edge top
    F(zoneLeadBot*2-1,1) = -40; %-X load on leading edge bottom
    F(zoneTrailBot*2-1,1) = -40; %-X load on trailing edge bottom
    
    %%%%%%%%%%%%%%%%%%%%%%%%%% ACCELERATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    F(NodeList_Spar1*2,2) = -100; % -Y load on spar 1
    F(NodeList_Spar1*2-1,2) = -10; % -X load on spar 1
    F(NodeList_Spar2*2, 2) = -100; % -Y load on spar 2
    F(NodeList_Spar2*2-1, 2) = -10; % -X load on spar 2
    
    F(zoneLeadTop*2,2) = -100; %-Y load on leading edge top
    F(zoneTrailTop*2,2) = -100; %-Y load on trailing edge top
    F(zoneLeadBot*2,2) = -100; %-Y load on leading edge bottom
    F(zoneTrailBot*2,2) = -100; %-Y load on trailing edge bottom
    
    F(zoneLeadTop*2-1,2) = -10; %-X load on leading edge top
    F(zoneTrailTop*2-1,2) = -10; %-X load on trailing edge top
    F(zoneLeadBot*2-1,2) = -10; %-X load on leading edge bottom
    F(zoneTrailBot*2-1,2) = -10; %-X load on trailing edge bottom
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%% BRAKING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    F(NodeList_Spar1*2,3) = -200; % -Y load on spar 1
    F(NodeList_Spar1*2-1,3) = +40; % -X load on spar 1
    F(NodeList_Spar2*2, 3) = -200; % -Y load on spar 2
    F(NodeList_Spar2*2-1, 3) = +40; % -X load on spar 2
    
    F(zoneLeadTop*2,3) = -800; %-Y load on leading edge top
    F(zoneTrailTop*2,3) = -800; %-Y load on trailing edge top
    F(zoneLeadBot*2,3) = -400; %-Y load on leading edge bottom
    F(zoneTrailBot*2,3) = -400; %-Y load on trailing edge bottom
    
    F(zoneLeadTop*2-1,3) = 100; %+X load on leading edge top
    F(zoneTrailTop*2-1,3) = 100; %+X load on trailing edge top
    F(zoneLeadBot*2-1,3) = 100; %+X load on leading edge bottom
    F(zoneTrailBot*2-1,3) = 100; %+X load on trailing edge bottom
    
    fixeddofs   = reshape([NodeList_FixedDOF*2 - 1; NodeList_FixedDOF*2], [], 1);
    alldofs     = [1:2*nN];
    freedofs    = setdiff(alldofs,fixeddofs);
    % SOLVING
    U(freedofs,:) = K(freedofs,freedofs) \ F(freedofs,:);      
    U(fixeddofs,:)= 0;
end

%%%%%%%%%% ELEMENT STIFFNESS MATRIX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function k = calcStiffness(nodePos, E, v)
    % k: individual stiffness matrix for generic quadrilateral element (8x8 matrix)
    % nodePos: X, Y coordinates of nodes (2x4 matrix)...
    %          [X1, X2, X3, X4]
    %          [Y1, Y2, Y3, Y4]
    % Nodes numbered counterclockwise starting with bottom left node         
    % 4-------3
    % |       |
    % |       |
    % 1-------2
    %
    % E: Young's modulus (scalar)
    % v: Poisson ratio (scalar)

    x1 = nodePos(1,1);
    x2 = nodePos(1,2);
    x3 = nodePos(1,3);
    x4 = nodePos(1,4);
    y1 = nodePos(2,1);
    y2 = nodePos(2,2);
    y3 = nodePos(2,3);
    y4 = nodePos(2,4);

    % Constituative Matrix for Hookes Law for Plane Stress 
    Emat = E/(1-v^2)*[1, v, 0;...
                      v, 1, 0;...
                      0, 0, (1-v)/2];
    
    % Determinant of Jacobian Matrix for Strain Displacement Matrix
    % (derived with MATLAB Symbolic toolbox for our geometry)
    detJacEq = @(e,n,x1,x2,x3,x4,y1,y2,y3,y4)(x1*y2)/8.0-(x2*y1)/8.0-(x1*y4)/8.0+(x2*y3)/8.0-(x3*y2)/8.0+(x4*y1)/8.0+(x3*y4)/8.0-(x4*y3)/8.0-(e*x1*y3)/8.0+(e*x3*y1)/8.0+(e*x1*y4)/8.0+(e*x2*y3)/8.0-(e*x3*y2)/8.0-(e*x4*y1)/8.0-(e*x2*y4)/8.0+(e*x4*y2)/8.0-(n*x1*y2)/8.0+(n*x2*y1)/8.0+(n*x1*y3)/8.0-(n*x3*y1)/8.0-(n*x2*y4)/8.0+(n*x4*y2)/8.0+(n*x3*y4)/8.0-(n*x4*y3)/8.0;

    % Strain Displacement Matrix
    % (derived with MATLAB Symbolic toolbox for our geometry)
    BEq = @(e,n,x1,x2,x3,x4,y1,y2,y3,y4)reshape([((n/4.0-1.0/4.0)*(y1+y2-y3-y4-e*y1+e*y2-e*y3+e*y4)*-2.0)/(x1*y2-x2*y1-x1*y4+x2*y3-x3*y2+x4*y1+x3*y4-x4*y3-e*x1*y3+e*x3*y1+e*x1*y4+e*x2*y3-e*x3*y2-e*x4*y1-e*x2*y4+e*x4*y2-n*x1*y2+n*x2*y1+n*x1*y3-n*x3*y1-n*x2*y4+n*x4*y2+n*x3*y4-n*x4*y3)+((e/4.0-1.0/4.0)*(y1-y2-y3+y4-n*y1+n*y2-n*y3+n*y4)*2.0)/(x1*y2-x2*y1-x1*y4+x2*y3-x3*y2+x4*y1+x3*y4-x4*y3-e*x1*y3+e*x3*y1+e*x1*y4+e*x2*y3-e*x3*y2-e*x4*y1-e*x2*y4+e*x4*y2-n*x1*y2+n*x2*y1+n*x1*y3-n*x3*y1-n*x2*y4+n*x4*y2+n*x3*y4-n*x4*y3),0.0,((n/4.0-1.0/4.0)*(x1+x2-x3-x4-e*x1+e*x2-e*x3+e*x4)*2.0)/(x1*y2-x2*y1-x1*y4+x2*y3-x3*y2+x4*y1+x3*y4-x4*y3-e*x1*y3+e*x3*y1+e*x1*y4+e*x2*y3-e*x3*y2-e*x4*y1-e*x2*y4+e*x4*y2-n*x1*y2+n*x2*y1+n*x1*y3-n*x3*y1-n*x2*y4+n*x4*y2+n*x3*y4-n*x4*y3)-((e/4.0-1.0/4.0)*(x1-x2-x3+x4-n*x1+n*x2-n*x3+n*x4)*2.0)/(x1*y2-x2*y1-x1*y4+x2*y3-x3*y2+x4*y1+x3*y4-x4*y3-e*x1*y3+e*x3*y1+e*x1*y4+e*x2*y3-e*x3*y2-e*x4*y1-e*x2*y4+e*x4*y2-n*x1*y2+n*x2*y1+n*x1*y3-n*x3*y1-n*x2*y4+n*x4*y2+n*x3*y4-n*x4*y3),0.0,((n/4.0-1.0/4.0)*(x1+x2-x3-x4-e*x1+e*x2-e*x3+e*x4)*2.0)/(x1*y2-x2*y1-x1*y4+x2*y3-x3*y2+x4*y1+x3*y4-x4*y3-e*x1*y3+e*x3*y1+e*x1*y4+e*x2*y3-e*x3*y2-e*x4*y1-e*x2*y4+e*x4*y2-n*x1*y2+n*x2*y1+n*x1*y3-n*x3*y1-n*x2*y4+n*x4*y2+n*x3*y4-n*x4*y3)-((e/4.0-1.0/4.0)*(x1-x2-x3+x4-n*x1+n*x2-n*x3+n*x4)*2.0)/(x1*y2-x2*y1-x1*y4+x2*y3-x3*y2+x4*y1+x3*y4-x4*y3-e*x1*y3+e*x3*y1+e*x1*y4+e*x2*y3-e*x3*y2-e*x4*y1-e*x2*y4+e*x4*y2-n*x1*y2+n*x2*y1+n*x1*y3-n*x3*y1-n*x2*y4+n*x4*y2+n*x3*y4-n*x4*y3),((n/4.0-1.0/4.0)*(y1+y2-y3-y4-e*y1+e*y2-e*y3+e*y4)*-2.0)/(x1*y2-x2*y1-x1*y4+x2*y3-x3*y2+x4*y1+x3*y4-x4*y3-e*x1*y3+e*x3*y1+e*x1*y4+e*x2*y3-e*x3*y2-e*x4*y1-e*x2*y4+e*x4*y2-n*x1*y2+n*x2*y1+n*x1*y3-n*x3*y1-n*x2*y4+n*x4*y2+n*x3*y4-n*x4*y3)+((e/4.0-1.0/4.0)*(y1-y2-y3+y4-n*y1+n*y2-n*y3+n*y4)*2.0)/(x1*y2-x2*y1-x1*y4+x2*y3-x3*y2+x4*y1+x3*y4-x4*y3-e*x1*y3+e*x3*y1+e*x1*y4+e*x2*y3-e*x3*y2-e*x4*y1-e*x2*y4+e*x4*y2-n*x1*y2+n*x2*y1+n*x1*y3-n*x3*y1-n*x2*y4+n*x4*y2+n*x3*y4-n*x4*y3),((n/4.0-1.0/4.0)*(y1+y2-y3-y4-e*y1+e*y2-e*y3+e*y4)*2.0)/(x1*y2-x2*y1-x1*y4+x2*y3-x3*y2+x4*y1+x3*y4-x4*y3-e*x1*y3+e*x3*y1+e*x1*y4+e*x2*y3-e*x3*y2-e*x4*y1-e*x2*y4+e*x4*y2-n*x1*y2+n*x2*y1+n*x1*y3-n*x3*y1-n*x2*y4+n*x4*y2+n*x3*y4-n*x4*y3)-((e/4.0+1.0/4.0)*(y1-y2-y3+y4-n*y1+n*y2-n*y3+n*y4)*2.0)/(x1*y2-x2*y1-x1*y4+x2*y3-x3*y2+x4*y1+x3*y4-x4*y3-e*x1*y3+e*x3*y1+e*x1*y4+e*x2*y3-e*x3*y2-e*x4*y1-e*x2*y4+e*x4*y2-n*x1*y2+n*x2*y1+n*x1*y3-n*x3*y1-n*x2*y4+n*x4*y2+n*x3*y4-n*x4*y3),0.0,((n/4.0-1.0/4.0)*(x1+x2-x3-x4-e*x1+e*x2-e*x3+e*x4)*-2.0)/(x1*y2-x2*y1-x1*y4+x2*y3-x3*y2+x4*y1+x3*y4-x4*y3-e*x1*y3+e*x3*y1+e*x1*y4+e*x2*y3-e*x3*y2-e*x4*y1-e*x2*y4+e*x4*y2-n*x1*y2+n*x2*y1+n*x1*y3-n*x3*y1-n*x2*y4+n*x4*y2+n*x3*y4-n*x4*y3)+((e/4.0+1.0/4.0)*(x1-x2-x3+x4-n*x1+n*x2-n*x3+n*x4)*2.0)/(x1*y2-x2*y1-x1*y4+x2*y3-x3*y2+x4*y1+x3*y4-x4*y3-e*x1*y3+e*x3*y1+e*x1*y4+e*x2*y3-e*x3*y2-e*x4*y1-e*x2*y4+e*x4*y2-n*x1*y2+n*x2*y1+n*x1*y3-n*x3*y1-n*x2*y4+n*x4*y2+n*x3*y4-n*x4*y3),0.0,((n/4.0-1.0/4.0)*(x1+x2-x3-x4-e*x1+e*x2-e*x3+e*x4)*-2.0)/(x1*y2-x2*y1-x1*y4+x2*y3-x3*y2+x4*y1+x3*y4-x4*y3-e*x1*y3+e*x3*y1+e*x1*y4+e*x2*y3-e*x3*y2-e*x4*y1-e*x2*y4+e*x4*y2-n*x1*y2+n*x2*y1+n*x1*y3-n*x3*y1-n*x2*y4+n*x4*y2+n*x3*y4-n*x4*y3)+((e/4.0+1.0/4.0)*(x1-x2-x3+x4-n*x1+n*x2-n*x3+n*x4)*2.0)/(x1*y2-x2*y1-x1*y4+x2*y3-x3*y2+x4*y1+x3*y4-x4*y3-e*x1*y3+e*x3*y1+e*x1*y4+e*x2*y3-e*x3*y2-e*x4*y1-e*x2*y4+e*x4*y2-n*x1*y2+n*x2*y1+n*x1*y3-n*x3*y1-n*x2*y4+n*x4*y2+n*x3*y4-n*x4*y3),((n/4.0-1.0/4.0)*(y1+y2-y3-y4-e*y1+e*y2-e*y3+e*y4)*2.0)/(x1*y2-x2*y1-x1*y4+x2*y3-x3*y2+x4*y1+x3*y4-x4*y3-e*x1*y3+e*x3*y1+e*x1*y4+e*x2*y3-e*x3*y2-e*x4*y1-e*x2*y4+e*x4*y2-n*x1*y2+n*x2*y1+n*x1*y3-n*x3*y1-n*x2*y4+n*x4*y2+n*x3*y4-n*x4*y3)-((e/4.0+1.0/4.0)*(y1-y2-y3+y4-n*y1+n*y2-n*y3+n*y4)*2.0)/(x1*y2-x2*y1-x1*y4+x2*y3-x3*y2+x4*y1+x3*y4-x4*y3-e*x1*y3+e*x3*y1+e*x1*y4+e*x2*y3-e*x3*y2-e*x4*y1-e*x2*y4+e*x4*y2-n*x1*y2+n*x2*y1+n*x1*y3-n*x3*y1-n*x2*y4+n*x4*y2+n*x3*y4-n*x4*y3),((n/4.0+1.0/4.0)*(y1+y2-y3-y4-e*y1+e*y2-e*y3+e*y4)*-2.0)/(x1*y2-x2*y1-x1*y4+x2*y3-x3*y2+x4*y1+x3*y4-x4*y3-e*x1*y3+e*x3*y1+e*x1*y4+e*x2*y3-e*x3*y2-e*x4*y1-e*x2*y4+e*x4*y2-n*x1*y2+n*x2*y1+n*x1*y3-n*x3*y1-n*x2*y4+n*x4*y2+n*x3*y4-n*x4*y3)+((e/4.0+1.0/4.0)*(y1-y2-y3+y4-n*y1+n*y2-n*y3+n*y4)*2.0)/(x1*y2-x2*y1-x1*y4+x2*y3-x3*y2+x4*y1+x3*y4-x4*y3-e*x1*y3+e*x3*y1+e*x1*y4+e*x2*y3-e*x3*y2-e*x4*y1-e*x2*y4+e*x4*y2-n*x1*y2+n*x2*y1+n*x1*y3-n*x3*y1-n*x2*y4+n*x4*y2+n*x3*y4-n*x4*y3),0.0,((n/4.0+1.0/4.0)*(x1+x2-x3-x4-e*x1+e*x2-e*x3+e*x4)*2.0)/(x1*y2-x2*y1-x1*y4+x2*y3-x3*y2+x4*y1+x3*y4-x4*y3-e*x1*y3+e*x3*y1+e*x1*y4+e*x2*y3-e*x3*y2-e*x4*y1-e*x2*y4+e*x4*y2-n*x1*y2+n*x2*y1+n*x1*y3-n*x3*y1-n*x2*y4+n*x4*y2+n*x3*y4-n*x4*y3)-((e/4.0+1.0/4.0)*(x1-x2-x3+x4-n*x1+n*x2-n*x3+n*x4)*2.0)/(x1*y2-x2*y1-x1*y4+x2*y3-x3*y2+x4*y1+x3*y4-x4*y3-e*x1*y3+e*x3*y1+e*x1*y4+e*x2*y3-e*x3*y2-e*x4*y1-e*x2*y4+e*x4*y2-n*x1*y2+n*x2*y1+n*x1*y3-n*x3*y1-n*x2*y4+n*x4*y2+n*x3*y4-n*x4*y3),0.0,((n/4.0+1.0/4.0)*(x1+x2-x3-x4-e*x1+e*x2-e*x3+e*x4)*2.0)/(x1*y2-x2*y1-x1*y4+x2*y3-x3*y2+x4*y1+x3*y4-x4*y3-e*x1*y3+e*x3*y1+e*x1*y4+e*x2*y3-e*x3*y2-e*x4*y1-e*x2*y4+e*x4*y2-n*x1*y2+n*x2*y1+n*x1*y3-n*x3*y1-n*x2*y4+n*x4*y2+n*x3*y4-n*x4*y3)-((e/4.0+1.0/4.0)*(x1-x2-x3+x4-n*x1+n*x2-n*x3+n*x4)*2.0)/(x1*y2-x2*y1-x1*y4+x2*y3-x3*y2+x4*y1+x3*y4-x4*y3-e*x1*y3+e*x3*y1+e*x1*y4+e*x2*y3-e*x3*y2-e*x4*y1-e*x2*y4+e*x4*y2-n*x1*y2+n*x2*y1+n*x1*y3-n*x3*y1-n*x2*y4+n*x4*y2+n*x3*y4-n*x4*y3),((n/4.0+1.0/4.0)*(y1+y2-y3-y4-e*y1+e*y2-e*y3+e*y4)*-2.0)/(x1*y2-x2*y1-x1*y4+x2*y3-x3*y2+x4*y1+x3*y4-x4*y3-e*x1*y3+e*x3*y1+e*x1*y4+e*x2*y3-e*x3*y2-e*x4*y1-e*x2*y4+e*x4*y2-n*x1*y2+n*x2*y1+n*x1*y3-n*x3*y1-n*x2*y4+n*x4*y2+n*x3*y4-n*x4*y3)+((e/4.0+1.0/4.0)*(y1-y2-y3+y4-n*y1+n*y2-n*y3+n*y4)*2.0)/(x1*y2-x2*y1-x1*y4+x2*y3-x3*y2+x4*y1+x3*y4-x4*y3-e*x1*y3+e*x3*y1+e*x1*y4+e*x2*y3-e*x3*y2-e*x4*y1-e*x2*y4+e*x4*y2-n*x1*y2+n*x2*y1+n*x1*y3-n*x3*y1-n*x2*y4+n*x4*y2+n*x3*y4-n*x4*y3),((n/4.0+1.0/4.0)*(y1+y2-y3-y4-e*y1+e*y2-e*y3+e*y4)*2.0)/(x1*y2-x2*y1-x1*y4+x2*y3-x3*y2+x4*y1+x3*y4-x4*y3-e*x1*y3+e*x3*y1+e*x1*y4+e*x2*y3-e*x3*y2-e*x4*y1-e*x2*y4+e*x4*y2-n*x1*y2+n*x2*y1+n*x1*y3-n*x3*y1-n*x2*y4+n*x4*y2+n*x3*y4-n*x4*y3)-((e/4.0-1.0/4.0)*(y1-y2-y3+y4-n*y1+n*y2-n*y3+n*y4)*2.0)/(x1*y2-x2*y1-x1*y4+x2*y3-x3*y2+x4*y1+x3*y4-x4*y3-e*x1*y3+e*x3*y1+e*x1*y4+e*x2*y3-e*x3*y2-e*x4*y1-e*x2*y4+e*x4*y2-n*x1*y2+n*x2*y1+n*x1*y3-n*x3*y1-n*x2*y4+n*x4*y2+n*x3*y4-n*x4*y3),0.0,((n/4.0+1.0/4.0)*(x1+x2-x3-x4-e*x1+e*x2-e*x3+e*x4)*-2.0)/(x1*y2-x2*y1-x1*y4+x2*y3-x3*y2+x4*y1+x3*y4-x4*y3-e*x1*y3+e*x3*y1+e*x1*y4+e*x2*y3-e*x3*y2-e*x4*y1-e*x2*y4+e*x4*y2-n*x1*y2+n*x2*y1+n*x1*y3-n*x3*y1-n*x2*y4+n*x4*y2+n*x3*y4-n*x4*y3)+((e/4.0-1.0/4.0)*(x1-x2-x3+x4-n*x1+n*x2-n*x3+n*x4)*2.0)/(x1*y2-x2*y1-x1*y4+x2*y3-x3*y2+x4*y1+x3*y4-x4*y3-e*x1*y3+e*x3*y1+e*x1*y4+e*x2*y3-e*x3*y2-e*x4*y1-e*x2*y4+e*x4*y2-n*x1*y2+n*x2*y1+n*x1*y3-n*x3*y1-n*x2*y4+n*x4*y2+n*x3*y4-n*x4*y3),0.0,((n/4.0+1.0/4.0)*(x1+x2-x3-x4-e*x1+e*x2-e*x3+e*x4)*-2.0)/(x1*y2-x2*y1-x1*y4+x2*y3-x3*y2+x4*y1+x3*y4-x4*y3-e*x1*y3+e*x3*y1+e*x1*y4+e*x2*y3-e*x3*y2-e*x4*y1-e*x2*y4+e*x4*y2-n*x1*y2+n*x2*y1+n*x1*y3-n*x3*y1-n*x2*y4+n*x4*y2+n*x3*y4-n*x4*y3)+((e/4.0-1.0/4.0)*(x1-x2-x3+x4-n*x1+n*x2-n*x3+n*x4)*2.0)/(x1*y2-x2*y1-x1*y4+x2*y3-x3*y2+x4*y1+x3*y4-x4*y3-e*x1*y3+e*x3*y1+e*x1*y4+e*x2*y3-e*x3*y2-e*x4*y1-e*x2*y4+e*x4*y2-n*x1*y2+n*x2*y1+n*x1*y3-n*x3*y1-n*x2*y4+n*x4*y2+n*x3*y4-n*x4*y3),((n/4.0+1.0/4.0)*(y1+y2-y3-y4-e*y1+e*y2-e*y3+e*y4)*2.0)/(x1*y2-x2*y1-x1*y4+x2*y3-x3*y2+x4*y1+x3*y4-x4*y3-e*x1*y3+e*x3*y1+e*x1*y4+e*x2*y3-e*x3*y2-e*x4*y1-e*x2*y4+e*x4*y2-n*x1*y2+n*x2*y1+n*x1*y3-n*x3*y1-n*x2*y4+n*x4*y2+n*x3*y4-n*x4*y3)-((e/4.0-1.0/4.0)*(y1-y2-y3+y4-n*y1+n*y2-n*y3+n*y4)*2.0)/(x1*y2-x2*y1-x1*y4+x2*y3-x3*y2+x4*y1+x3*y4-x4*y3-e*x1*y3+e*x3*y1+e*x1*y4+e*x2*y3-e*x3*y2-e*x4*y1-e*x2*y4+e*x4*y2-n*x1*y2+n*x2*y1+n*x1*y3-n*x3*y1-n*x2*y4+n*x4*y2+n*x3*y4-n*x4*y3)],[3,8]);

    k = 0;
    % Calculate Stiffness Matrix using Gauss Quadrature
    for i = [-1/sqrt(3), 1/sqrt(3)]
        for j = [-1/sqrt(3), 1/sqrt(3)]
            detJac = detJacEq(i, j, x1, x2, x3, x4, y1, y2, y3, y4);
            B = BEq(i, j, x1, x2, x3, x4, y1, y2, y3, y4);

            BEB = B'*Emat*B*detJac;
            k = k + BEB;
        end
    end
end

function [strain, stress, stressVM] = evalElement(nodePos, nodeQ, E, v)
    % 
    % nodePos: X, Y coordinates of nodes (2x4 matrix)...
    %          [X1, X2, X3, X4]
    %          [Y1, Y2, Y3, Y4]
    % Nodes numbered counterclockwise starting with bottom left node         
    % 4-------3
    % |       |
    % |       |
    % 1-------2
    %
    % nodeQ: X,Y displacements (8x1 matrix)
    %        [dX1; dY1; dX2; dY2; ...]
    %
    % E: Young's modulus (scalar)
    % v: Poisson ratio (scalar)

    x1 = nodePos(1,1);
    x2 = nodePos(1,2);
    x3 = nodePos(1,3);
    x4 = nodePos(1,4);
    y1 = nodePos(2,1);
    y2 = nodePos(2,2);
    y3 = nodePos(2,3);
    y4 = nodePos(2,4);

    % Constituative matrix for ookes law for plane stress 
    Emat = E/(1-v^2)*[1, v, 0;...
                      v, 1, 0;...
                      0, 0, (1-v)/2];

    % Strain Displacement Matrix Evaluated at Isoparametric Center of Element            
    B = [ (y1 + y2 - y3 - y4)/(2*(x1*y2 - x2*y1 - x1*y4 + x2*y3 - x3*y2 + x4*y1 + x3*y4 - x4*y3)) - (y1 - y2 - y3 + y4)/(2*(x1*y2 - x2*y1 - x1*y4 + x2*y3 - x3*y2 + x4*y1 + x3*y4 - x4*y3)), 0, - (y1 + y2 - y3 - y4)/(2*(x1*y2 - x2*y1 - x1*y4 + x2*y3 - x3*y2 + x4*y1 + x3*y4 - x4*y3)) - (y1 - y2 - y3 + y4)/(2*(x1*y2 - x2*y1 - x1*y4 + x2*y3 - x3*y2 + x4*y1 + x3*y4 - x4*y3)), 0, (y1 - y2 - y3 + y4)/(2*(x1*y2 - x2*y1 - x1*y4 + x2*y3 - x3*y2 + x4*y1 + x3*y4 - x4*y3)) - (y1 + y2 - y3 - y4)/(2*(x1*y2 - x2*y1 - x1*y4 + x2*y3 - x3*y2 + x4*y1 + x3*y4 - x4*y3)), 0,   (y1 + y2 - y3 - y4)/(2*(x1*y2 - x2*y1 - x1*y4 + x2*y3 - x3*y2 + x4*y1 + x3*y4 - x4*y3)) + (y1 - y2 - y3 + y4)/(2*(x1*y2 - x2*y1 - x1*y4 + x2*y3 - x3*y2 + x4*y1 + x3*y4 - x4*y3)), 0;...
        0, (x1 - x2 - x3 + x4)/(2*(x1*y2 - x2*y1 - x1*y4 + x2*y3 - x3*y2 + x4*y1 + x3*y4 - x4*y3)) - (x1 + x2 - x3 - x4)/(2*(x1*y2 - x2*y1 - x1*y4 + x2*y3 - x3*y2 + x4*y1 + x3*y4 - x4*y3)), 0,   (x1 + x2 - x3 - x4)/(2*(x1*y2 - x2*y1 - x1*y4 + x2*y3 - x3*y2 + x4*y1 + x3*y4 - x4*y3)) + (x1 - x2 - x3 + x4)/(2*(x1*y2 - x2*y1 - x1*y4 + x2*y3 - x3*y2 + x4*y1 + x3*y4 - x4*y3)), 0, (x1 + x2 - x3 - x4)/(2*(x1*y2 - x2*y1 - x1*y4 + x2*y3 - x3*y2 + x4*y1 + x3*y4 - x4*y3)) - (x1 - x2 - x3 + x4)/(2*(x1*y2 - x2*y1 - x1*y4 + x2*y3 - x3*y2 + x4*y1 + x3*y4 - x4*y3)), 0, - (x1 + x2 - x3 - x4)/(2*(x1*y2 - x2*y1 - x1*y4 + x2*y3 - x3*y2 + x4*y1 + x3*y4 - x4*y3)) - (x1 - x2 - x3 + x4)/(2*(x1*y2 - x2*y1 - x1*y4 + x2*y3 - x3*y2 + x4*y1 + x3*y4 - x4*y3));...
        (x1 - x2 - x3 + x4)/(2*(x1*y2 - x2*y1 - x1*y4 + x2*y3 - x3*y2 + x4*y1 + x3*y4 - x4*y3)) - (x1 + x2 - x3 - x4)/(2*(x1*y2 - x2*y1 - x1*y4 + x2*y3 - x3*y2 + x4*y1 + x3*y4 - x4*y3)), (y1 + y2 - y3 - y4)/(2*(x1*y2 - x2*y1 - x1*y4 + x2*y3 - x3*y2 + x4*y1 + x3*y4 - x4*y3)) - (y1 - y2 - y3 + y4)/(2*(x1*y2 - x2*y1 - x1*y4 + x2*y3 - x3*y2 + x4*y1 + x3*y4 - x4*y3)),   (x1 + x2 - x3 - x4)/(2*(x1*y2 - x2*y1 - x1*y4 + x2*y3 - x3*y2 + x4*y1 + x3*y4 - x4*y3)) + (x1 - x2 - x3 + x4)/(2*(x1*y2 - x2*y1 - x1*y4 + x2*y3 - x3*y2 + x4*y1 + x3*y4 - x4*y3)), - (y1 + y2 - y3 - y4)/(2*(x1*y2 - x2*y1 - x1*y4 + x2*y3 - x3*y2 + x4*y1 + x3*y4 - x4*y3)) - (y1 - y2 - y3 + y4)/(2*(x1*y2 - x2*y1 - x1*y4 + x2*y3 - x3*y2 + x4*y1 + x3*y4 - x4*y3)), (x1 + x2 - x3 - x4)/(2*(x1*y2 - x2*y1 - x1*y4 + x2*y3 - x3*y2 + x4*y1 + x3*y4 - x4*y3)) - (x1 - x2 - x3 + x4)/(2*(x1*y2 - x2*y1 - x1*y4 + x2*y3 - x3*y2 + x4*y1 + x3*y4 - x4*y3)), (y1 - y2 - y3 + y4)/(2*(x1*y2 - x2*y1 - x1*y4 + x2*y3 - x3*y2 + x4*y1 + x3*y4 - x4*y3)) - (y1 + y2 - y3 - y4)/(2*(x1*y2 - x2*y1 - x1*y4 + x2*y3 - x3*y2 + x4*y1 + x3*y4 - x4*y3)), - (x1 + x2 - x3 - x4)/(2*(x1*y2 - x2*y1 - x1*y4 + x2*y3 - x3*y2 + x4*y1 + x3*y4 - x4*y3)) - (x1 - x2 - x3 + x4)/(2*(x1*y2 - x2*y1 - x1*y4 + x2*y3 - x3*y2 + x4*y1 + x3*y4 - x4*y3)),   (y1 + y2 - y3 - y4)/(2*(x1*y2 - x2*y1 - x1*y4 + x2*y3 - x3*y2 + x4*y1 + x3*y4 - x4*y3)) + (y1 - y2 - y3 + y4)/(2*(x1*y2 - x2*y1 - x1*y4 + x2*y3 - x3*y2 + x4*y1 + x3*y4 - x4*y3))];

    %Strain-dsiplacement relationship
    strain = B * nodeQ;

    %Hookes Law
    stress = Emat*strain;

    %von-Mises Stress
    stressVM = sqrt(stress(1)^2 - stress(1)*stress(2) + stress(2)^2 + 3*stress(3)^2);
end           

function [areas, centroid] = calcArea(ElementList, NodeList)

    %Get nodes for each element
    nodes = ElementList(:,[2 3 4 5 2]);
    nodesX = NodeList(:,2);
    nodesY = NodeList(:,3);
    
    %Get x and y coordinates for each node for all elements
    xx = nodesX(nodes);
    yy = nodesY(nodes);
    
    %Calculate lengths of sides of elements
    dx = diff(xx,1,2);
    dy = diff(yy,1,2);
    L = sqrt(dx.^2 + dy.^2);
    
    %Calculate length of diagonal of element
    dxMid = diff(xx(:, [1,3]), 1, 2);
    dyMid = diff(yy(:, [1,3]), 1, 2);
    LMid = sqrt(dxMid.^2 + dyMid.^2);
    
    %Use Herons formula for each half of quadrilateral
    s1 = (L(:,1) + L(:,2) + LMid)/2;
    s2 = (L(:,3) + L(:,4) + LMid)/2;
    
    area1 = sqrt(s1.*(s1-L(:,1)).*(s1-L(:,2)).*(s1-LMid));
    area2 = sqrt(s2.*(s2-L(:,3)).*(s2-L(:,4)).*(s2-LMid));
    
    areas = area1 + area2;
    
    %use someone else's modified code to get centroids of each element
    [Geom, ~] = polygeom(xx, yy);
    centroid = Geom(:,2:3);    
end
% End of entire function
end