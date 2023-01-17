clear all
close all
clc
%
% -------- Enter inputs -------------
N1min = 1;                % Min dwell time of Healthy mode (H) after STATE poisoning recovery
N1max = 0;                % Max dwell time of STATE Attack mode (A) 
N2min = 1;                % Min dwell time of Healthy mode (H) after INPUT poisoning recovery
N2max = 0;               % Max dwell time of INPUT Attack mode (A)
epsilon = 1e-5;           % Accuracy for Set computations
max_iter = 120;           % Maximum iterations for backward reachability
% 
% -------- End of inputs ------------
% nominal system 
A = [0.9 0.1;0.1 0.8];
B = [0.1; 0];
% LQR control
K = dlqr(A,B,10*eye(2),.1);
Ac1 = A-B*K; % nominal closed loop
% box constraints
u_max = 1;
state_max = 1;
%
% state-constraint set
Z0 = Polyhedron([eye(2);-eye(2);-K;K],...
                [state_max*ones(2,1);state_max*ones(2,1);u_max;u_max]);
% disturbance set
% attack set
au_lim = .3;
ay_lim = .3;
Au = Polyhedron([1;-1],[au_lim;au_lim]);
Ay = Polyhedron([1;-1],[ay_lim;ay_lim]);
%
% construction of system matrices
A0 = A-B*K;
B1 = -B*K*[1;0]; % attack on sensor1: x(t+1) = A0x +B1*ay
B2 = -B*K*[0;1]; % attack on sensor2: x(t+1) = A0x +B2*ay
B3 = -B*K*[1;1]; % attack on sensor12: x(t+1) = A0x +B3*ay
B4 = B;          % attack on actuator: x(t+1) = A0x +B*au
%
% construction of incidence matrix
%
if N1min == 0 && N1max == 0 || N2min == 0 && N2max == 0
    fprintf('Invalid Attack pattern!!!\n')
    return
elseif N1min > 1 && N1max == 0
    N1min = 1;
elseif N2min > 1  && N2max == 0
    N2min = 1;
else
    fprintf('Valid attack pattern has been given\n')
end
%
% -------------- Construction of graphs 
% incidence matrices
Ix = incid_mat(N1min, N1max);
Iu = incid_mat(N2min, N2max);
% edge labels
labelsx = make_labelsxu(Ix);
labelsu = make_labelsxu(Iu);
%
Ixu = [kron(ones(length(Iu(:,1)),1),Ix) kron(Iu, ones(length(Ix(:,1)),1))];
Ixu = [Ixu(:,1) Ixu(:,4) Ixu(:,2) Ixu(:,5) Ixu(:,3) Ixu(:,6)];
%
Ixu_sort = sortrows(Ixu);
%
ids_nodes = [1 Ixu_sort(1,1:2)];
for i = 2:length(Ixu_sort(:,1))
    if Ixu_sort(i,1) ~= Ixu_sort(i-1,1) || Ixu_sort(i,2) ~= Ixu_sort(i-1,2) 
        ids_nodes = [ids_nodes;ids_nodes(length(ids_nodes(:,1)),1)+1 Ixu_sort(i,1:2)];
    end
end
%
multi_edges = [];
Itot_labeled = [];
for i = 1:length(Ixu_sort(:,1))
    for j = 1:length(ids_nodes(:,1))
        if ids_nodes(j,2) == Ixu_sort(i,1) && ids_nodes(j,3) == Ixu_sort(i,2)
            multi_edges(i,1) = ids_nodes(j,1);
            Itot_labeled(i,1) = multi_edges(i,1);
            Itot_labeled(i,3) = Ixu_sort(i,5);
        end
        if ids_nodes(j,2) == Ixu_sort(i,3) && ids_nodes(j,3) == Ixu_sort(i,4)
            multi_edges(i,2) = ids_nodes(j,1);
            Itot_labeled(i,2) = multi_edges(i,2);
            Itot_labeled(i,4) = Ixu_sort(i,6);
        end
    end
end
labels_tot_xu = make_labels_tot(Itot_labeled);

% -------------- Plots of Graphs ------------------------------------
figure(1)
Gx = digraph(Ix(:,1), Ix(:,2))
Px = plot(Gx, 'MarkerSize', 10, 'LineWidth', 2, 'ArrowSize', 10, 'NodeFontSize',10)
Gx.Edges.labels = labelsx;
Px.EdgeLabel = Gx.Edges.labels
title('State poisoning attack pattern')
%
figure(2)
Gu = digraph(Iu(:,1), Iu(:,2))
Pu = plot(Gu, 'MarkerSize', 10, 'LineWidth', 2, 'ArrowSize', 10, 'NodeFontSize',10)
Gu.Edges.labels = labelsu;
Pu.EdgeLabel = Gu.Edges.labels
title('Input poisoning attack pattern')
%
figure(3)
Gxu = digraph(Itot_labeled(:,1), Itot_labeled(:,2))
Pxu = plot(Gxu, 'MarkerSize', 10, 'LineWidth', 2, 'ArrowSize', 10, 'NodeFontSize',10)
Gxu.Edges.labels = labels_tot_xu;
Pxu.EdgeLabel = Gxu.Edges.labels
title('State-Input poisoning attack pattern')
% ------- Construction of Multi-Sets and Maximal Safe set -----------------
nodes = max(Itot_labeled(:,1));
rows = length(Itot_labeled(:,1));
for i = 1:nodes
    S{i}(1) = Z0;
    G{i,1} = Z0.A;
    w{i,1} = Z0.b;
end
%
check = 0;
convergence = 0;
i = 1;
% 1st label refers to BI on sensor 1, 2nd label refers to DoS on sensor 2
while convergence == 0 && i <= max_iter
    for k = 1:nodes
        G{k,i+1} = Z0.A;
        w{k,i+1} = Z0.b;
    end
    for j = 1:rows
        if Itot_labeled(j,3) == 1 && Itot_labeled(j,4) == 1 % attack-free
            G{Itot_labeled(j,1),i+1} = [G{Itot_labeled(j,1),i+1}; G{Itot_labeled(j,2),i}*A0];
            w{Itot_labeled(j,1),i+1} = [w{Itot_labeled(j,1),i+1}; w{Itot_labeled(j,2),i}];
        elseif Itot_labeled(j,3) == 2 && Itot_labeled(j,4) == 1 % attack on sensor
            G{Itot_labeled(j,1),i+1} = [G{Itot_labeled(j,1),i+1}; G{Itot_labeled(j,2),i}*A0];
            w{Itot_labeled(j,1),i+1} = [w{Itot_labeled(j,1),i+1}; w{Itot_labeled(j,2),i}-maxH(G{Itot_labeled(j,2),i}, B1,Ay)];
        elseif Itot_labeled(j,3) == 1 && Itot_labeled(j,4) == 2 % attack on actuator
            G{Itot_labeled(j,1),i+1} = [G{Itot_labeled(j,1),i+1}; G{Itot_labeled(j,2),i}*A0];
            w{Itot_labeled(j,1),i+1} = [w{Itot_labeled(j,1),i+1}; w{Itot_labeled(j,2),i}-maxH(G{Itot_labeled(j,2),i}, B4,Au)];
        else % attack on both sensor and actuator
            G{Itot_labeled(j,1),i+1} = [G{Itot_labeled(j,1),i+1}; G{Itot_labeled(j,2),i}*A0];
            w{Itot_labeled(j,1),i+1} = [w{Itot_labeled(j,1),i+1}; w{Itot_labeled(j,2),i}-maxH(G{Itot_labeled(j,2),i}, B1,Ay)-maxH(G{Itot_labeled(j,2),i}, B4,Au)];
        end
    end
    for j = 1:nodes
        S{j}(i+1) = Polyhedron(G{j,i+1},  w{j,i+1});
        S{j}(i+1).minHRep;
        if S{j}(i+1).volume <= epsilon
            fprintf('Node %d has zero volume\n', j)
            return
        end
        G{j,i+1} = S{j}(i+1).A;
        w{j,i+1} = S{j}(i+1).b;
    end
    for ii = 1:nodes
        if max(max(S{ii}(i+1).A*S{ii}(i).V'-(1+epsilon)*S{ii}(i+1).b ))<=1e-10
        check = check + 1;
        end
    end
    if check == nodes
        disp('convergence')
        convergence = 1;
    end
    i = i+1
    check = 0;
end
%
for iii = 1:nodes
    Smax{iii} = S{iii}(end);
end
Ssafe = Polyhedron();
for iiii = 1:nodes
    Ssafe = Polyhedron( [Ssafe.A; Smax{iiii}.A], [Ssafe.b; Smax{iiii}.b] );
end
% 
%
X0 = Z0;
Xsafe = Ssafe;
figure(4)
X0.plot('alpha', 0.2, 'color', 'blue')
hold on
Xsafe.plot('alpha', 1)
% legend('Constraints set', 'Maximal Safe Set under attack')
title('Poisoning Attacks')
hold off
Ssafe.volume
% metrics
I1_Lebesgue = (X0.volume - Ssafe.volume)/X0.volume
lambda_minkowski = Minkowski_fun(X0, Ssafe);
I2_Minkowski = 1 - lambda_minkowski

%%
% erosion
function hstar = maxH(Gk, E0, Hset)
A = Gk;
B = E0;
rows = length(A(:,1));
hstar = zeros(rows,1);
for i = 1:rows
    [x hstar(i)] = linprog(-(A(i,:)*B)', Hset.A, Hset.b);
end
hstar = -hstar;
end
%
function lambda = Minkowski_fun(Sn, Sa)
epsilon = 1e-5;
lambda_min = 0;
lambda_max = 1;
j = 1;
convergence = 0;
while convergence == 0
    lambda = (lambda_min + lambda_max)/2;
    Sn_new = lambda * Sn;
    vertices = Sn_new.V';
    number_of_vertices = size(vertices,2);
    count = 0;
    for i = 1:number_of_vertices
        if Sa.contains(vertices(:,i))
            count = count + 1;
        end
    end
    if count == number_of_vertices
        lambda_min = lambda;
    else
        lambda_max =lambda;
    end
    if abs(lambda_min - lambda_max) <= epsilon
        convergence = 1;
    end
    j = j + 1;
end
end

%%
function I = incid_mat(Nmin, Nmax)
%---------------------
I = []; 
% 'H' for Healthy
% 'A' for Attack
if Nmin ~= 0
    for i = 1:(Nmin + Nmax)
        if i == 1
            I = [I; 1 1 1];
        elseif i > 1 && i <= Nmin
            I = [I; i-1 i 1];
        elseif i > Nmin
            I = [I; i-1 i 2; i 1 1];
        end
    end
else
    for i = 1:(Nmax + 1)
        if i == 1
            I = [I; 1 1 1];
        else 
            I = [I; i-1 i 2; i 1 1];
        end
    end
end
end
%
%%
function  labels = make_labelsxu(I)
for i = 1:length(I(:,1))
    if I(i,3) == 1
        labels{i,1} = 'H'; % H for Healthy
    else
        labels{i,1} = 'A'; % A for Attack
    end
end
end
%
function labels_tot = make_labels_tot(I)
% H for Healthy
% A for Attack
for i = 1:length(I(:,1))
    if I(i,3) == 1 && I(i,4) == 1
        labels_tot{i,1} = 'HH';
    elseif I(i,3) == 1 && I(i,4) == 2
        labels_tot{i,1} = 'HA'; 
    elseif I(i,3) == 2 && I(i,4) == 1
        labels_tot{i,1} = 'AH';
    else
        labels_tot{i,1} = 'AA';
    end
end
end

%
