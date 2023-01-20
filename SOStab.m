classdef SOStab < handle
    %SOSTAB is a class to approximate ROA using SOS programming
    %   This class create instances that represents dynamical systems. The
    %   methods enables to calculate inner and outer approximations of the
    %   "constrained finite horizon region of attraction" of a target set K
    %   A_M^T(K) = {x_0 | s.t. \exists x, \dot{x} = f(x), x(0) = x_0, 
    %   x(t) \in M, x(T) \in K}
    
    properties
        dimension; % dimension of the problem
        phi_eq; % equilibrium state
        phi_str; % name of the variables (for plotting)
        delta_phi; % distance from the equilibirum to the max attainable value
        angle_eq=[]; % equilibrium angle (empty if no variables are angle)
        angle_str=[]; % name of the angles
        phi; % sdpvar for solving the problem
        t = sdpvar(1,1); % time variable
        polynomial_dynamic; % dynamic of the system : \dot(\phi) = p(\phi)
        D; % matrix of variable change : x = D(\phi - \phi_eq)
        invD; % inverse of the previous matrix
        d; % degree of polynomials, can change
        vcoef_inner; % coefficients of the solution v for the last calculated optimization
        wcoef_inner; % coefficients of the solution w for the last calculated optimization
        vcoef_outer; % coefficients of the solution v for the last calculated optimization
        wcoef_outer; % coefficients of the solution w for the last calculated optimization
        solution; % volume of the current ROA calculation
        A; % matrix s. t. the target set is defined by A(\phi-\phi_eq) < epsilon
        epsilon; % radius of the target set
        T; % Time horizon
        solver = 'mosek'; % can be changed manually
        verbose = 2; % can be changed manually
    end
    
    methods
        function obj = SOStab(phi_eq, delta_phi, phi_str)
            %SOSTAB Initialize an object for calculating ROA approximations
            %
            obj.dimension = size(phi_eq,1);
            obj.phi_eq = phi_eq;
            obj.delta_phi = reshape(delta_phi, obj.dimension, 1);
            obj.phi_str = reshape(phi_str, obj.dimension, 1);
            obj.phi = sdpvar(obj.dimension,1);
            invD = zeros(obj.dimension,1);
            for k=1:obj.dimension
                invD(k) = obj.delta_phi(k); % max([(phi_max(k) - phi_eq(k)), (phi_eq(k) - phi_min(k))]);% 
            end
            obj.invD = diag(invD);
            obj.D = obj.invD ^(-1);

        end
        
        function [sol, vc, wc] = solveoptim_outer(obj, A, epsilon, d, T)
            %SOLVEOPTIM_INNER solves the outer approximation problem
            
            n = obj.dimension;
            s = sdpvar(1,1); % s= t / T;
            x = sdpvar(n,1);
            obj.d = d;
            obj.A = A;
            obj.T = T;
            obj.epsilon = epsilon;
        
            % definition of K = hypersphere of radius eps
            g = epsilon^2 - x' * ((A*obj.invD)' * A*obj.invD) * x;
            
            % definition of X = hypercube [-1,1]^n
            h = sdpvar(n,1);
            for k=1:n
                h(k) = 1 - x(k)^2;
            end
            [w, wc] = polynomial(x, d);
            [v, vc] = polynomial([s; x], d);
            
            zz = monolist(x,d);
            % Integration of zz on X, optimized way:
            y = ones(size(zz));
            fz = sdisplay(zz);
            for k=1:size(zz,1)
                f = fz{k,1};
                iszero = false;
                if contains(f, ')*') || f(end)==')' % search for a monomial of degree 1
                    y(k) =0;
                    continue
                end
                for j=3:2:d
                    if contains(f, ['^' num2str(j)]) % search for odd degree monomials
                        y(k) =0;
                        iszero = true;
                        continue
                    end
                end
                if iszero
                    continue
                end
                for j=1:size(x,1) % integration for other situations
                    iscons = true;
                    for deg = 2:2:d
                       if contains(f, [num2str(j) ')^' num2str(deg)]) % integration of even degree monomials
                           y(k)= y(k)*2/(deg+1);
                           iscons = false;
                           continue
                       end
                    end
                    if iscons % integration of constant over [-1,1]
                        y(k)= y(k)*2;
                    end
                end
                        
            end
            objectif = wc'* y;
            var = [wc; vc];
            % Constraints
            % v(1 ,.) >= 0 on K
            [q1, qc1] = polynomial(x, d -2);
            con = [sos(q1), sos(replace(v, s, 1) - g * q1)];
            var = [var; qc1];
            
            % w >= v(0 ,.) + 1 on X
            hs = 0;
            for i=1:n
                [q, qc] = polynomial(x, d-2);
                con = [con, sos(q)];
                hs = hs + h(i) * q;
                var = [var; qc];
            end
            con = [con, sos(w-replace(v, s , 0) - 1 - hs)];
            
            % dot(v) <= 0 on X x [0,T]
            l = (1 - s) * s; % definition of [0,1]
            hs = 0;
            for i=1:n
                [q, qc] = polynomial([s; x], d-2);
                con = [con, sos(q)];
                hs = hs + h(i) * q;
                var = [var; qc];
            end
            [q2, qc2] = polynomial([s; x], d-2);
            con = [con, sos(q2)];
            f = (obj.D) * replace(obj.polynomial_dynamic, [obj.t;obj.phi], [s*T;(obj.invD)*x+obj.phi_eq]); % f defined thanks to variable change
        
            con = [con, sos(-(max(1, 1/T) * jacobian(v,s) + max(1,T) * jacobian(v,x) * f) - hs - l * q2)];
            var = [var; qc2];
            
            % w >= 0 on X
            hs = 0;
            for i=1:n
                [q, qc] = polynomial(x,d-2);
                con = [con, sos(q)];
                hs = hs + h(i) * q;
                var = [var; qc];
            end
            con = [con, sos(w-hs)];
            
            % Optimize
            ops = sdpsettings('solver', obj.solver, 'verbose', obj.verbose);
            solvesos(con, objectif, ops, var);
            sol = double(wc' * y);
            obj.wcoef_outer = wc;
            obj.vcoef_outer = vc;
        end

        function [sol, vc, wc] = solveoptim_inner(obj, A, epsilon, d, T)
            %SOLVEOPTIM_INNER Solves the inner approximation problem
            
            n = obj.dimension;
            s = sdpvar(1,1); % s= t / T;
            x = sdpvar(n,1);
            obj.d = d;
            obj.A = A;
            obj.T = T;
            obj.epsilon = epsilon;
        
            % definition of K = hypersphere of radius eps
            g = epsilon^2 - x' * ((A*obj.invD)' * A*obj.invD) * x;
            
            % definition of X = hypercube [-1,1]^n
            h = sdpvar(n,1);
            for k=1:n
                h(k) = 1 - x(k)^2;
            end
            [w, wc] = polynomial(x, d);
            [v, vc] = polynomial([s; x], d);
            
            zz = monolist(x,d);
            % Integration of zz on X, optimized way:
            y = ones(size(zz));
            fz = sdisplay(zz);
            for k=1:size(zz,1)
                f = fz{k,1};
                iszero = false;
                if contains(f, ')*') || f(end)==')' % search for a monomial of degree 1
                    y(k) =0;
                    continue
                end
                for j=3:2:d
                    if contains(f, ['^' num2str(j)]) % search for odd degree monomials
                        y(k) =0;
                        iszero = true;
                        continue
                    end
                end
                if iszero
                    continue
                end
                for j=1:size(x,1) % integration for other situations
                    iscons = true;
                    for deg = 2:2:d
                       if contains(f, [num2str(j) ')^' num2str(deg)]) % integration of even degree monomials
                           y(k)= y(k)*2/(deg+1);
                           iscons = false;
                           continue
                       end
                    end
                    if iscons % integration of constant over [-1,1]
                        y(k)= y(k)*2;
                    end
                end
                        
            end
            objectif = wc'* y;
            var = [wc; vc];
            con = [];
            % Constraints
            % v(1,.) >= 0 on X\K = {g(x) <=0, h(x)>=0}
            [q1, qc1] = polynomial(x, d -2);
            var = [var; qc1];
            con = [con, sos(q1)];
            hs = 0;
            for i=1:n
                [q, qc] = polynomial(x, d-2);
                con = [con, sos(q)];
                hs = hs + h(i) * q;
                var = [var; qc];
            end
            con = [con, sos(replace(v, s, 1) + g * q1 - hs)];

            % v(1,.) <= 0 on K = {g(x) >=0}
            %[q1, qc1] = polynomial(x, d-2);
            %var = [var; qc1];
            %con = [con, sos(q1)];
            %con = [con, sos(-replace(v, s, 1) - g * q1)];
            
            % w >= v(0 ,.) + 1 on X
            hs = 0;
            for i=1:n
                [q, qc] = polynomial(x, d-2);
                con = [con, sos(q)];
                hs = hs + h(i) * q;
                var = [var; qc];
            end
            con = [con, sos(w-replace(v, s , 0) - 1 - hs)];
            
            % dot(v) <= 0 on X x [0,T]
            l = (1 - s) * s; % definition of [0,1]
            hs = 0;
            for i=1:n
                [q, qc] = polynomial([s; x], d-2);
                con = [con, sos(q)];
                hs = hs + h(i) * q;
                var = [var; qc];
            end
            [q2, qc2] = polynomial([s; x], d-2);
            con = [con, sos(q2)];
            f = (obj.D) * replace(obj.polynomial_dynamic, [obj.t;obj.phi], [s*T;(obj.invD)*x+obj.phi_eq]); % f defined thanks to variable change
        
            con = [con, sos(-(min(1,1/T) * jacobian(v,s) + min(1,T) * jacobian(v,x) * f) - hs - l * q2)];
            var = [var; qc2];
            
            % w >= 0 on X
            hs = 0;
            for i=1:n
                [q, qc] = polynomial(x,d-2);
                con = [con, sos(q)];
                hs = hs + h(i) * q;
                var = [var; qc];
            end
            con = [con, sos(w-hs)];

            % v(t,x) >= 0 on X_frontier x [0,T]
            hs = 0;
            prod_hi = 1;
            for i=1:n
                [q, qc] = polynomial([s; x], d-2);
                con = [con, sos(q)];
                hs = hs + h(i) * q;
                prod_hi = h(i) * prod_hi;
                var = [var; qc];
            end
            [p, pc] = polynomial([s; x], d-2);
            [q2, qc2] = polynomial([s; x], d-2);
            con = [con, sos(q2), sos(v - hs - l * q2 - p * prod_hi)];
            var = [var; qc2; pc];
            
            % Optimize
            ops = sdpsettings('solver', obj.solver, 'verbose', obj.verbose);
            solvesos(con, objectif, ops, var);
            sol = double(wc' * y);
            obj.wcoef_inner = wc;
            obj.vcoef_inner = vc;
        end
        
        function [] = plot_roa(obj, ind1, ind2, approximation, mesh_size)
            % PLOT_ROA plots the graph of the ROA approximation
            % as a function of the two variables indicated by the indices
            % ind1 and ind2. The other variables are taken to be at
            % equilibrium.

                        %%%%%
            % INPUTS
            %   ind1, ind2: indices of the variables to be plotted.
            %       If an angle is plotted, the variable should be an array
            %       of three elements. The first two indicate the indices
            %       of respectively sin(angle) and (1-cos(angle))/2, the
            %       last one indicates the indice of the plotted angle in
            %       the angle array
            %   approximation: 'i' or 'inner' for inner approximation of
            %      the ROA, 'o' or 'outer' for the outer approximation
            %   mesh_size (optionnal): size of the mesh for plotting,
            %   2-values array
            %%%%%

            if nargin < 5
                mesh_size = [40,40];
            end
            if all(approximation == 'i') || all(approximation == 'inner')
                curvcoef = obj.vcoef_inner;
                col = '-b';
                app_str = 'inner';
                border_value = 1.2;
            else
                curvcoef = obj.vcoef_outer;
                col = '-r';
                app_str = 'outer';
                border_value = 0.8;
            end
            
            angle1 =  any(size(ind1) > 1);
            angle2 =  any(size(ind2) > 1);
            if angle1
                if any(size(ind1) > 2)
                    ind1 = reshape(ind1,1,3);
                else
                    ind1 = [reshape(ind1,1,2) 1];
                end
                vec_x = linspace(obj.angle_eq(ind1(3))-pi, obj.angle_eq(ind1(3))+pi, mesh_size(1));
            else
                vec_x = linspace(obj.phi_eq(ind1)-obj.delta_phi(ind1), obj.phi_eq(ind1)+obj.delta_phi(ind1),mesh_size(1));
            end
            if angle2
                if any(size(ind2) > 2)
                    ind2 = reshape(ind2,1,3);
                else
                    ind2 = [reshape(ind2,1,2) 1];
                end
                vec_y = linspace(obj.angle_eq(ind2(3))-pi, obj.angle_eq(ind2(3))+pi, mesh_size(2));
            else
                vec_y = linspace(obj.phi_eq(ind2)-obj.delta_phi(ind2),obj.phi_eq(ind2)+obj.delta_phi(ind2), mesh_size(2));
            end
            [plot_x,plot_y] = meshgrid(vec_x, vec_y);
            
            X = sdpvar(obj.dimension,1);
            vv = monolist([obj.t;X], obj.d);
            k = replace(curvcoef' * vv, obj.t, 0) + 1;
            k = replace(k, curvcoef, double(curvcoef));
            for j=1:obj.dimension
                if not(any([ind1, ind2] == j))
                    k = replace(k, X(j), 0);
                end
            end
            k = vectorize(sdisplay(k));
            Z = zeros(mesh_size);
            for i=1:mesh_size(1)
                for j = 1:mesh_size(2)
                    if angle1
                        [X(ind1(1)) ,X(ind1(2))] = deal(sin(plot_x(i,j) - obj.angle_eq(ind1(3))) ,...
                        (1 - cos(plot_x(i,j) - obj.angle_eq(ind1(3)))) / 2);
                    else
                        X(ind1) = (plot_x(i,j) - obj.phi_eq(ind1)) * obj.D(ind1,ind1);
                    end
                    if angle2
                        [X(ind2(1)) ,X(ind2(2))] = deal(sin(plot_y(i,j) - obj.angle_eq(ind2(3))) ,...
                        (1 - cos(plot_y(i,j) - obj.angle_eq(ind2(3)))) / 2);
                    else
                        X(ind2) = (plot_y(i,j) - obj.phi_eq(ind2)) * obj.D(ind2,ind2);
                    end
                    Z(i,j) = eval(k);
                end
            end
            Z(1,:) = border_value*ones(1,mesh_size(2));
            Z(mesh_size(1),:) = border_value*ones(1,mesh_size(2));
            Z(:,1) = border_value*ones(mesh_size(1),1);
            Z(:,mesh_size(2)) = border_value*ones(mesh_size(1),1);
            figure(1)
            hold on
            contour(plot_x, plot_y, Z, [1 1], col,'linewidth', 2, 'DisplayName', ['Region of attraction at T : ' app_str]);
            if angle1
                coef1 = 1;
                xeq = obj.angle_eq(ind1(3));
            else
                coef1 = obj.A(ind1,ind1);
                xeq = obj.phi_eq(ind1);
            end
            if angle2
                coef2 = 1;
                yeq = obj.angle_eq(ind2(3));
            else
                coef2 = obj.A(ind2,ind2);
                yeq = obj.phi_eq(ind2);
            end
            parametric_angle = linspace(0,2*pi,100);
            if all(app_str == 'inner')
                plot(xeq + obj.epsilon/coef1 * cos(parametric_angle), yeq + obj.epsilon/coef2 * sin(parametric_angle), '-k','linewidth', 2, 'DisplayName','Target set');
            end
            %ZZ = (coef1*plot_x).^2+(coef2*plot_y).^2;
            %contour(plot_x, plot_y, ZZ, [obj.epsilon obj.epsilon], '-k','linewidth', 2);
            % rectangle('position',[-1+ rho_eq -pi+theta_eq 2 2* pi]);
            hold off
            %axis equal
            title('Region of attraction over X = [- 1, 1]')
            legend('Location','southwest')
            
            if angle1
                x_lab = obj.angle_str(ind1(3));
            else
                x_lab = obj.phi_str(ind1);
            end
            if angle2
                y_lab = obj.angle_str(ind2(3));
            else
                y_lab = obj.phi_str(ind2);
            end
            title(strjoin(['ROA(',x_lab,',',y_lab,')']))
            xlabel(x_lab);
            ylabel(y_lab);
        end
    
        function [] = plot_w(obj, ind1, ind2, approximation, mesh_size)
            % PLOT_W plots the graph of w
            % as a function of the two variables indicated by the indices
            % ind1 and ind2. The other variables are taken to be at
            % equilibrium.

            %%%%%
            % INPUTS
            %   ind1, ind2: indices of the variables to be plotted. All other
            %       variables will be taken at equilibrium
            %       If an angle is plotted, the variable should be an array
            %       of three elements. The first two indicate the indices
            %       of respectively sin(angle) and (1-cos(angle))/2, the
            %       last one indicates the indice of the angle in the angle
            %       array
            %   approximation: 'i' or 'inner' for inner approximation of
            %      the ROA, 'o' or 'outer' for the outer approximation
            %   mesh_size (optionnal): size of the mesh for plotting,
            %   2-values array
            %%%%%
            if nargin < 5
                mesh_size = [30,30];
            end

            if all(approximation == 'i') || all(approximation == 'inner')
                curwcoef = obj.wcoef_inner;
            else
                curwcoef = obj.wcoef_outer;
            end
            
            angle1 =  any(size(ind1) > 1);
            angle2 =  any(size(ind2) > 1);
            if angle1
                if any(size(ind1) > 2)
                    ind1 = reshape(ind1,1,3);
                else
                    ind1 = [reshape(ind1,1,2) 1];
                end
                vec_x = linspace(obj.angle_eq(ind1(3))-pi, obj.angle_eq(ind1(3))+pi, mesh_size(1));
            else
                vec_x = linspace(obj.phi_eq(ind1)-obj.delta_phi(ind1), obj.phi_eq(ind1)+obj.delta_phi(ind1),mesh_size(1));
            end
            if angle2
                if any(size(ind2) > 2)
                    ind2 = reshape(ind2,1,3);
                else
                    ind2 = [reshape(ind2,1,2) 1];
                end
                vec_y = linspace(obj.angle_eq(ind2(3))-pi, obj.angle_eq(ind2(3))+pi, mesh_size(2));
            else
                vec_y = linspace(obj.phi_eq(ind2)-obj.delta_phi(ind2),obj.phi_eq(ind2)+obj.delta_phi(ind2), mesh_size(2));
            end
            [plot_x,plot_y] = meshgrid(vec_x, vec_y);
            
            X = sdpvar(obj.dimension,1);
            ww = monolist(X, obj.d);
            k = replace(curwcoef'*ww, curwcoef, double(curwcoef));
            for j=1:obj.dimension
                if not(any([ind1, ind2] == j))
                    k = replace(k, X(j), 0);
                end
            end
            k = vectorize(sdisplay(k));
            K = zeros(mesh_size); L = zeros(mesh_size);
            for i=1:mesh_size(1)
                for j = 1:mesh_size(2)
                    if angle1
                        [X(ind1(1)) ,X(ind1(2))] = deal(sin(plot_x(i,j) - obj.angle_eq(ind1(3))) ,...
                        (1 - cos(plot_x(i,j) - obj.angle_eq(ind1(3)))) / 2);
                    else
                        X(ind1) = (plot_x(i,j) - obj.phi_eq(ind1)) * obj.D(ind1,ind1);
                    end
                    if angle2
                        [X(ind2(1)) ,X(ind2(2))] = deal(sin(plot_y(i,j) - obj.angle_eq(ind2(3))) ,...
                        (1 - cos(plot_y(i,j) - obj.angle_eq(ind2(3)))) / 2);
                    else
                        X(ind2) = (plot_y(i,j) - obj.phi_eq(ind2)) * obj.D(ind2,ind2);
                    end
                    K(i,j) = eval(k);
                    L(i,j) = 1;
                end
            end
            figure(2)
            hold on
            mesh(plot_x, plot_y, double(K));
            plot3(plot_x, plot_y, double(L) ,'k');
            view([37.5 30]);
            hold off
            if angle1
                x_lab = obj.angle_str(ind1(3));
            else
                x_lab = obj.phi_str(ind1);
            end
            if angle2
                y_lab = obj.angle_str(ind2(3));
            else
                y_lab = obj.phi_str(ind2);
            end
            title(strjoin(['w(',x_lab,',',y_lab,')']))
            xlabel(x_lab);
            ylabel(y_lab);
        end

        function [] = plot_v(obj, ind1, ind2, tau, approximation, mesh_size)
            % PLOT_V plots the graph of v 
            % as a function of the time tau and the two variables indicated
            % by the indices ind1 and ind2. The other variables are taken 
            % at equilibrium.

            %%%%%
            % INPUTS
            %   ind1, ind2: indices of the variables to be plotted. All other
            %       variables will be taken at equilibrium
            %       If an angle is plotted, the variable should be an array
            %       of three elements. The first two indicate the indices
            %       of respectively sin(angle) and (1-cos(angle))/2, the
            %       last one indicates the indice of the angle in the angle
            %       array
            %   approximation: 'i' or 'inner' for inner approximation of
            %      the ROA, 'o' or 'outer' for the outer approximation
            %   mesh_size (optionnal): size of the mesh for plotting,
            %   2-values array
            %%%%%
            if nargin < 5
                mesh_size = [30,30];
            end

            if all(approximation == 'i') || all(approximation == 'inner')
                curwcoef = obj.wcoef_inner;
            else
                curwcoef = obj.wcoef_outer;
            end
            
            angle1 =  any(size(ind1) > 1);
            angle2 =  any(size(ind2) > 1);
            if angle1
                if any(size(ind1) > 2)
                    ind1 = reshape(ind1,1,3);
                else
                    ind1 = [reshape(ind1,1,2) 1];
                end
                vec_x = linspace(obj.angle_eq(ind1(3))-pi, obj.angle_eq(ind1(3))+pi, mesh_size(1));
            else
                vec_x = linspace(obj.phi_eq(ind1)-obj.delta_phi(ind1), obj.phi_eq(ind1)+obj.delta_phi(ind1),mesh_size(1));
            end
            if angle2
                if any(size(ind2) > 2)
                    ind2 = reshape(ind2,1,3);
                else
                    ind2 = [reshape(ind2,1,2) 1];
                end
                vec_y = linspace(obj.angle_eq(ind2(3))-pi, obj.angle_eq(ind2(3))+pi, mesh_size(2));
            else
                vec_y = linspace(obj.phi_eq(ind2)-obj.delta_phi(ind2),obj.phi_eq(ind2)+obj.delta_phi(ind2), mesh_size(2));
            end
            [plot_x,plot_y] = meshgrid(vec_x, vec_y);
            
            X = sdpvar(obj.dimension,1);
            vv = monolist([obj.t;X], obj.d);
            k = replace(curvcoef' * vv, obj.t, tau/obj.T);
            k = replace(k, curvcoef, double(curvcoef));
            for j=1:obj.dimension
                if not(any([ind1, ind2] == j))
                    k = replace(k, X(j), 0);
                end
            end
            k = vectorize(sdisplay(k));
            K = zeros(mesh_size); L = zeros(mesh_size);
            for i=1:mesh_size(1)
                for j = 1:mesh_size(2)
                    if angle1
                        [X(ind1(1)) ,X(ind1(2))] = deal(sin(plot_x(i,j) - obj.angle_eq(ind1(3))) ,...
                        (1 - cos(plot_x(i,j) - obj.angle_eq(ind1(3)))) / 2);
                    else
                        X(ind1) = (plot_x(i,j) - obj.phi_eq(ind1)) * obj.D(ind1,ind1);
                    end
                    if angle2
                        [X(ind2(1)) ,X(ind2(2))] = deal(sin(plot_y(i,j) - obj.angle_eq(ind2(3))) ,...
                        (1 - cos(plot_y(i,j) - obj.angle_eq(ind2(3)))) / 2);
                    else
                        X(ind2) = (plot_y(i,j) - obj.phi_eq(ind2)) * obj.D(ind2,ind2);
                    end
                    K(i,j) = eval(k);
                    L(i,j) = 1;
                end
            end
            figure(2)
            hold on
            mesh(plot_x, plot_y, double(K));
            plot3(plot_x, plot_y, double(L) ,'k');
            view([37.5 30]);
            hold off
            if angle1
                x_lab = obj.angle_str(ind1(3));
            else
                x_lab = obj.phi_str(ind1);
            end
            if angle2
                y_lab = obj.angle_str(ind2(3));
            else
                y_lab = obj.phi_str(ind2);
            end
            title(strjoin(['v(',x_lab,',',y_lab,',',tau,')']))
            xlabel(x_lab);
            ylabel(y_lab);
        end
    end
end
