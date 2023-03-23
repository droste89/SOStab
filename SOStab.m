classdef SOStab < handle
    %SOSTAB is a class to approximate ROA using SOS programming
    %   This class create instances that represents dynamical systems. The
    %   methods enables to calculate inner and outer approximations of the
    %   "constrained finite horizon region of attraction" of a target set K
    %   A_M^T(K) = {x_0 | s.t. \exists x, \dot{x} = f(x), x(0) = x_0, 
    %   x(t) \in M, x(T) \in K}
    
    properties
        dimension; % dimension of the problem
        x_eq; % equilibrium state
        delta_x; % distance from the equilibirum to the max attainable value
        angle_eq=[]; % equilibrium angle (empty if no variables are angle)
        angle_ind; % indices of each angle
        x; % sdpvar for solving the problem
        t = sdpvar(1,1); % time variable
        dynamics; % dynamic of the system : \dot(x) = p(x)
        D; % matrix of variable change : x = D(x - x_eq)
        invD; % inverse of the previous matrix
        d; % degree of polynomials, can change
        vcoef_inner; % coefficients of the solution v for the last calculated optimization
        wcoef_inner; % coefficients of the solution w for the last calculated optimization
        vcoef_outer; % coefficients of the solution v for the last calculated optimization
        wcoef_outer; % coefficients of the solution w for the last calculated optimization
        solution; % volume of the current ROA calculation
        A; % matrix s. t. the target set is defined by A(x-x_eq) < epsilon
        epsilon; % radius of the target set
        T; % Time horizon
        solver = 'mosek'; % can be changed manually
        verbose = 2; % can be changed manually
    end
    
    methods
        function obj = SOStab(x_eq, delta_x, parse_matrix)
            %SOSTAB Initialize an object for calculating ROA approximations
            %
            %
            %%%%%
            % INPUTS
            %   x_eq: equilibrium point of the system, must have the size of the 
            %   dimension of the problem
            %   delta_x: range around the equilibrium defining the admissible set,
            %   must have the same size as x_eq
            %   parse_matrix (optionnal): if the problem involves angles, this 
            %   matrix should contain rows defining each angles and relating them
            %   to the corresponding variables in the first two vectors:
            %   each row should be of the type [theta_eq, ind1, ind2] where ind1
            %   (resp. ind2) indicates the place of the variables corresponding
            %   to sin(theta) (resp. cos(theta)).
            %%%%%

            if nargin < 3
                parse_matrix = [];
            end
            obj.dimension = max(size(x_eq));
            obj.x_eq = reshape(x_eq, obj.dimension, 1);
            obj.delta_x = reshape(delta_x, obj.dimension, 1);
            obj.x = sdpvar(obj.dimension,1);
            invD = zeros(obj.dimension,1);
            for j=1:size(parse_matrix,1) % equilibrium, sine, cosine
                obj.angle_eq(j, 1) = parse_matrix(j,1);
                obj.angle_ind(j, :) = parse_matrix(j, 2:3);
                obj.x_eq(parse_matrix(j, 2)) = 0;
                obj.x_eq(parse_matrix(j, 3)) = 0;
                if any(obj.delta_x(parse_matrix(j, 2:3))>1)
                    error("Sin/Cosine range can't be superior to 1")
                end
                if parse_matrix(j,2) == parse_matrix(j,3)
                    error("Sin and cosine can't have the same indice")
                end
            end
            for k=1:obj.dimension
                invD(k) = obj.delta_x(k);
            end
            obj.invD = diag(invD);
            obj.D = obj.invD ^(-1);

        end
        
        function [sol, vc, wc] = SoS_out(obj, d, T, epsilon, A)
            %SOS_OUT solves the outer approximation problem

            n = obj.dimension;
            if isempty(obj.dynamics)
                error("You forgot to input the dynamics")
            end
            if nargin < 5
                A = eye(n)
            end
            s = sdpvar(1,1); % s= t / T;
            z = sdpvar(n,1);
            obj.d = d;
            obj.A = A;
            obj.T = T;
            obj.epsilon = epsilon;
        
            % definition of K = hypersphere of radius eps
            g = epsilon^2 - z' * ((A*obj.invD)' * A*obj.invD) * z;
            
            % definition of X = hypercube [-1,1]^n
            h = sdpvar(n,1);
            for k=1:n
                h(k) = 1 - z(k)^2;
            end
            [w, wc] = polynomial(z, d);
            [v, vc] = polynomial([s; z], d);
            
            zz = monolist(z,d);
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
                for j=1:size(z,1) % integration for other situations
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
            [q1, qc1] = polynomial(z, d -2);
            con = [sos(q1), sos(replace(v, s, 1) - g * q1)];
            var = [var; qc1];
            
            % w >= v(0 ,.) + 1 on X
            hs = 0;
            for i=1:n
                [q, qc] = polynomial(z, d-2);
                con = [con, sos(q)];
                hs = hs + h(i) * q;
                var = [var; qc];
            end
            con = [con, sos(w-replace(v, s , 0) - 1 - hs)];
            
            % dot(v) <= 0 on X x [0,T]
            l = (1 - s) * s; % definition of [0,1]
            hs = 0;
            for i=1:n
                [q, qc] = polynomial([s; z], d-2);
                con = [con, sos(q)];
                hs = hs + h(i) * q;
                var = [var; qc];
            end
            [q2, qc2] = polynomial([s; z], d-2);
            con = [con, sos(q2)];
			
			if isempty(obj.angle_eq)
				f = (obj.D) * replace(obj.dynamics, [obj.t;obj.x], [s*T;(obj.invD)*z+obj.x_eq]); % f defined thanks to variable change
			else
				f = replace(obj.dynamics, obj.t, s*T);
				for j=1:size(obj.angle_eq, 1)
					% obj.x(obj.angle_ind(j,1)) is sin(theta), obj.x(obj.angle_ind(j,2)) is cos(theta)
					% z(obj.angle_ind(j,1)) is sin(theta - theta_eq), z(obj.angle_ind(j,2)) is 0.5-0.5cos(theta - theta_eq)
					% and sin(theta) = sin(theta-theta_eq)cos(theta_eq) + sin(theta_eq)(1-2(0.5-0.5cos(theta-theta_eq)))
					% similar for cos(theta)
                    costhetaeq = cos(obj.angle_eq(j));
                    sinthetaeq = sin(obj.angle_eq(j));
                    
					f = replace(f, [obj.x(obj.angle_ind(j,1)); obj.x(obj.angle_ind(j,2))], ...
						[ sinthetaeq * (1 - 2 * z(obj.angle_ind(j,2))) + costhetaeq * z(obj.angle_ind(j,1)); ...
						costhetaeq * (1 - 2 * z(obj.angle_ind(j,2))) - sinthetaeq * z(obj.angle_ind(j,1))]);
                    % sin(theta-theta_eq) = sin(theta)cos(theta_eq)-cos(theta)sin(theta_eq)
                    % d(sin(theta-theta_eq)) = d(sin(theta))cos(theta_eq)-d(cos(theta))sin(theta_eq)
                    % 
                    % -0.5d(cos(theta-theta_eq)) = -0.5(cos(theta)cos(theta_eq) + sin(theta)sin(theta_eq)))
                    [f(obj.angle_ind(j,1)), f(obj.angle_ind(j,2))] = deal(costhetaeq*f(obj.angle_ind(j,1))-sinthetaeq*f(obj.angle_ind(j,2)),...
                    -0.5*(costhetaeq*f(obj.angle_ind(j,2))+sinthetaeq*f(obj.angle_ind(j,1))));
				end
				f = (obj.D) * replace(f, [obj.x], [(obj.invD)*z+obj.x_eq]); % f defined thanks to variable change, note that variables corresponding to angles are already replaced
            end

            con = [con, sos(-(max(1, 1/T) * jacobian(v,s) + max(1,T) * jacobian(v,z) * f) - hs - l * q2)];
            var = [var; qc2];
            
            % w >= 0 on X
            hs = 0;
            for i=1:n
                [q, qc] = polynomial(z,d-2);
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

        function [sol, vc, wc] = SoS_in(obj, d, T, epsilon, A)
            %SOS_IN Solves the inner approximation problem
            n = obj.dimension;
            if isempty(obj.dynamics)
                error("You forgot to input the dynamics")
            end
            if nargin < 5
                A = eye(n)
            end
            
            n = obj.dimension;
            s = sdpvar(1,1); % s= t / T;
            z = sdpvar(n,1);
            obj.d = d;
            obj.A = A;
            obj.T = T;
            obj.epsilon = epsilon;
        
            % definition of K = hypersphere of radius eps
            g = epsilon^2 - z' * ((A*obj.invD)' * A*obj.invD) * z;
            
            % definition of X = hypercube [-1,1]^n
            h = sdpvar(n,1);
            for k=1:n
                h(k) = 1 - z(k)^2;
            end
            [w, wc] = polynomial(z, d);
            [v, vc] = polynomial([s; z], d);
            
            zz = monolist(z,d);
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
                for j=1:size(z,1) % integration for other situations
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
            % v(1,.) >= 0 on X\K = {g(z) <=0, h(z)>=0}
            [q1, qc1] = polynomial(z, d -2);
            var = [var; qc1];
            con = [con, sos(q1)];
            hs = 0;
            for i=1:n
                [q, qc] = polynomial(z, d-2);
                con = [con, sos(q)];
                hs = hs + h(i) * q;
                var = [var; qc];
            end
            con = [con, sos(replace(v, s, 1) + g * q1 - hs)];

            % v(1,.) <= 0 on K = {g(z) >=0}
            %[q1, qc1] = polynomial(z, d-2);
            %var = [var; qc1];
            %con = [con, sos(q1)];
            %con = [con, sos(-replace(v, s, 1) - g * q1)];
            
            % w >= v(0 ,.) + 1 on X
            hs = 0;
            for i=1:n
                [q, qc] = polynomial(z, d-2);
                con = [con, sos(q)];
                hs = hs + h(i) * q;
                var = [var; qc];
            end
            con = [con, sos(w-replace(v, s , 0) - 1 - hs)];
            
            % dot(v) <= 0 on X x [0,T]
            l = (1 - s) * s; % definition of [0,1]
            hs = 0;
            for i=1:n
                [q, qc] = polynomial([s; z], d-2);
                con = [con, sos(q)];
                hs = hs + h(i) * q;
                var = [var; qc];
            end
            [q2, qc2] = polynomial([s; z], d-2);
            con = [con, sos(q2)];
			
			if isempty(obj.angle_eq)
				f = (obj.D) * replace(obj.dynamics, [obj.t;obj.x], [s*T;(obj.invD)*z+obj.x_eq]); % f defined thanks to variable change
			else
				f = replace(obj.dynamics, obj.t, s*T);
				for j=1:size(obj.angle_eq, 1)
					% obj.x(obj.angle_ind(j,1)) is sin(theta), obj.x(obj.angle_ind(j,2)) is cos(theta)
					% z(obj.angle_ind(j,1)) is sin(theta - theta_eq), z(obj.angle_ind(j,2)) is 0.5-0.5cos(theta - theta_eq)
					% and sin(theta) = sin(theta-theta_eq)cos(theta_eq) + sin(theta_eq)(1-2(0.5-0.5cos(theta-theta_eq)))
					% similar for cos(theta)
                    costhetaeq = cos(obj.angle_eq(j));
                    sinthetaeq = sin(obj.angle_eq(j));
                    
					f = replace(f, [obj.x(obj.angle_ind(j,1)); obj.x(obj.angle_ind(j,2))], ...
						[ sinthetaeq * (1 - 2 * z(obj.angle_ind(j,2))) + costhetaeq * z(obj.angle_ind(j,1)); ...
						costhetaeq * (1 - 2 * z(obj.angle_ind(j,2))) - sinthetaeq * z(obj.angle_ind(j,1))]);
                    % sin(theta-theta_eq) = sin(theta)cos(theta_eq)-cos(theta)sin(theta_eq)
                    % d(sin(theta-theta_eq)) = d(sin(theta))cos(theta_eq)-d(cos(theta))sin(theta_eq)
                    % 
                    % -0.5d(cos(theta-theta_eq)) = -0.5(cos(theta)cos(theta_eq) + sin(theta)sin(theta_eq)))
                    [f(obj.angle_ind(j,1)), f(obj.angle_ind(j,2))] = deal(costhetaeq*f(obj.angle_ind(j,1))-sinthetaeq*f(obj.angle_ind(j,2)),...
                    -0.5*(costhetaeq*f(obj.angle_ind(j,2))+sinthetaeq*f(obj.angle_ind(j,1))));
				end
				f = (obj.D) * replace(f, [obj.x], [(obj.invD)*z+obj.x_eq]); % f defined thanks to variable change, note that variables corresponding to angles are already replaced
            end
        
            con = [con, sos(-(min(1,1/T) * jacobian(v,s) + min(1,T) * jacobian(v,z) * f) - hs - l * q2)];
            var = [var; qc2];
            
            % w >= 0 on X
            hs = 0;
            for i=1:n
                [q, qc] = polynomial(z,d-2);
                con = [con, sos(q)];
                hs = hs + h(i) * q;
                var = [var; qc];
            end
            con = [con, sos(w-hs)];

            % v(t,z) >= 0 on X_frontier x [0,T]
            hs = 0;
            prod_hi = 1;
            for i=1:n
                [q, qc] = polynomial([s; z], d-2);
                con = [con, sos(q)];
                hs = hs + h(i) * q;
                prod_hi = h(i) * prod_hi;
                var = [var; qc];
            end
            [p, pc] = polynomial([s; z], d-2);
            [q2, qc2] = polynomial([s; z], d-2);
            con = [con, sos(q2), sos(v - hs - l * q2 - p * prod_hi)];
            var = [var; qc2; pc];
            
            % Optimize
            ops = sdpsettings('solver', obj.solver, 'verbose', obj.verbose);
            solvesos(con, objectif, ops, var);
            sol = double(wc' * y);
            obj.wcoef_inner = wc;
            obj.vcoef_inner = vc;
        end
        
        function [] = plot_roa(obj, ind1, ind2, approximation, plot_target, str1, str2, mesh_size)
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
			%	str1, str2: names to be displayed in the plot for respectively
			%	ind1 and ind2 variables.
			%	plot_target: 0 or 1, to indicate whether to plot the target set
			%	or not
            %   mesh_size (optionnal): size of the mesh for plotting,
            %   2-values array
            %%%%%

            if nargin < 8
                mesh_size = [40,40];
			end
			if nargin < 5
				plot_target = 0;
			end
			if nargin < 7
                str1 = ['x_{' num2str(ind1) '}'];
                str2 = ['x_{' num2str(ind2) '}'];
            end
            if all(approximation == 'i') || all(approximation == 'inner')
                curvcoef = obj.vcoef_inner;
                col = '-b';
                app_str = 'Inner';
                border_value = 1.2;
            else
                curvcoef = obj.vcoef_outer;
                col = '-r';
                app_str = 'Outer';
                border_value = 0.8;
            end
            
            angle1 = any(size(ind1) > 1);
            angle2 = any(size(ind2) > 1);
            if angle1
                if any(size(ind1) > 2)
                    ind1 = reshape(ind1,1,3);
                else
                    ind1 = [reshape(ind1,1,2) 1];
                end
                vec_x = linspace(obj.angle_eq(ind1(3))-pi, obj.angle_eq(ind1(3))+pi, mesh_size(1));
            else
                vec_x = linspace(obj.x_eq(ind1)-obj.delta_x(ind1), obj.x_eq(ind1)+obj.delta_x(ind1),mesh_size(1));
            end
            if angle2
                if any(size(ind2) > 2)
                    ind2 = reshape(ind2,1,3);
                else
                    ind2 = [reshape(ind2,1,2) 1];
                end
                vec_y = linspace(obj.angle_eq(ind2(3))-pi, obj.angle_eq(ind2(3))+pi, mesh_size(2));
            else
                vec_y = linspace(obj.x_eq(ind2)-obj.delta_x(ind2),obj.x_eq(ind2)+obj.delta_x(ind2), mesh_size(2));
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
                        X(ind1) = (plot_x(i,j) - obj.x_eq(ind1)) * obj.D(ind1,ind1);
                    end
                    if angle2
                        [X(ind2(1)) ,X(ind2(2))] = deal(sin(plot_y(i,j) - obj.angle_eq(ind2(3))) ,...
                        (1 - cos(plot_y(i,j) - obj.angle_eq(ind2(3)))) / 2);
                    else
                        X(ind2) = (plot_y(i,j) - obj.x_eq(ind2)) * obj.D(ind2,ind2);
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
            contour(plot_x, plot_y, Z, [1 1], col,'linewidth', 2, 'DisplayName', [app_str ' RoA approximation: T=' num2str(obj.T) ', d=' num2str(obj.d)]);
            if angle1
                coef1 = 1;
                xeq = obj.angle_eq(ind1(3));
            else
                coef1 = obj.A(ind1,ind1);
                xeq = obj.x_eq(ind1);
            end
            if angle2
                coef2 = 1;
                yeq = obj.angle_eq(ind2(3));
            else
                coef2 = obj.A(ind2,ind2);
                yeq = obj.x_eq(ind2);
            end
            parametric_angle = linspace(0,2*pi,100);
            if plot_target
                plot(xeq + obj.epsilon/coef1 * cos(parametric_angle), yeq + obj.epsilon/coef2 * sin(parametric_angle), '-k','linewidth', 2, 'DisplayName','Target set');
            end
            %ZZ = (coef1*plot_x).^2+(coef2*plot_y).^2;
            %contour(plot_x, plot_y, ZZ, [obj.epsilon obj.epsilon], '-k','linewidth', 2);
            % rectangle('position',[-1+ rho_eq -pi+theta_eq 2 2* pi]);
            hold off
            %axis equal
            %title('Region of attraction over X = [- 1, 1]')
            legend('Location','southwest')
            
            title(strjoin(['ROA(',str1,',',str2,')']))
            xlabel(str1);
            ylabel(str2);
        end
    
        function [] = plot_w(obj, ind1, ind2, approximation, str1, str2, mesh_size)
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

            if nargin < 7
                mesh_size = [40,40];
            end
			if nargin < 5
                str1 = ['x_{' num2str(ind1) '}'];
                str2 = ['x_{' num2str(ind2) '}'];
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
                vec_x = linspace(obj.x_eq(ind1)-obj.delta_x(ind1), obj.x_eq(ind1)+obj.delta_x(ind1),mesh_size(1));
            end
            if angle2
                if any(size(ind2) > 2)
                    ind2 = reshape(ind2,1,3);
                else
                    ind2 = [reshape(ind2,1,2) 1];
                end
                vec_y = linspace(obj.angle_eq(ind2(3))-pi, obj.angle_eq(ind2(3))+pi, mesh_size(2));
            else
                vec_y = linspace(obj.x_eq(ind2)-obj.delta_x(ind2),obj.x_eq(ind2)+obj.delta_x(ind2), mesh_size(2));
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
                        X(ind1) = (plot_x(i,j) - obj.x_eq(ind1)) * obj.D(ind1,ind1);
                    end
                    if angle2
                        [X(ind2(1)) ,X(ind2(2))] = deal(sin(plot_y(i,j) - obj.angle_eq(ind2(3))) ,...
                        (1 - cos(plot_y(i,j) - obj.angle_eq(ind2(3)))) / 2);
                    else
                        X(ind2) = (plot_y(i,j) - obj.x_eq(ind2)) * obj.D(ind2,ind2);
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

            title(strjoin(['w(',str1,',',str2,')']))
            xlabel(str1);
            ylabel(str2);
        end

        function [] = plot_v(obj, ind1, ind2, tau, approximation, str1, str2, mesh_size)
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

            if nargin < 8
                mesh_size = [40,40];
            end
			if nargin < 6
                str1 = ['x_{' num2str(ind1) '}'];
                str2 = ['x_{' num2str(ind2) '}'];
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
                vec_x = linspace(obj.x_eq(ind1)-obj.delta_x(ind1), obj.x_eq(ind1)+obj.delta_x(ind1),mesh_size(1));
            end
            if angle2
                if any(size(ind2) > 2)
                    ind2 = reshape(ind2,1,3);
                else
                    ind2 = [reshape(ind2,1,2) 1];
                end
                vec_y = linspace(obj.angle_eq(ind2(3))-pi, obj.angle_eq(ind2(3))+pi, mesh_size(2));
            else
                vec_y = linspace(obj.x_eq(ind2)-obj.delta_x(ind2),obj.x_eq(ind2)+obj.delta_x(ind2), mesh_size(2));
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
                        X(ind1) = (plot_x(i,j) - obj.x_eq(ind1)) * obj.D(ind1,ind1);
                    end
                    if angle2
                        [X(ind2(1)) ,X(ind2(2))] = deal(sin(plot_y(i,j) - obj.angle_eq(ind2(3))) ,...
                        (1 - cos(plot_y(i,j) - obj.angle_eq(ind2(3)))) / 2);
                    else
                        X(ind2) = (plot_y(i,j) - obj.x_eq(ind2)) * obj.D(ind2,ind2);
                    end
                    K(i,j) = eval(k);
                    L(i,j) = 1;
                end
            end
            figure(3)
            hold on
            mesh(plot_x, plot_y, double(K));
            plot3(plot_x, plot_y, double(L) ,'k');
            view([37.5 30]);
            hold off

            title(strjoin(['v(',str1,',',str2,',',tau,')']))
            xlabel(str1);
            ylabel(str2);
        end
    end
end