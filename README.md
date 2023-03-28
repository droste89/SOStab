# SOStab

A sum-of-squares Matlab toolbox for stability analysis. The goal of the toolbox is to facilitate the use of SoS programming for calculating an approximation of the region of attraction (RoA) of a dynamical system. The only input needed are an equilibrium point, the admissible set and the dynamics of the system.

SOStab requires [YALMIP](yalmip.github.io/), as well as a semidefinite solver. [Mosek](mosek.com/) is used by default, but it can be replaced by any other solver, provided they are installed and interfaced through YALMIP. 

This Readme illustrates how to use the class with several examples and then describe the properties and the method of the class in the [last section](#properties-and-methods-of-the-class).


# Installation

Only the `SOStab.m` file is needed for the toolbox to work. You can put it in your working directory or add it to your Matlab path, to have it available everywhere:

```
addpath SOStab
```

# Simple example 
## Initialize the class

Initialize with an equilibrium point $x_{eq}$ and a range $\Delta x$, that defines the admissible set: $[x_{eq}-\Delta x, x_{eq}+\Delta x]$. The two vectors must have the same size. You can change the default values for the solver used and the verbose parameter of the optimization calls.

```
toy = SOStab([0;0],[1;1]);
%toy.solver = 'mosek';
%toy.verbose = 2;
```

## Define the dynamics of the system

Define $f$ such that $\dot{x} = f(x)$, the vector must have the same size as $x_{eq}$. The property `x` of the class is a YALMIP sdpvar, which is a symbolic representation of the state variable. It has the same size as the problem. The property `t` is also an sdpvar, representing time that can be useful here.
For example, here we implement a simple radial system:

```
toy.dynamics = (toy.x'*toy.x - 0.25)*toy.x;
```

## Solve the outer approximation problem

Solve the outer approximation of the finite horizon RoA for a target set $\{x : ||x-x_{eq}||< \varepsilon\}$, with $\varepsilon=0.1$ and a time horizon $T=20$. The degree of the polynomials is defined by $d=6$. $d$ above 20 is useless and for high dimensions ($n > 6$), a standard computer RAM can't cope with $d$ above 10.

```
toy.SoS_out(6,20,0.1);
```


## Solve the inner approximation problem

Solve the inner approximation of the RoA. Note that the results of this approximation are very sensitive to the relative size of the sets and the degree of the polynomials. You can check the relevance of the results by verifying the coefficients of $w$. If they are null for most of them, then this means that the solver was not able to find the optimum of the problem. You can try to increase the size of the target set or the degree of the polyomials.

```
toy.SoS_in(6,20,0.1);
witoy = value(toy.wcoef_inner);
witoy(1) = witoy(1)-1;
if witoy'*witoy > 0.00001
    toy.plot_roa(1,2,'inner',0,"x","y");
else
    disp("w not relevant")
end
```

## Plot the RoA

Plot the projection of the RoA in two dimensions (for the last calculated approximation).

```
toy.plot_roa(1,2,'outer');
axis('equal');
```
The first two arguments indicate the indices of the variables to plot. The fourth (optionnal) argument indicates to plot the target set.

## Plot the surface of $w$

Plot the surface of $w$ for the last calculated inner or outer approximation.

```
toy.plot_w(1,2, 'o');
```



# Various other examples

## Van der Pol

The Van der Pol oscillator is a simple 2 dimensional system:
$$ \begin{pmatrix} \dot{x}_1 \\ \dot{x}_2 \end{pmatrix} = \begin{pmatrix} -2x_2 \\ 0.8x_1 + 10(x_1^2-0.21)x_2 \end{pmatrix} $$

```
vdp = SOStab([0;0],[1.1;1.1]);
vdp.dynamics = [-2*vdp.x(2); 0.8*vdp.x(1) + 10*(vdp.x(1)^2-0.21)*vdp.x(2)];
vdp.SoS_out(12,1,0.5);
vovdp = value(vdp.vcoef_outer);
wovdp = value(vdp.wcoef_outer);
vdp.plot_roa(1,2,'outer',1,"x_1","x_2");
vdp.SoS_in(12,1,0.5);
vivdp = value(vdp.vcoef_inner);
wivdp = value(vdp.wcoef_inner);
wivdp(1) = wivdp(1)-1;
if wivdp'*wivdp>0.00001
    vdp.plot_roa(1,2,'inner',0,"x_1","x_2");
end
```

## Scaled pendulum

The pendulum is a simple example of "polynomialization" of a system.
$$\dot{\theta} = \sin(2\theta)=2\sin\theta\cos\theta $$
The system used as input of the toolbox will be $x =(\sin\theta, \cos\theta, \omega) $ (here $\omega=0$ is used for plotting in 2 dimensions), then we have:

$$\dot{x}= \begin{pmatrix}-2\sin\theta\cos^2\theta \\ 2\sin^2\theta\cos\theta \\ 0 \end{pmatrix}$$

Hence, the dynamic is indeed polynomial in the variables.

```
pen = SOStab([0;1;0],[1;1;1],[1,2]);
pen.dynamics = [-2*pen.x(1)*pen.x(2)^2;2*pen.x(1)^2*pen.x(2);0];
pen.SoS_out(6,40,0.5);
vopen = value(pen.vcoef_outer);
wopen = value(pen.wcoef_outer);
pen.plot_roa([1,2,1],3,'outer',1,"\theta","\omega");
pen.SoS_in(6,40,0.5);
if double(pen.wcoef_outer(abs(double(pen.wcoef_outer)) > 0.00001))
    vipen = value(pen.vcoef_inner);
    wipen = value(pen.wcoef_inner);
    pen.plot_roa([1,2],3,'inner',0,"\theta","\omega");
else print("w not relevant")
end
```

## Power system

Consider an electrical power network made of synchronous
machines connected in a cycle, along with its second reduced
order model (see M. Anghel, F. Milano, and A. Papachristodoulou, "Algorithmic construction of Lyapunov functions for power system stability analysis"):

$$    \begin{align*}
        \dot \theta_1 &= \omega_1, \qquad \dot \theta_2 = \omega_2 \\
        \dot \omega_1 &= -\sin\theta_1 - 0.5\sin(\theta_1-\theta_2) - 0.4 \, \omega_1 \\
        \dot \omega_2 &= -0.5\sin\theta_2 - 0.5\sin(\theta_2-\theta_1) - 0.5 \, \omega_2 + 0.05
    \end{align*}
$$

Similarly to the pendulum, it can be polynomialized we the change of variable $x= (\sin\theta_1, \cos\theta_1, \sin\theta_2, \cos\theta_2,\omega_1, \omega_2)$. The RoA of the system is calculated by:

```
eq = [sin(0.02);cos(0.02);sin(0.06);cos(0.06);0;0];
dev = [1;1;1;1;pi;pi];
ang_ind = [1,2;3,4];
pow_sys = SOStab(eq,dev,ang_ind);
pow_sys.dynamics = [pow_sys.x(5)*pow_sys.x(2);-pow_sys.x(5)*pow_sys.x(1);...
    pow_sys.x(6)*pow_sys.x(4);-pow_sys.x(6)*pow_sys.x(3);
    -pow_sys.x(1)-0.5*(pow_sys.x(1)*pow_sys.x(4)-pow_sys.x(2)*pow_sys.x(3))-0.4*pow_sys.x(5);
    -0.5*pow_sys.x(3)+0.5*(pow_sys.x(1)*pow_sys.x(4)-pow_sys.x(2)*pow_sys.x(3))-0.5*pow_sys.x(6)+0.05];
pow_sys.SoS_out(10,8,0.1);
vopow = value(pow_sys.vcoef_outer);
wopow = value(pow_sys.wcoef_outer);
pow_sys.plot_roa([1,2],[3,4],'o','1',"\theta_1","\theta_2");
```

# Properties and methods of the class

Properties of the class are the following:
- `dimension`: dimension of the problem (number of variables)
- `x_eq`: first argument of the call, the equilibrium state of the system
- `delta_x`: the range around the equilibrium defining the feasible set
- `angle_eq`: the vector $\theta_{eq}$ - the equilibrium angles - recalculated from the values of their sine and cosine (empty if no phases are involved)
- `angle_ind`: the indices of sine and cosine functions in the recast variable $x$, given as last input of `SOStab` call (empty if no phases are involved)
- `x`: a YALMIP sdpvar polynomial object, of the dimension of the problem. It represents the variable $x$ and is used by the user to define the dynamics of the system
- `t`: sdpvar polynomial of size 1, representing the time variable, which can be needed to define the dynamics of the system (if non-autonomous)
- `D`, the matrix $D$ of the variable change used in the toolbox to normalize the system
- `invD`, the inverse $D^{-1}$ of the matrix $D$
- `solver`, the choice of the solver used in the optimization, defined as Mosek by default
- `verbose`, the value of the verbose parameters of the YALMIP optimization calls, defined at 2 by default
- `dynamics`, a YALMIP polynomial defining the polynomial dynamics $f$ of the system

The class has seven methods:

- The class initialization `SOStab` takes as input the equilibrium point $x_{eq}$, the distance $\Delta x$ as two vectors, that must have the same length -- the dimension of the problem. It initializes an object of the class `SOStab`.
- `moments` is a method called inside the other methods, which calculates the integration of all the monomials on the normalized feasible set $[-1,1]^n$. The function has only one argument, the maximum degree of the monomials. It returns the vector $y$ such that the integration of a polynomial $p$ on the admissible set $X$ would be $\int_X p = y^\top v_p$ where $v_p$ is the vector of coefficients of $p$ in the basis of monomials (with Yalmip ordering convention if $n\geq 2$).
- `SoS_out` takes as input an optional matrix $A$ (by default $A = I_n$) and a real positive value $\varepsilon$ - defining the target $\mathbb{M} = \{x \in \mathbb{R}^n : \| A(x-x_{eq})\|^2 \leq \varepsilon^2 \}$ -, an even degree $d=2k$ for the polynomials and a time horizon $T$. It solves the strengthening of the optimization problem giving the RoA in bounded degree $d$ and returns its optimal value (an over approximation of the volume of the calculated RoA) and the coefficients of the polynomials $v^\star_k$ and $w^\star_k$ defining this approximation.
- `SoS_in` takes the same inputs and solves the counterpart of the previous problem, which is the inner approximation optimization problem. It returns the same outputs as the previous method. Currently, the inner approximation has some limitations on some of the instances tested, due to the algebraic representation of the boundary of the admissible state.
- `plot_roa` takes as inputs the two indices $i,j$ of the variables on which to project the ROA, a string which states which approximation to plot (outer or inner) and four optional arguments: an `int` indicating whether or not to plot the target set (1 or 0) - 0 by default -, two strings indicating the name of the axis (respectively for the first and second indices) - $x_i$ and $x_j$ by default - and a vector giving the size of the plotting mesh - $(40,40)$ by default. It plots the expected slice of the ROA, with all other variables at equilibrium. 

    If both inner and outer approximations are called sequentially for the same variables, the two plots will appear on the same figure.
    
    If an angle $\theta$ is one of the plotted variable, the index argument should be a vector of 2 values: the index of $\sin\theta$ in the variables and the one of $\cos\theta$.
- `plot_v` and `plot_w` take the same inputs. They respectively plot the graph of $v^\star_k$ and $w^\star_k$ in 3D (with non-selected variables at equilibrium).
\end{itemize}
