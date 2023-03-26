# SOStab

A sum-of-squares Matlab toolbox for stability analysis. The goal of the toolbox is to facilitate the use of SoS programming for calculating an approximation of the region of attraction of a dynamical system. The only input needed are an equilibrium point, the admissible set and the dynamic of the system.

SOStab requires [YALMIP](yalmip.github.io/), as well as a semidefinite solver. [Mosek](mosek.com/) is used by default, but it can be replaced by any other solver, provided they are installed and interfaced through YALMIP. 


# Installation 

## Simple example 
### Initialize the class

Initialize with an equilibrium point `x_eq` and a range `delta_x`, that defines the admissible set: `[x_eq-delta_x, x_eq+delta_x]`. The two vectors must have the same size. You can change the default values for the solver used and the verbose parameter of the optimization calls.

```
toy = SOStab([0;0],[1;1]);
%toy.solver = 'mosek';
%toy.verbose = 2;
```

### Define the dynamic of the system

Define f such that \dot{x} = f(x), the vector must have the same size as x_eq. The property `x` of the class is a Yalmip sdpvar, which is a symbolic representation of the state variable. It has the same size as the problem. The property `t` is also an sdpvar, representing time that can be useful here.
For example, here we implement a simple radial system:

```
toy.dynamics = (toy.x'*toy.x - 0.25)*toy.x;
```

### Solve the outer approximation problem

Solve the outer approximation of the finite horizon ROA for a target set {x : ||x-x_eq||< epsilon}, with `epsilon=0.1` and a time horizon `T=20`. The degree of the polynomials is defined by `d=6`. `d` above 20 is useless and for high dimensions (n > 6), a standard computer RAM can't cope with `d` above 10.

```
toy.SoS_out(6,20,0.1);
```


### Solve the inner approximation problem

Solve the inner approximation of the ROA. Note that the results of this approximation are very sensitive to the relative size of the sets and the degree of the polynomials. You can check the relevance of the results by verifying the coefficients of `w`. If they are null for most of them, then the computation didn't worked. You can try to increase the size of the target set or the degree of the polyomials.

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

### Plot the ROA

Plot the projection of the ROA on two dimensions (for the last calculated approximation).

```
toy.plot_roa(1,2,'outer');
axis('equal');
```
The first two arguments indicate the indices of the variables to plot. The fourth (optionnal) argument indicates to plot the target set.

### Plot the surface of w

Plot the surface of w for the last calculated inner or outer approximation.

```
toy.plot_w(1,2, 'o');
```