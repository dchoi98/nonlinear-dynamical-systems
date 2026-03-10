# Nonlinear Dynamical Systems

Derrick Choi

Differential equations relate functions to their derivatives. They show up in fields ranging from engineering and physics to economics and biology. When the function or its derivatives appear raised to powers greater than one, or are multiplied together, the equation is nonlinear. Nonlinear differential equations often describe chaotic systems, where tiny differences in initial conditions lead to very different outcomes over time. This notebook walks through several examples of systems governed by nonlinear differential equations, deriving the equations of motion from the Lagrangian and solving them numerically in Mathematica.

## Numerical Methods for Differential Equations

Most nonlinear differential equations can only be solved analytically for special initial conditions, so we usually have to solve them numerically. The simplest approach is Euler's method. Starting from a known point, you compute the slope of the tangent line there, step forward by some small increment *h*, and use the slope to estimate the function's value at the new point. Then you repeat. More steps over a given interval means a smaller *h* and a better approximation. The interactive plot in the notebook shows how the approximation improves as the number of steps increases for the ODE *y'(x) = x²y(x) − (3/2)y(x)*, *y(0) = 1*.

<p align="center">
  <img src="imgs/eulers\_method.gif" alt="Euler's method converging to the exact solution" width="500">
</p>

Euler's method is fast but not particularly accurate. The higher-order methods used by `NDSolve` are more precise, but the underlying idea of stepping forward using local derivative information is the same.

## Double Pendulum

A double pendulum is a pendulum with a second pendulum hanging from its end. The motion is governed by a pair of coupled nonlinear differential equations, and the system is chaotic.

To derive the equations of motion, we use the Lagrangian formalism. We define the positions of the centers of mass of the two pendulum rods in terms of the generalized coordinates *θ₁* and *θ₂* (the angles each rod makes from the vertical, with *θ* = 0 pointing straight down). From these positions we compute the velocities by differentiating with respect to time, then build the kinetic energy *T* as the sum of translational and rotational kinetic energy for each rod:

$$T = \\tfrac{1}{2} m v\_1^2 + \\tfrac{1}{2} I \\dot{\\theta}\_1^2 + \\tfrac{1}{2} m v\_2^2 + \\tfrac{1}{2} I \\dot{\\theta}\_2^2$$

where *I = (1/12)ml²* is the moment of inertia of a uniform rod about its center of mass. The potential energy is solely due to gravity:

$$V = m g y\_1 + m g y\_2$$

The Lagrangian is *L = T − V*. Applying the Euler–Lagrange equations for each generalized coordinate gives two coupled second-order ODEs in *θ₁* and *θ₂*, which are solved numerically with `NDSolve`.

<p align="center">
  <img src="imgs/double\_pendulum.gif" alt="Double pendulum simulation" width="450">
</p>

The motion of the double pendulum is highly sensitive to initial conditions. The two simulations below use identical parameters and starting conditions for the inner pendulum. The outer pendulums start with the same angular velocity and a difference in angle of only *Δθ* = 0.0001. At first the two systems move in sync, but they quickly diverge until their motions look completely unrelated.

<p align="center">
  <img src="imgs/diverging\_double\_pendulums.gif" alt="Two double pendulums diverging from nearly identical initial conditions" width="450">
</p>

## Swinging Atwood's Machine

A swinging Atwood's machine consists of two masses connected by an inextensible string that passes over two pulleys. One mass (the "swinging" mass, *m* = 1) can swing freely in a plane, while the other (the "counterweight" mass, *m* = *μ*) is constrained to move vertically.

The system has two degrees of freedom: *r*, the distance of the swinging mass from its pulley, and *θ*, the angle the string makes with the downward vertical. The kinetic energy has contributions from both masses:

$$T = \\tfrac{1}{2} m (\\dot{r}^2 + r^2 \\dot{\\theta}^2) + \\tfrac{1}{2} \\mu \\dot{r}^2$$

The swinging mass moves in polar coordinates, so its speed squared is *ṙ² + r²θ̇²*. The counterweight moves along the string, so its speed is just |*ṙ*|. The potential energy is:

$$V = -m g r \\cos\\theta + \\mu g r$$

Forming *L = T − V* and applying the Euler–Lagrange equations for *r* and *θ* gives two coupled second-order ODEs, solved numerically with `NDSolve`.

<p align="center">
  <img src="imgs/swinging\_atwoods\_machine.gif" alt="Swinging Atwood's machine animation for μ = 4.5" width="450">
</p>

For certain mass ratios the swinging mass traces out complex but periodic orbits. The figures below show the trajectories for various values of *μ*, each producing a distinct pattern.

<p align="center">
  <img src="imgs/mu\_3.jpg" alt="Swinging Atwood's machine trajectory for μ = 3" width="450">
</p>

<p align="center">
  <img src="imgs/mu\_5\_6\_16\_19.jpg" alt="Swinging Atwood's machine trajectories for μ = 5, 6, 16, and 19" width="600">
</p>

## Three-Body Problem

The three-body problem is the problem of taking the initial positions and velocities of three point masses and calculating their trajectories using Newton's laws of motion and gravitation. There is no general closed-form analytic solution, and the system is generically chaotic.

Each mass has three degrees of freedom, giving nine for the system in total. The kinetic energy is the sum of the individual kinetic energies:

$$T = \\tfrac{1}{2} m\_1 |\\dot{\\mathbf{r}}\_1|^2 + \\tfrac{1}{2} m\_2 |\\dot{\\mathbf{r}}\_2|^2 + \\tfrac{1}{2} m\_3 |\\dot{\\mathbf{r}}\_3|^2$$

The potential energy is the sum of the gravitational potential energies between each pair:

$$V = -\\frac{G m\_1 m\_2}{|\\mathbf{r}\_1 - \\mathbf{r}\_2|} - \\frac{G m\_1 m\_3}{|\\mathbf{r}\_1 - \\mathbf{r}\_3|} - \\frac{G m\_2 m\_3}{|\\mathbf{r}\_2 - \\mathbf{r}\_3|}$$

With *L = T − V*, the Euler–Lagrange equation for each coordinate produces nine coupled second-order ODEs, solved numerically with `NDSolve`.

<p align="center">
  <img src="imgs/three\_body.gif" alt="Three-body problem simulation" width="450">
</p>

