
# Two Dimensional Motion, and Quadratic Air Resistance
## Investigation of falling bodies using raylib, and my own RK4 Numerical Integrator
==================================================================================
In my Falling_Bodies program, I demonstrate how in a vaccuum, objects fall with the same acceleration while no resistance is present

In this program, there is a non-conservative force at play, (drag) that changes the terminal velocity based on the mass (and shape) of the object

Here, spheres will be used, and we are assuming Standard Temperature and Pressure (STP) which is the standard temperature and pressure at the surface of the Earth.

## Why is the air resistance quadratic here? Why not Linear? Even when |v| == 0? 

One can see which force (linear or quadratic drag) dominates, or if both need to be considered by taking the ratio and making a decision.

From Taylor's Classical Mechanics, page 44 and 45.

Assuming the drag force is always in the exact opposite direction of the velocity, and the object is a sphere, we can make the following statements

$\large \vec{f} = -f(v)\mathbf{\hat{v}}$

$\large f(v) = bv + c v^2 = f_{lin} + f_{quad}$

$\large b = \beta D$ $\hspace{1cm}$ $\large c = \gamma D^2$

Where

$\large \beta = 1.6 \times 10^{-4}   N \dot s / m^2$  $\hspace{1cm}$  $\large \gamma = 0.25   N \dot s^2 / m^4$

$\large \frac{f_{quad}}{f_{lin}} = \frac{cv^{2}}{bv} = \frac{\gamma D}{\beta}v = (1.6 \times 10^{3} \frac{s}{m^2}) D v$ 

Since the ratio is such a large number for essentially all values (except where |v| = 0), this means that the quadratic force dominates

So for the split second where the velocity in the y direction is 0, this should be treated as linear, but we are interested in smooth, continuous equations for now,
we will be treating these as smooth functions, with quadratic drag only.

And for most of the flight of the sphere, this will be quadratic.

## The Equations of Motion for Quadratic Drag

$$

m \ddot{x} = -c \sqrt{\dot{x^2} + \dot{y^2}} \dot{x} \\

m \ddot{y} = -c \sqrt{\dot{x^2} + \dot{y^2}} \dot{x}

$$

## Varying Mass

Here, all of the diameters of the spheres are 1 meter, and the masses vary.

$\m_{b} = 0.1 kg$ $\hspace{1cm}$ $\m_{r} = 10.0 kg$ $\hspace{1cm}$ $\m_{o} = 1000.0 kg$

![Varying Mass](plots/Pos_Quad_AR_VaryingMass.png)

As you can see, the orange sphere, which is much more massive, takes a longer time to reach it's terminal velocity
