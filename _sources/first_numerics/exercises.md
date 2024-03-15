# Exercises

## Exercise 0

Draw a `unit_cube` using the `webgui` of NGSolve utilizing one of the methods to run NGSolve described [here](first_numerics)

## Exercise 1 

Read {numref}`mol` and {numref}`basic_fe_wave` (more information can be found in the [iFEM](https://jschoeberl.github.io/iFEM) tutorial (especially [the basic example](https://jschoeberl.github.io/iFEM/primal/first_example.html) )) and answer the following:
  - For a given strong formulation of a time-domain wave problem (e.g., {eq}`wave_2o`, describe in at least 6 steps how one approximates the solution of said problem using a method of lines approach with finite elements (assuming you have a finite element code like NGSolve readily available).
  - Which of these steps are inherent to the fact that we treat a time-domain problem using the method of lines (compared to finite elements for stationary problems)?
  - Which of these steps are inherent to the fact that we specifically use a finite element method (compared to a general Galerkin method)?
  - Which choices do you have to make in each step?

## Exercise 2

Add a parameter `order = n` for $n\in\mathbb N$ in the initialization of the `H1` finite element space and look at the resulting basis functions (default `order = 1` if the parameter is ommited). Write down at least 4 characterizing properties of the basis functions. Indicate which of this properties characterize the space $V$ and which specifically the choice of the basis.

## Exercise 3

Create a geometry (two or three-dimensional) and mesh where you would expect interesting acoustic effects (e.g., some (very simplified) musical instrument, ...).
Simulate the wave problem on your geometry/mesh. 


## Exercise 3

Verify that the implementation of the Newmark time-stepping scheme in {numref}`basic_fe_wave` is actually consistent with the formulae presented above.
How does it have to be modified if $B\neq 0$ or $r\neq 0$?

Change your code to do the following
  - Add time-harmonic source term of the form $f(t,x)=\cos(\omega t) f_0(x)$ for some $\omega>0$ and $f_0:\Omega\to \mathbb R$ with small support.
  - Add a first order absorbing boundary condition to (at least) one of the boundaries.

## Exercise 4

Derive a weak formulation for the 2nd order time domain Maxwell system {eq}`maxwell_2o`. What are the natural, homogeneous boundary conditions for this system?
Look at the basis functions of the space `HCurl` (instead of `H1`). Use the option `vectors = True` in the `Draw` command. How do you think these basis functions are constructed?
Simulate a time-domain electromagnetic wave on the  `unit_cube` with initial data of your choice and zero right hand side and homogeneous natural boundary conditions.
