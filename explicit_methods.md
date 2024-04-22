(explicit_methods)=
# Explicit methods for time-domain waves

In {numref}`time_integration` we have defined **explicit** time-stepping methods as methods where in each time step merely the inverse of one or two mass matrices has to be applied. These methods are particularly efficient when the inverse mass matrice(s) are cheap to compute/apply. To construct spaces featuring such matrices is the goal of this section.

We present two different approaches to construct such spaces. Notably both of these methods leave the setting of classical Galerkin methods.
