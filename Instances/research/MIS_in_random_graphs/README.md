These problem instances are intended for further research of hard problems within the class of NP-hard problems.
It seems that some NP-hard problems are harder then other. These instances could lead to a significant slowdonw of (exact) optimization solvers.
The idea here is to identify properties of such problem instances.

Such instances could be from the class of random graphs with edge density 0.5 for the maximum independent set problem (MIS) and for the fixed size independet set problem (IS). 
Two such instances seems to be promissing, so far. The computations were performed on a laptop with 13th Gen Intel(R) Core(TM) i7-1370P CPU.

- G_145_0.5_0
    - An instance of a random graph on 145 vertices with edge density 0.5 (Each edge was generated independently with a probability 0.5). MIS defined as BQP problem with an "artificial" equality constrain.
- K_10_G_145_0.5_0
    - This the problem instance on the same graph as in G_145_0.5_0 but with IS of size exactly 10 as BQP problem.





