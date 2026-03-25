File Included:

INSTRUCTIONS.md
setup_new.jl - script that will install and activate all relevant modules
Project.toml - file that tracks julia packages. Use by setup_new.jl
q1.jl - Julia implementation of Q1 - not covered in report
q2.jl - Julia file used to produce Q2 graph
q3.jl - Julia file with implementations of Q3. Contains q3_a(), q3_b() and q3_c() described below.
q3_tfunctions.jl - Julia file with functions used to test optimisation methods
q3_tests.jl - Julia file that has various tests and benchmarks of optimisation methods



============================ Question 1 ============================




============================ Question 2 ============================

Contains an implementation of our extension to the labour-leisure model.
Code is included that shows how the graph in q2 was generated. It was run via plot_utility_variation(param = :epsilon, values = [0.0, 0.5, 1.0, 1.5, 2.0])
The default parameters were hand picked so the plot clearly showed how it varied with different values.
The plotting function includes code that optimises the value of h*, and shows evaluating utility at ĥ for different epsilon values, which could be implemented in solving the model.

============================ Question 3 ============================

------------------------ Python Implementation ------------------------


------------------------ Julia Implementation ------------------------

Step 1 - run julia in the shell to enter the REPO

Step 2 - run include("setup_new.jl") to install and activate all used packages and include all relevant modules.

Step 3 - run q3_a() - this will print the results found by 2 different methods in the command line.

Step 4 - run q3_b() - this will print h(2) & h(3) in the command line and generate the required graph. M is a keyword argument and can be
         altered if desired by using q3_b(M = 1000). The default value is 100.

Step 5 - run q3_c(function) - this requires the chosen function as an argument. Z can be passed in here as it is defined in the module e.g. q3_c(Z).
         The function will print the global max found using 2 different algorithms. grid_size is a keyword argument and can be varied if desired - 1000
         is the default value e.g. q3_c(Z, grid_size = 100)

**Please note this is simply how to replicate the questions assigned. Other tests, benchmarks and interesting are included and described in the files that can be viewed to better understand the development and testing process. They have shown significant robustness of the methods. 