## All code

>functions.py

The primary script containing most of the work done. Contains all kinds of different functions used in all other python scripts.

Every function has accompanying documentation and variable explanation, plus comments wherever I felt it was necessary.

>main.py

Not really used that much in this project. It contains two different functions for calculating and solving a) and e) respectively. 

The reasoning as to why c, d and f are not called from here is explained in a comment at the top of the script.

Run from the console:

Example: proj_a()

>probability.py

Essentially all of task e). Takes an input for T, temperature, L, lattice length and cycles for Monte Carlo cycles

Plots a normalized histogram for probabilities of finding different energies after it has stabilized/reach equilibrium

Call either in main.py or just directly here.

>runner.py

This is the big bad function. Takes care of task f). Parallizes the code appropriately and dynamically based on assigned processor threads and solves.

When it is done, it returns a set of plots. Has a set of booleans that dictate what to return and calculate.

Run normally in the console

>stabilizer.py

Similar to runner.py, but instead cycles through different numbers of Monte Carlo cycles to examine when we reach a steady state/equilibrium.

>timer.py

For reading our data and plotting the time spent calculating as a function of CPU threads.