## This is all the code used in project 1. It's all written in Python 3.8.
>main.py 

main.py is the most flexible "all-around" program. It's neither built for speed or very short and to the point, but it can be used for anything regarding the Tridiagonal Matrix Algorithm

>general_solution.py

general_solution.py is a program made to solve a TMAD problem generally. Goes from n = 10 to n = 10000000

>special_solution.py

special_solution.py solves our special symmetrical TMAD problem. Goes from n = 10 to n = 10000000

>LU.py

LU.py solves our system on linear equations using LU-factorization. Goes from n = 10 to n = 10000

>plotter.py

plotter.py plots the benchmark times and the error. Requires full runs through n = 10 to n = 10000 / 10000000. Because this was a last minute addition it's not very flexible and thus not compatible with more than thoss full runs.
