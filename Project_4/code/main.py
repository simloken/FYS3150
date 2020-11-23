from functions import analytics
from probability import proba

"""
Here are a few function that directly solve a) and e) respectively.
As Python multiprocessing requires the code to be within __main__ on windows to work,
I couldn't find a way to parse arguments from here to the appropriate script
Therefore, for demonstrations on how to run/examine c, d and f, see either it's
appropriate readme on github or just change variables in the code and run it
like so. It's pretty intuitive.
"""

def proj_a(T=1):
    L = 2
    anas = analytics(T)/L**2
    print('Analytical Values')
    print('E: %g\n|M|: %g\nCv: %g\nX: %g\n' %(anas[0],anas[1],anas[2],anas[3]))    
    
    
def proj_e(T=1, L=20, cycles=100000):
    proba(T,L,cycles)