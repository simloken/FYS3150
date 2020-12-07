import numpy as np
import matplotlib.pyplot as plt

class Simulation:
    
    def __init__(self,
            N,
            NN,
            MCS,
            plotting=False):
        
        self.N = N
        self.NN = NN
        self.MCS = MCS
        self.plotting = plotting
        
        
    def run(self):
        Alst = []; Blst = []; ABlst = [];
        self.V = np.random.choice([-1,1], size=(self.N)) #either up or down
        if self.plotting==True:
            plt.plot((range(len(self.V))), self.V)
            plt.show()
       
        
        
        self.V, self.m = self.dynamicRules(self.V, self.MCS) #run over with MCS montecarlo steps
       
        
        if self.plotting==True:
            plt.plot((range(len(self.V))), self.V)
            plt.show()
        A = 0; B = 0; AB = 0
        for i in range(self.NN): #pick random spot NN times to gauge spin consensus
            ridx = np.random.randint(self.N-2)
            lst = []
            lst.append(self.V[ridx-1])
            lst.append(self.V[ridx-2])
            lst.append(self.V[ridx+1])
            lst.append(self.V[ridx+2])
            if 1 not in lst: #if all -1
                B += 1
            elif -1 not in lst: #if all +1
                A += 1
            elif 1 in lst and -1 in lst: #if both
                AB += 1
    
        Alst.append(A); Blst.append(B); ABlst.append(AB)
            
        A = np.mean(Alst)
        B = np.mean(Blst)
        AB = np.mean(ABlst)
        A /= self.NN
        B /= self.NN
        AB /= self.NN
        print("State 1: {:.2f} (AAAA)\nState 2: {:.2f} (BBBB)\nState 3: {:.2f} (ABAB)".format(A, B, AB))
        
        return A, B, AB, self.m
        
        
        
        
        
        
        
        
    def dynamicRules(self, V, MCS):
        m = []
        for i in range(MCS):
            N = len(V)
            ridx = np.random.randint(N-2)
            if ridx != 0:
                if V[ridx] == V[ridx+1]:
                    V[ridx-1] = V[ridx]; V[ridx+2] = V[ridx]
                elif V[ridx] == -V[ridx+1]:
                    V[ridx-1] = V[ridx+1]; V[ridx+2] = V[ridx]
            
            elif ridx == 0:
                if V[ridx]*V[ridx+1] == 1:
                    V[ridx+2] = V[ridx]
                elif V[ridx]*V[ridx+1] == -1:
                    V[ridx+2] = V[ridx]
                    
            m.append(self.magnet())
        
        return V, m
    
    
    def magnet(self):
        m = (1/len(self.V)) * np.sum(self.V)
        return m