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
        self.times = np.zeros(N)
        self.dT = []
        self.time = 0
        
    def run(self):
        Alst = []; Blst = []; ABlst = [];
        self.V = np.random.choice([-1,1], size=(self.N)) #either up or down
        if self.plotting==True:
            plt.figure(figsize=(8, 1))
            plt.title('Initial Spread')
            plt.plot((range(len(self.V))), self.V, color='k')
            plt.show()
        
        
        self.V, self.m = self.dynamicRules(self.V, self.MCS) #run over with MCS montecarlo steps
       
        
        if self.plotting==True:
            plt.figure(figsize=(8, 1))
            plt.title('Finished Spread')
            plt.plot((range(len(self.V))), self.V, color='k')
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
        
        tenth = MCS//10
        tenths = []
        for i in range(9):
                tenths.append(tenth*(i+1))
        j = 1
        for i in range(MCS):
            N = len(V)
            ridx = np.random.randint(N-2)
            if ridx != 0:
                if V[ridx] == V[ridx+1]:
                    V[ridx-1] = V[ridx]; V[ridx+2] = V[ridx]
                elif V[ridx] == -V[ridx+1]:
                    V[ridx-1] = V[ridx+1]; V[ridx+2] = V[ridx]
                
                self.dT.append(i - self.times[ridx+2])
                self.dT.append(i - self.times[ridx-1])
                
                self.times[ridx+2] = i
                self.times[ridx-1] = i
            
            elif ridx == 0:
                if V[ridx] == V[ridx+1]:
                    V[ridx-1] = V[ridx]
                elif V[ridx] == -V[ridx+1]:
                    V[ridx-1] = V[ridx+1]
                    
                
                self.dT.append(i - self.times[ridx+2])
                self.times[ridx+2] = i
                
            m.append(self.magnet())
            
            if self.plotting==True and i in tenths:
                plt.figure(figsize=(8, 1))
                plt.title('Midway spread (%i%%)' %(j*10))
                j += 1
                
                plt.plot((range(len(V))), V, color='k')
                plt.show()
                
            self.time += 1
            
        
        return V, m
    
    
    def magnet(self):
        m = (1/len(self.V)) * np.sum(self.V)
        return m
    
    
class person:
    
    def __init__(self,
             opinion):
        self.opinion = opinion
        self.changed = False
        
    def change(self, s, i, ct):
        if self.changed == False:
            self.opinion = s
            self.changed = True
            self.change_time = i
            
            
        return self.opinion
        
        
    