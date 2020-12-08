import numpy as np
import matplotlib.pyplot as plt
class Simulation:
    
    def __init__(self,
            N,
            NN,
            MCS,
            cB=0.5,
            clusterd=False,
            random = False,
            plotting=False):
        
        self.N = N
        self.NN = NN
        self.MCS = MCS
        self.plotting = plotting
        self.times = np.zeros(N)
        self.dT = []
        self.time = 0
        self.clusterd = clusterd
        self.cB = cB
        self.random = random
        self.times = np.zeros(1000)
        
        
    def run(self):
        Alst = []; Blst = []; ABlst = [];
        if self.clusterd == False:
            self.V = np.random.choice([-1,1], size=(self.N)) #either up or down
        else:
            self.V = np.ones(self.N) #array of all A
            if self.random == True:
                k = 0
                while k < int(self.N*self.cB):
                    ridx = np.random.randint(self.N)
                    if self.V[ridx] != -1:
                        self.V[ridx] = -1
                        k += 1
            else:
                for i in range(int(self.N*self.cB)):
                    self.V[i] = -1
                
            
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
            ridx = np.random.randint(N-1)
            if ridx not in [0, N-2]:
                if V[ridx] == V[ridx+1]:
                    V[ridx-1] = self.opinion(V[ridx-1], V[ridx], i, ridx, -1)
                    V[ridx+2] = self.opinion(V[ridx+2], V[ridx], i, ridx, +2)
                elif V[ridx] == -V[ridx+1]:
                    V[ridx-1] = self.opinion(V[ridx-1], V[ridx+1], i, ridx, -1)
                    V[ridx+2] = self.opinion(V[ridx+2], V[ridx], i, ridx, +2)
                
            
                
            
            elif ridx == N-2: #upper bound
                if V[ridx] == V[ridx+1]:
                    V[ridx-1] = self.opinion(V[ridx-1], V[ridx], i, ridx, -1)
                elif V[ridx] == -V[ridx+1]:
                    V[ridx-1] = self.opinion(V[ridx-1], V[ridx+1], i, ridx, -1)
            
            elif ridx == 0: #lower bound
                if V[ridx] == V[ridx+1]:
                    V[ridx+2] = self.opinion(V[ridx+2], V[ridx], i, ridx, +2)
                elif V[ridx] == -V[ridx+1]:
                    V[ridx+2] = self.opinion(V[ridx+2], V[ridx], i, ridx, +2)
                    
                
            m.append(self.magnet(V))
            
            if self.plotting==True and i in tenths:
                plt.figure(figsize=(8, 1))
                plt.title('Midway spread (%i%%)' %(j*10))
                j += 1
                
                plt.plot((range(len(V))), V, color='k')
                plt.show()
                
            self.time += 1
            
        
        return V, m
    
    
    def magnet(self, V):
        m = (1/len(V)) * np.sum(V)
        return m
    
    def opinion(self, reciever, flipper, i, ridx, tun):
        
        if reciever != flipper:
            reciever = flipper
            self.dT.append(i-self.times[ridx+tun]) #causes the weird plot spread in c
            #self.dT.append(i)
            
            self.times[ridx+tun] = i
            
        
        return reciever
