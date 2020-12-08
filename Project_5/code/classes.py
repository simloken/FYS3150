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
            plotting=False,
            switch_proba=None):
        
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
        self.switch_proba = switch_proba
        
        
    def run(self):
        Alst = []; Blst = []; ABlst = [];
        if self.clusterd == False:
            self.V = np.random.choice([-1,1], size=(self.N)) #either up or down
        else:
            self.V = np.ones(self.N) #array of all A
            if self.random == True: #random placements of B, total N*cB
                k = 0
                while k < int(self.N*self.cB):
                    ridx = np.random.randint(self.N)
                    if self.V[ridx] != -1:
                        self.V[ridx] = -1
                        k += 1
            else: #first cB*N elements are B
                for i in range(int(self.N*self.cB)):
                    self.V[i] = -1
                
            
        if self.plotting==True: #initial distribution
            plt.figure(figsize=(8, 1))
            plt.ylim([-1.05,1.05])
            plt.title('Initial Spread')
            plt.plot((range(len(self.V))), self.V, color='k')
            plt.xlabel('Element i')
            plt.ylabel('Opinion')
            plt.show()
        
        
        self.V, self.m = self.dynamicRules(self.V, self.MCS) #run over with MCS montecarlo steps
       
        
        if self.plotting==True: #final distribution
            plt.figure(figsize=(8, 1))
            plt.ylim([-1.05,1.05])
            plt.title('Finished Spread')
            plt.plot((range(len(self.V))), self.V, color='k')
            plt.xlabel('Element i')
            plt.ylabel('Opinion')
            plt.show()
            
        A = 0; B = 0; AB = 0
        for i in range(self.NN): #pick random spot NN times to gauge spin consensus
            ridx = np.random.randint(self.N-2)
            lst = []
            lst.append(self.V[ridx-1]) #add closest neighbours (2)
            lst.append(self.V[ridx-2]) #on both sides
            lst.append(self.V[ridx+1]) #then gauge
            lst.append(self.V[ridx+2]) #contents of list
            if 1 not in lst: #if all -1
                B += 1
            elif -1 not in lst: #if all +1
                A += 1
            elif 1 in lst and -1 in lst: #if both
                AB += 1
    
        Alst.append(A); Blst.append(B); ABlst.append(AB) #add states to list, remnant from when I would do multiple runs
            
        A = np.mean(Alst) #find means
        B = np.mean(Blst)
        AB = np.mean(ABlst)
        
        A /= self.NN #percentages
        B /= self.NN
        AB /= self.NN
        
        print("State 1: {:.2f} (AAAA)\nState 2: {:.2f} (BBBB)\nState 3: {:.2f} (ABAB)".format(A, B, AB))
        
        return A, B, AB, self.m
        
        
        
        
        
        
        
        
    def dynamicRules(self, V, MCS):
        m = []
        
        tenth = MCS//10
        tenths = []
        for i in range(9): #for plotting midway
                tenths.append(tenth*(i+1))
        j = 1
        N = len(V)
        has_happened = False
        for i in range(MCS):
            ridx = np.random.randint(N-1) #random element chosen
            if self.switch_proba == None and has_happened == False: #stops sim early, only active when chance 
                if np.equal(V, -1).all() == True: #to switch if off (naturally)
                    print('Consensus Reached: No') 
                    has_happened=True #activates, extra MCS*0.05 steps before stopping sim (looks better for m, etc.)
                    self.consensus = 'No'
                    timer = 0
                elif np.equal(V, 1).all() == True:
                    print('Consensus Reached: Yes')
                    has_happened=True
                    self.consensus = 'Yes'
                    timer = 0
                elif np.equal(V[::2], -1).all() == True or np.equal(V[1::2], -1).all() == True:
                    if np.equal(V, 1).all() == False and np.equal(V, -1).all() == False:
                        print('Consensus Reached: Stalemate')
                        has_happened=True
                        self.consensus = 'Stalemate'
                        timer = 0
                
            if self.switch_proba != None: #whether to have some random chance for 
                chance = np.random.uniform() #element ridx to flip
                if chance <= self.switch_proba:
                        V[ridx] *= -1
                    
            if ridx not in [0, N-2]: #boundaries, must be treated as special cases
                if V[ridx] == V[ridx+1]: #else errors arise
                    V[ridx-1] = self.opinion(V[ridx-1], V[ridx], i, ridx, -1)
                    V[ridx+2] = self.opinion(V[ridx+2], V[ridx], i, ridx, +2)
                elif V[ridx] == -V[ridx+1]:
                    V[ridx-1] = self.opinion(V[ridx-1], V[ridx+1], i, ridx, -1)
                    V[ridx+2] = self.opinion(V[ridx+2], V[ridx], i, ridx, +2)
                
            
                
            
            elif ridx == N-2: #upper bound special case
                if V[ridx] == V[ridx+1]:
                    V[ridx-1] = self.opinion(V[ridx-1], V[ridx], i, ridx, -1)
                elif V[ridx] == -V[ridx+1]:
                    V[ridx-1] = self.opinion(V[ridx-1], V[ridx+1], i, ridx, -1)
            
            elif ridx == 0: #lower bound special case
                if V[ridx] == V[ridx+1]:
                    V[ridx+2] = self.opinion(V[ridx+2], V[ridx], i, ridx, +2)
                elif V[ridx] == -V[ridx+1]:
                    V[ridx+2] = self.opinion(V[ridx+2], V[ridx], i, ridx, +2)
                    
                
            m.append(self.magnet(V)) #find m for current MCS
            
            if has_happened==True:
                timer += 1
                
                if timer >= np.ceil(MCS*0.05):
                    break
            
            if self.plotting==True and i in tenths: #print distribution for each 10th of MCS
                plt.figure(figsize=(8, 1))
                plt.title('Midway spread (%i%%)' %(j*10))
                plt.xlabel('Element i')
                plt.ylabel('Opinion')
                plt.ylim([-1.05,1.05])
                j += 1
                
                plt.plot((range(len(V))), V, color='k')
                plt.show()
                
            self.time += 1
            
        
        return V, m
    
    
    def magnet(self, V):
        m = (1/len(V)) * np.sum(V)
        return m
    
    def opinion(self, reciever, flipper, i, ridx, tun):
        
        if reciever != flipper: #if opinions are different
            reciever = flipper #change opinion of reciever
            self.dT.append(i-self.times[ridx+tun]) #causes the weird plot spread in c (why?!)
            
            self.times[ridx+tun] = i #track time for when opinion was changed
        return reciever #returns (potentially new) opinion
