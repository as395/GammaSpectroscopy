
# coding: utf-8

# In[1]:


get_ipython().magic('matplotlib inline')
import numpy as np
import matplotlib.pyplot as plt
import scipy as sp


# In[2]:


'''Function to determine full count ratios WITH summing (S)'''
def FullSum(x, F, Energies, et, ep, Q, alpha, debug, show):        
    
    '''Normalizing X and F'''
    # Normalize branching ratios
    for i in range(0,len(x)):
        xsum = x[i].sum()
        if np.equal(xsum,0)==True:
            for j in range(0,len(x)):
                x[i,j] = x[i,j]
                j=j+1
        else:
            if np.not_equal(xsum,1)==True:
                for j in range(0,len(x)):
                    x[i,j] = x[i,j]/float(xsum)   #float(xsum)
                    j=j+1
        i=i+1  
    
    # Normalize feeding fractions
    fsum=F.sum()
    # print(fsum)
    if F.sum()>1:
        for i in range(0,len(F)):
            F[i] = F[i]/float(fsum)
            
    '''Probabilities'''
        # probability of branch decay occuring (accounting for loss through conversion)
    c = x/(1+alpha) 
        # probability of peak detected at expected photopeak (decay * peak efficiency)
    a = c*ep/Q
        # probability of peak detected anywhere in spectrum (decay * total efficiency)
    e = c*et/Q
        # probability of no gamma ray detection (used for summing out)
    b = x-e

    E = np.identity(w)
    
    '''FULL SUMMING'''
    '''These while loops currently compare float values to zero without a defined tolerance for equality.
    I can set the absolute tolerance using np.isclose() or math.isclose(), but I need to decide on a threshold value.'''

        # summing in 
    aa = a.dot(a)
    asum = aa.sum
    asum = aa.sum(dtype=float)
    while (asum > 0):
        a = a + aa
        aa = np.dot(aa,a)
        asum = aa.sum(dtype=float)
    A = Q*a

        # summing out
    bb = b.dot(b)
    bsum = bb.sum(dtype=float)
    while (bsum > 0):
        b = b + bb
        bb = np.dot(bb,b)
        bsum = bb.sum(dtype=int)
    B = E + b #kth power

    '''Calculating full measured peaks (S)'''

    # EQUATION 6

    n = F.dot(B)
    N = E*n


    m = B[0:w,0]
    M = E*m
    
    S = R*N.dot(A).dot(M)
    
    if debug==True:
        print("Energy levels:",'\n',Energies)
        print("Total efficiency:",'\n',et)
        print('Tline', Tline)
        print("Peak efficiencies:",'\n',ep) 
        print("Conversion coefficients:",'\n',alpha)
        print("alpha caps", astart, aend)
        print("Source disintegration rate:", '\n', R)
        print("Number of detectors:", '\n', Q)
        print("A", '\n', A)
        print("B", '\n', B)
    
    if show==True:
        print("Branching fractions:",'\n',x)
        print("Feeding fractions:",'\n',F)
        print("Summed peaks", '\n', S)
    
    return(S)


# In[3]:


'''Function to determine count ratios WITHOUT summing (So)'''
def NoSum(x, F, Energies, et, ep, Q, alpha, debug, show):        
    
    '''Normalizing X and F'''
    # Normalize branching ratios
    for i in range(0,len(x)):
        xsum = x[i].sum()
        if np.equal(xsum,0)==True:
            for j in range(0,len(x)):
                x[i,j] = x[i,j]
                j=j+1
        else:
            if np.not_equal(xsum,1)==True:
                for j in range(0,len(x)):
                    x[i,j] = x[i,j]/float(xsum)   #float(xsum)
                    j=j+1
        i=i+1  
    
    # Normalize feeding fractions
    fsum=F.sum()
    # print(fsum)
    if F.sum()>1:
        for i in range(0,len(F)):
            F[i] = F[i]/float(fsum)
            
    '''Probabilities'''
        # probability of branch decay occuring (accounting for loss through conversion)
    c = x/(1+alpha) 
        # probability of peak detected at expected photopeak (decay * peak efficiency)
    a = c*ep/Q
        # probability of peak detected anywhere in spectrum (decay * total efficiency)
    e = c*et/Q
        # probability of no gamma ray detection (used for summing out)
    b = x-e

    E = np.identity(w)

    '''FULL SUMMING'''
    '''These while loops currently compare float values to zero without a defined tolerance for equality.
    I can set the absolute tolerance using np.isclose() or math.isclose(), but I need to decide on a threshold value.'''

        # summing in 
    aa = a.dot(a)
    asum = aa.sum
    asum = aa.sum(dtype=float)
    while (asum > 0):
        a = a + aa
        aa = np.dot(aa,a)
        asum = aa.sum(dtype=float)
    A = Q*a

        # summing out
    bb = b.dot(b)
    bsum = bb.sum(dtype=float)
    while (bsum > 0):
        b = b + bb
        bb = np.dot(bb,b)
        bsum = bb.sum(dtype=int)
    B = E + b #kth power
    
    '''NO SUMMING'''
    Ao = a
    Bo = E+x

    no = F.dot(Bo)
    No = E*no

    mo = Bo[0:w,0]
    Mo = E

    So = R*No.dot(Ao).dot(Mo)
    
    if debug==True:
        print("Energy levels:",'\n',Energies)
        print("Feeding fractions:",'\n',F)
        print("Branching fractions:",'\n',x)
        print("Total efficiency:",'\n',et)
        print('Tline', Tline)
        print("Peak efficiencies:",'\n',ep) 
        print("Conversion coefficients:",'\n',alpha)
        print("alpha caps", astart, aend)
        print("Source disintegration rate:", '\n', R)
        print("Number of detectors:", '\n', Q)
        print("Ao", '\n', Ao)
        print("Bo",'\n', Bo)
        print("(F . B_o)", '\n', no)
        print("No",'\n', No)
        print("Mo",'\n', Mo)
    
    if show==True:
        print("Peaks without summing:",'\n',So)
    
    return(So)


# In[4]:


'''READ IN SOURCE FILE (Forwards)'''
if __name__ == '__main__':
    '''INCLUDE TEXT FILE'''
    Data = np.genfromtxt("Co-60", comments="#", delimiter=",", filling_values='0')
    Efficiencies = np.genfromtxt("efficiency.txt", comments="#", delimiter=",", filling_values='0')

    # Order of matrix (number of energy levels, including the ground state)
    w = len(Data[0])
    '''Defining energy drops (keV)
    This section reads in energy levels and converts them to a matrix containing the energies of gamma rays emitted from all possible energy drops.
    ''' 
    # Creates empty matrix of appropriate dimensions, and reads in data from first line of text file
    Energies = np.zeros((w,w))
    length_end = w+1
    # print(length_end)
    for i in range(0,w):
    # creates empty vector
        v=0
        for j in range(0,len(Data[0])):
    #         Defines energy drops
            k = Data[0,i]-Data[0,j]
            j=j+1
            v=np.append(v,k)
    #     Ensures that matrix only contains one copy of each value
        for d in range(0,len(v)):
            if v[d] < 0:
                v[d]=0
    #     Adds each vector to empty matrix
        Energies[i]=np.array(v[1:length_end])
    # Final energy levels (keV)

    '''Feeding Fractions'''
    # Create vector from file
    F = np.array([Data[1]])

    '''Branching Fractions (NORMALIZE)'''
    # Create matrix from file
    x = np.array(Data[2:2+w])  
    
    '''Total Efficiency'''
    Tline = 2+w
    # If debug=True, print Tline
    T = Data[Tline]
    et = T.sum() #for total array

    '''Peak Efficiencies'''
    # Round energy levels to nearest integer
    ep = np.rint(Energies)
    # print(P)
    for i in range(0,len(ep[0])):
        for j in range(0,len(ep[0])):
            if np.not_equal(ep[i,j],0)==True:
                d=ep[i,j]
    #             Change D type to integer and value to reflect line number in TRIUMF data 
                D=d.astype(int)-9
    #             Change array value to peak efficiency for specified energy
                ep[i,j] = Efficiencies[D,1]

    '''Conversion Coefficients'''
    astart = Tline+1
    aend = astart + w
    # if debug=True, print astart and aend ("alpha caps")
    alpha = np.array(Data[astart:aend])
    ###alpha = np.ones((4,4))*0 #no conversion coefficients for toy system

    '''Source disintegration rate'''
    rline = aend
    r = np.array(Data[rline])
    R = r.sum()

    '''Number of detectors'''
    nline = rline + 1
    q = (Data[nline])
    Q=1


# In[5]:


'''USE FUNCTIONS'''    
S=FullSum(x, F, Energies, et, ep, Q, alpha, debug=False, show=True)
So=NoSum(x, F, Energies, et, ep, Q, alpha, debug=False, show=True)


# In[6]:


'''READ IN PEAK FILES'''
'''Inputs: Measured peaks (S), peak efficiencies (ep), total efficiency (et)'''
PeakData = np.genfromtxt("Co-60 PEAKS.txt", comments="#", delimiter=",", filling_values='0')
Efficiencies = np.genfromtxt("efficiency.txt", comments="#", delimiter=",", filling_values='0')

def PeaktoSource(PeakData, Efficiencies, debug, show):
    # Order of matrices (number of energy levels, including the ground state)
    w = len(PeakData[0])

    '''PROVIDED VALUES'''
    '''Source disintegration rate'''
    rline = 0
    r = np.array(PeakData[rline])
    R = r.sum()
    # R = 100

    '''Total Efficiency'''
    Tline = rline+1
    # print('Tline', Tline)
    T = PeakData[Tline]
    et = T.sum() #for total array
    # et = 1

    '''Conversion Coefficients'''
    astart = Tline+1
    aend = astart + w
    # print("Acaps", astart, aend)
    alpha = np.array(PeakData[astart:aend])
    ###alpha = np.ones((4,4))*0 #no conversion coefficients for toy system

    '''Number of detectors'''
    nline = aend
    q = (PeakData[nline])
    Q = q.sum()

    '''Energy Levels (keV)'''
    # Creates empty matrix of appropriate dimensions, and reads in data from line of text file
    Energies = np.zeros((w,w))
    length_end = w+1
    eline = nline+1
    # print(length_end)
    for i in range(0,w):
    # creates empty vector
        v=0
        for j in range(0,w):
    #         Defines energy drops
            k = PeakData[eline,i]-PeakData[eline,j]
            j=j+1
            v=np.append(v,k)
    #     Ensures that matrix only contains one copy of each value
        for d in range(0,len(v)):
            if v[d] < 0:
                v[d]=0
    #     Adds each vector to empty matrix
        Energies[i]=np.array(v[1:length_end])
    # Final energy levels (keV)

    '''Measured Peaks'''
    MPstart = eline+1
    MPend = MPstart + w
    # print("MPcaps", MPstart, MPend)
    S = np.array(PeakData[MPstart:MPend])

    '''Peak Efficiencies'''
    # Round energy levels to nearest integer
    ep = np.rint(Energies)
    # print(P)
    for i in range(0,len(ep[0])):
        for j in range(0,len(ep[0])):
            if np.not_equal(ep[i,j],0)==True:
                d=ep[i,j]
    #             Change D type to integer and value to reflect line number in TRIUMF data 
                D=d.astype(int)-9
    #             Change array value to peak efficiency for specified energy
                ep[i,j] = Efficiencies[D,1]
    
    if debug==True:
        print('Number of Detectors', '\n', Q)
        print("Conversion coefficients:",'\n',alpha)
        print("Source disintegration rate:", '\n', R)
        print("Total efficiency:",'\n',et)
    
    if show==True:
        print("Energy levels:",'\n', Energies)
        print("Peak efficiencies:",'\n',ep)
        print("Measured Peaks:",'\n',S)        
    
    return(S)


# In[7]:


Sm=PeaktoSource(PeakData, Efficiencies, debug=False, show=True)


# In[8]:


'''Guess branching ratios (Xg)'''
def findX(w, S, ep, bebug, show):
    Xg = np.zeros((w,w))
    for i in range(0,len(ep[0])):
        for j in range(0,len(ep[0])):
            if np.not_equal(S[j,i],0)==True:
                Xg[j,i]=S[j,i]/ep[j,i]
# Attempt to remove summing out:
#         for j in range(1,len(ep[0])):
#             if np.not_equal(S[j,i],0)==True:
#                 Xg[j,i]=S[j,i]/ep[j,i]
#         for j in range(0,0):
#             if np.not_equal(S[0,i],0)==True:
#                 Xg[0,i]=(S[0,i]/ep[0,i])-sum(S[0:len(S[0]),i]/ep[0:len(ep[0]),i])
                
    if bebug==True:
        print("Peaks",'\n', S)
        print("peak efficiencies",'\n', ep)
        
    if show==True:
        print("Branching ratio guess:",'\n', Xg)
    
    return(Xg)

'''Guess feeding fractions (Fg)'''
def findF(w, S, ep, bebug, show):
    Xg = np.zeros((w,w))
    for i in range(0,len(ep[0])):
        for j in range(0,len(ep[0])):
            if np.not_equal(S[i,j],0)==True:
                Xg[i,j]=S[i,j]/ep[i,j]
     
    Fg = np.zeros(w)
    for i in range(0,len(Xg[0])):
        for j in range(0,len(Xg[0])):
            fg = sum(Xg[i,0:len(Xg[0])])-sum(Xg[0:len(Xg[0]),i])
        Fg[i]=fg
        for d in range(0,len(Fg)):
            if Fg[d] < 0:
                Fg[d]=0
            
    if bebug==True:
        print("Peaks",'\n', S)
        print("peak efficiencies",'\n', ep)
        print("Branching ratio guess:",'\n', Xg)
        
    if show==True:
        print("Feeding ratio guess:",'\n', Fg)
    
    return(Fg)

'''Guess ONLY F (Fgonly)'''
def findonlyF(Xg, bebug, show):
    Fgonly = np.zeros(w)
    for i in range(0,len(Xg[0])):
        for j in range(0,len(Xg[0])):
            fg = sum(Xg[i,0:len(Xg[0])])-sum(Xg[0:len(Xg[0]),i])
        Fgonly[i]=fg
        for d in range(0,len(Fgonly)):
            if Fgonly[d] < 0:
                Fgonly[d]=0
            
    if bebug==True:
        print("Branching ratio guess:",'\n', Xg)
        
    if show==True:
        print("Feeding ratio guess:",'\n', Fg)
    
    return(Fgonly)


# In[9]:


Xg=findX(w, Sm, ep, bebug=False, show=True)
Fg=findF(w, Sm, ep, bebug=False, show=True)
Fgonly=findonlyF(Xg, bebug=False, show=True)


# In[10]:


Sg=FullSum(Xg, Fg, Energies, et, ep, Q, alpha, debug=False, show=True)
# print("Actual Summed Peaks",'\n', S)
Sdiff=Sg-S
# print("Difference", '\n', Sdiff)
SR0=Sdiff/S
# print(SR0)


# In[11]:


'''PLOT RESULTS OF LOOP'''
#  The goal here is to adjust Xg values to make Sg --> S.
#  Unfortunately, Xg impacts many values of Sg.

def plotSconvergence(p, delta, show):
    Sratio = np.zeros((w,w))
    '''Generalize this process:'''
#      Set up empty arrays for each desired matrix element to plot evolution of Sg:
    y10=[]
    y20=[]
    y21=[]
    y30=[]
    y31=[]
    y32=[]
    z10=0
    z20=0
    z21=0
    z30=0
    z31=0
    z32=0
    for x in range(0, p):        
        Fgonly=findonlyF(Xg, bebug=False, show=False)
        Sg=FullSum(Xg, Fg, Energies, et, ep, Q, alpha, debug=False, show=False)
        Sdiff=Sg-S
        for i in range(0,len(ep[0])):
            for j in range(0,len(ep[0])):
                if np.not_equal(S[i,j],0)==True:
                    Sratio[i,j]=Sdiff[i,j]/S[i,j]
        z10=Sdiff[1,0]
        z20=Sdiff[2,0]
        z21=Sdiff[2,1]
        z30=Sdiff[3,0]
        z31=Sdiff[3,1]
        z32=Sdiff[3,2]
        y10.append(z10)
        y20.append(z20)
        y21.append(z21)
        y30.append(z30)
        y31.append(z31)
        y32.append(z32)
        for i in range(0,len(Xg[0])):
            for j in range(0,len(Xg[0])):
                if (np.absolute(Sdiff[i,j])>delta):
                    Xg[i,j]=Xg[i,j]-Xg[i,j]*Sratio[i,j]
        
    if show==True:
        print("Final Guess (Sg)", '\n', Sg)
        print("Goal (S)" '\n', S)
        print("New Ratio (Sg/S)", '\n', Sratio)
#         print('Branching Fractions', '\n', Xg)
        
    return np.array(y10), np.array(y20), np.array(y21), np.array(y30), np.array(y31), np.array(y32)

# Recover arrays for each Sg[i,j] and plot:
Y10, Y20, Y21, Y30, Y31, Y32 = plotSconvergence(10000, 0.001, show=False)
plt.plot(Y10, 'violet')
plt.plot(Y20, 'b:')
plt.plot(Y21, 'g--')
plt.plot(Y30, 'r')
plt.plot(Y31, 'orange')
plt.plot(Y32, 'yellow')
