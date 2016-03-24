import numpy as np
import math as math
import os
import sys

N=int(sys.argv[1]) # number of people with disease in database
mg=int(sys.argv[2])  # how many people with disease have mutations in this gene (usually only 2, but could be 3)
gname=sys.argv[3]  #Gene name
recessive=int(sys.argv[4]) # 0 for Dominant, 1 for compound het, 2 for recessive.
pmiss=float(sys.argv[5]) # probability of person with disease not having a mutation in any causal gene
inv=sys.argv[6]  #  list of Inbreeding coefficients- 0 if absent- Per person with disease in database. Can input 0 for whole vector

InbreedingF=np.zeros(N)
if len(inv)>5:
    arr=inv.split(",")
    if not (len(arr)==N):
        print "\n Number of members in inbreeding list should be N !! \n "
        assert len(arr)==N
    InbreedingF[0]=float(arr[0][1:])
    InbreedingF[-1]=float(arr[-1][:-1])
    for i in range(1,N-1):
        InbreedingF[i]=float(arr[i])


os.chdir(r'/Users/sebastianakle/Documents/Matchmaker')
filename='global_prior_genes_disease_and_non_disease_v3.csv'
g=open(filename,'r')
w=g.readlines()[0]
g.close()
genes=[]
q=w.split("\r")
ll=len(q)-1
psqND=np.zeros(ll)
pND=np.zeros(ll)
psqD=np.zeros(ll)
pD=np.zeros(ll)
c=0
for ver in q[1:]:
    arr=ver.split(",")
    genes.append(arr[0])
    psqND[c]=arr[2]
    pND[c]=arr[1]
    psqD[c]=arr[4]
    pD[c]=arr[3]
    c+=1
pNDord=np.sort(pND)[::-1]
pDn0=pD[pD!=0.0]

order=pND.argsort()[::-1]
psqNDord=psqND[order]
pDord=pD[order]
geneDict={}
revgeneDict={}
for i in range(len(pND)):
    geneDict[genes[order[i]]]=i
    revgeneDict[i]=genes[order[i]]

def choose(N,k):
    p=1
    for j in range(N-k+1,N+1):
        p=p*j
    return (p)/math.factorial(k)

def bino(k,N,p):
    return choose(N,k)*pow(p,k)*pow((1-p),N-k)

def realpval(N,mg,fr):
    pj=0
    for j in range(mg):
        #pj+=sta.binom.pmf(j,N,fr)
        pj+=bino(j,N,fr)
    return 1-pj

pDis=pDord*pDord  # This would be the frequencies of disease mutations- !!!This will be tricky to pick!!!!! many issues with just using pD
# make sure the vector here is prdered according to pND values- for ranking
sp=sum(pDis)
mp=min(pDis)
#print sum(np.sort(pDis)[-50:])/sp  # top 50 genes- 80% of density
cutoff=np.sort(pDis)[-50]
#pick 10% of genes with 0 variants  and top 50 genes to have mp (minimum density)
pnew=np.zeros(len(pDis))
for i in range(len(pDis)):
    if(pDis[i]==0):
        if not np.random.randint(10):
            pnew[i]=mp
    else:
        if(pDis[i]>cutoff):
            pnew[i]=mp
        else:
            pnew[i]=pDis[i]

pDis=pnew
sp=sum(pDis)
pDisdensity=pDis/sp # normalized- to refelct probabilities
pDisCum=np.zeros(len(pDis))
cump=0
for i in range(len(pDis)):
    pDisCum[i]=cump+pDisdensity[i]
    cump=pDisCum[i]
    
def MatchesSimAll(N,G,pDis,pDisCum,pfi,reps):
    data=[]
    maxmatch=0
    for rep in range(reps):
        data.append([])
        idg=[]   # Picks G genes uniformly accros pDis density- using cumulative probability vector pDisCum
        for g in range(G):
            idg.append( np.searchsorted(pDisCum,np.random.uniform()) )
        idg=np.array(idg)  
        fg=pDis[idg]  
        sfg=fg.sum()
        pg=fg*(1-pfi)/sfg
        pg=np.append(pg,pfi)
        patients=np.random.multinomial(N,pg)
        success=0
        for g in range(G):
            if (patients[g]>1): #match 
                if(patients[g]>maxmatch):
                    maxmatch=patients[g]
                data[rep].append((patients[g],idg[g],fg[g])) 
    return maxmatch, data

def MatchesSimPower(N,G,k,pDis,pDisCum,pfi,reps):
    data=[]
    c=0
    for rep in range(reps):
        data.append([])
        idg=[]   # Picks G genes uniformly accros pDis density- using cumulative probability vector pDisCum
        for g in range(G):
            idg.append( np.searchsorted(pDisCum,np.random.uniform()) )
        idg=np.array(idg)  
        fg=pDis[idg]  
        sfg=fg.sum()
        pg=fg*(1-pfi)/sfg
        pg=np.append(pg,pfi)
        patients=np.random.multinomial(N,pg)
        success=0
        for g in range(G):
            if (patients[g]>1): #match 
                success+=1
        if (success>=k):
            c+=1
    return (c+0.0)/reps

def Bonferpower(N,G,alpha,fk,pDis,pDisCum,pfi,reps):  # Bonferroni - use bonferroni corrected alpha
    maxmatch,data=MatchesSimAll(N,G,pDis,pDisCum,pfi,reps)
    outp=np.arange(reps)
    tots=np.arange(reps)
    for rep in range(reps):
        matches=0
        rej=0
        for hits in data[rep]:
            pv=realpval(N,hits[0],fk[hits[1]])
            if (pv<alpha):
                matches+=1
        outp[rep]=matches
        tots[rep]=len(data[rep])
    return outp,tots

def GetN(G,k,alpha,fk,pDis,pDisCum,pmiss,power,reps):
    N=2
    alphaB=alpha/len(fk)
    while(1):
        ms=MatchesSimPower(N,G,k,pDis,pDisCum,pmiss,reps)
        if(ms>power):
            break
        else:
            N+=1
        if N==2000:
            print "N>2000"
            return np.nan
    #print N
    while(1):
        ms=Bonferpower(N,G,alphaB,fk,pDis,pDisCum,pmiss,reps)
        match=len(ms[0][ms[0]>=k])+0.0
        if((match/reps)>power):
            break
        else:
            N+=1
        if N==2000:
            print "N>2000"
            return np.nan
    return N+1

T=len(pND)  #20000  # Total number of possible associated genes
pk=pNDord  #np.array([(1.0/10)*np.exp(-0.005*x) for x in range(T)])
# This will be empirical- array of frequencies of alleles with rare mutations per gene. 
#These will be ranked pk[0] the most frequent- probably Titan
priF=0.0005 # avergae human inbreeding coefficient- use as prior
InbreedingF=[i if not i==0 else priF for i in InbreedingF]

if (recessive):
    if(recessive==1):   # Compound hets
        fk=pk*pk  # for each gene, expected freqeuncy of compund hets w rare mutations in our sample
    else:   #Recessive
        Fm=sum(InbreedingF)/N
        fk=Fm*(pk-psqNDord)+psqNDord 
else:  # Dominant
    fk=2*pk   # for each gene, expected freqeuncy of rare mutations in our (diploid) sample
    
# In the recesive case, this will not be very precise, it will assume that every individual has the same 
#probability of being a homozygote for a particular gene- so that the distribution of mutations will be binomial


matchid=geneDict[gname] # id of gene that matched- ID is by rank of frequency pk below.  If it is simpler get rid of this, use matchfreq below instead

pv=realpval(N,mg,fk[geneDict[gname]])

print "Bonferroni corrected p value"
print 1.0-pow((1-pv),T)
