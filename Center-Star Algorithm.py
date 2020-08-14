import sys
import collections
from Bio import SeqIO, Seq, SeqRecord
from Bio.SeqRecord import SeqRecord

def sim(c1,c2, match, mismatch):
    if c1==c2:
        return match
    else:
        return mismatch

def rev(s): 
  st = "" 
  for i in s: 
    st = i + st
  return st

def NWscore(x,y, match, mismatch, gap):
    #print("NWscor",x,y)
    #store just two rows as score(0) & score(1)
    score = []
    score.append([None]*(len(y)+1))
    score.append([None]*(len(y)+1))
    
    score[0][0] = 0
    for j in range(1,len(y)+1):
        score[0][j] = score[0][j-1] + gap     #first row    
        
    for i in range(1,len(x)+1):             #for each row
        score[1][0] = score[0][0] + gap
        for j in range(1,len(y)+1):
            if x[i-1]==y[j-1]:
                score[1][j] = score[0][j-1] + match
            else:
                scoreD = score[0][j-1] + mismatch
                scoreT = score[0][j] + gap
                scoreL = score[1][j-1] + gap
                score[1][j] = min(scoreD, scoreT, scoreL)
        for j in range((len(y)+1)):
            score[0][j] = score[1][j]
        
    last_row = score[1]

    #print("last",last_row)
    return last_row
            
        

def NW(x,y, match, mismatch, gap):
    #length of x or y == 1
    #Calculate the first two rows as score(0) & score(1)
    score = []
    for i in range(len(x)+1):
        score.append([None]*(len(y)+1))
    
    score[0][0] = 0
    for j in range(1,len(y)+1):
        score[0][j] = score[0][j-1] + gap     #first row
        
    for i in range(1,len(x)+1):             #for each row
        score[i][0] = score[i-1][0] + gap
        for j in range(1,len(y)+1):
            if x[i-1]==y[j-1]:
                score[i][j] = score[i-1][j-1] + match
            else:
                scoreD = score[i-1][j-1] + mismatch
                scoreT = score[i-1][j] + gap
                scoreL = score[i][j-1] + gap
                score[i][j] = min(scoreD, scoreT, scoreL)
    

    #backtrace
    m = len(x) #side
    n = len(y) #top
    z=""  
    w=""
    while( m>0 or n>0):
        if m>0 and n>0 and score[m][n]== score[m-1][n-1] + sim(x[m-1],y[n-1], match, mismatch):
            z = x[m-1] + z
            w = y[n-1] + w
            m = m - 1
            n = n - 1
        elif m>0 and score[m][n]==score[m-1][n] + gap: # check horizental
            z = x[m-1] + z
            w = '-' + w
            m = m - 1
        else:
            z = '-' + z
            w = y[n-1] + w
            n = n - 1
        
    #print(z,w)
    return z,w


def linear_space(x,y, match, mismatch, gap):
    #print("\nlinear_space",x,y)
    z=""   #Align x
    w=""   #Align y
    if len(x)==0:
        for i in range(len(y)):
            z = z + '-'
            w = w + y[i]
    elif len(y)==0:
        for i in range(len(x)):
            z = z + x[i]
            w = w + '-'
    elif len(x)==1 or len(y)==1:
        z,w = NW(x,y, match, mismatch, gap)
    else:
        xlen=len(x)
        xmid=len(x)//2
        ylen=len(y)
        
        scoreT = NWscore(x[:xmid],y, match, mismatch, gap)
        scoreB = NWscore(rev(x[xmid:]),rev(y), match, mismatch, gap)
        scoreB.reverse()
        sum_sc = list(scoreT)
        for i in range(len(scoreB)):
            sum_sc[i]+=scoreB[i]
            
        ymid = sum_sc.index(min(sum_sc))   #min for match 0 
        
        z1 , w1= linear_space(x[:xmid], y[:ymid], match, mismatch, gap)
        z2 , w2= linear_space(x[xmid:],y[ymid:], match, mismatch, gap)
        z = z1 + z2
        w = w1 + w2
        
    return z,w
        
                 


def input_fun(file_in, file_out):
    match=0
    mismatch=1
    gap=1
    #reading the input file
    seqs=SeqIO.parse(open(file_in, mode='r'), 'fasta')
    stri=[]
    x=""
    i=0
    for se in seqs:
        stri.append(str(se.seq))
        i+=1

    #Center-Star algorithm
    dp=[0]*len(stri)
    for i in range(len(stri)): #for each k seq
        s=stri[i]
        for j in range(len(stri)): #find dp
            if i != j:
                t = stri[j]
                align_score = NWscore(s,t, match, mismatch, gap)
                dp[i] += align_score[len(t)]
                
    print(*dp,sep=",")
    
    #finding the center seq
    center=dp.index(min(dp))
    
    #align center seq with other seqs

#multi_align=[""]*len(stri)
    for j in range(len(stri)):
        if center != j :
            a=""
            b=""
            a,b=linear_space(stri[center],stri[j], match, mismatch, gap)
            stri[j] = b
            if a != stri[center]:
                stri[center] = a
                for k in range(0,j): #realign the previous aligned seqs
                    a,stri[k]=linear_space(stri[center],stri[k], match, mismatch, gap)
                     
        
    #writing in a file
    f= open(file_out,"w+")

    #writing assembled seq in the out put file
    
    for num, s in enumerate(stri):
        f.write(">%d\n" %(num+1))
        f.write("%s\n"%stri[num])
        


if __name__ == '__main__':
    input_fun(sys.argv[1], sys.argv[2])
