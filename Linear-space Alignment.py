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
                score[1][j] = max(scoreD, scoreT, scoreL)
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
                score[i][j] = max(scoreD, scoreT, scoreL)
    

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


def linear_space(x,y, match, mismatch, gap,flag):
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
            
        ymid = sum_sc.index(max(sum_sc))   #min for match 0 
        if(flag==True):
            flag = False
            print("alignment score =",sum_sc[ymid])
        
        z1 , w1= linear_space(x[:xmid], y[:ymid], match, mismatch, gap,flag)
        z2 , w2= linear_space(x[xmid:],y[ymid:], match, mismatch, gap,flag)
        z = z1 + z2
        w = w1 + w2
        
    return z,w
        
                 


def input_fun(file_in, file_out, match, mismatch, gap):
    #print("input:",file_in, file_out, match, mismatch, gap)

    seqs=SeqIO.parse(open(file_in, mode='r'), 'fasta')
    stri=["", ""]
    for i, se in enumerate(seqs):
        stri[i]=str(se.seq)
        #print(i,stri[i])

    s=stri[0]
    t=stri[1]

        
    #output
    '''align_score = NWscore(s,t, match, mismatch, gap)
    print("alignment score=",align_score[len(t)])'''
    
    a=""
    b=""
    flag=True
    a,b=linear_space(s,t, match, mismatch, gap,flag)
    print(a)
    print(b)

    #writing in a file
    f= open(file_out,"w+")

    #writing assembled seq in the out put file
    num=1
    f.write(">%d\n" %num)
    f.write("%s\n"%a)
    num+=1
    f.write(">%d\n" %num)
    f.write("%s\n"%b)

    
    


if __name__ == '__main__':
    input_fun(sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5]))
