#Foroutan(2577847)-Swarna(2577834).-Assignment1-Exercise3

#Shift Or Algorithm

import sys

####################################################
def shiftor(text, pattern):
    #creating masks of characters in pattern
    m=len(pattern)
    masks=dict()
    bit=1
    for c in pattern:
        if c not in masks:
            masks[c]=0
        masks[c] |= bit
        bit *=2

    accept_state=bit//2
    accept_state=((1<<m)-1 )-accept_state


    #complement of masks of characters in pattern
    for c in masks:
        masks[c]=((1<<m)-1 )- masks[c]

    for c in text:
        if c not in pattern:
            masks[c]=(1<<m)-1    



    #cnosider 1 for each bit of active state D variable 
    D=0
    D=((1<<m)-1 )- D  #complement D=00000 to D=11111

    #If a state is active, the corrospending bit is zero
    i=0
    flag=0
    for c in text:
        #print("i=", i,"C=",c)
        #shift - Or
        temp=bin(D<<1)   #shift
        if(len(temp)>=(m+3)):
           temp=temp[3:]    #omit extra bit
        D=int(temp,2)
        D=D | masks[c]    #Or

        #print(bin(D), len(bin(D)))
        #If the last digit is zero, it is not considered in the digits of D,since it is after valuble digits,
        #EXample: D=15=11110 , the saved D will be D=1111
        #since only when we reach the final state, the last digit of D will be zero,
        #so if the length of D be one less than length of pattern (m), it shows the last digit is zero.
        #length of pattern=m=5 , D=1111 , bin(D)=0b1111 , lenth of bin(D) is 6=m+1=4 digits + 2 of (0b)
        #that result in finding the pattern.
        if len(bin(D))<=(m+1): 
            if flag==0:
                flag=1
                print(i-m+1, end=" ")
            else:
                print(",",i-m+1, end=" ")
        i+=1
        
    if flag==0:
        print("The pattern is not in the text file.")


if __name__ == '__main__':
    shiftor(sys.argv[1], open(sys.argv[2]).read().strip())
