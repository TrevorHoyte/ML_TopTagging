import tensorflow as tf
import IPython.display as display
from PIL import Image
import numpy as np
import matplotlib.pyplot as plt
import os
import pathlib
from tensorflow.keras import layers
import csv
import math


#errors
e4to3=0
e3to4=0

file_read1="/home/thoyte/dataset/ml/vect3.csv"
file_read2="/home/thoyte/dataset/ml/vect3_test.csv"

file_create1='train_preprocessed.csv'
file_create2='test_preprocessed.csv'

fieldnames1 = ["qcd=0, ttbar=1"]
for i in range(1,601):
        fieldnames1.append(i)

def takeFirst(elem):
    return elem[0]




#conversions between coordinates ref https://arxiv.org/pdf/1709.04533.pdf
def convert4to3(vector):
    E=vector[0]
    px=vector[1]
    py=vector[2]
    pz=vector[3]

    pt= math.sqrt(px*px+py*py)
    phi=math.atan2(px,py)
    try:
        #theta=math.acos(pz/E)
        #eta=math.log(1/(math.tan(theta/2)))
        eta=0.5*math.log((E+pz)/(E-pz))
    except:
        print("skip this file 4 to 3 error")
        print(E," ",px," ",py," ",pz)
        global e4to3
        e4to3+=1
        global writetofile
        writetofile=False
        pt=0
        phi=0
        eta=0
    return (pt,eta,phi)

def convert3to4(vector):
    pt=vector[0]
    eta=vector[1]
    phi=vector[2]
    try:
        theta=math.atan(1/math.exp(eta))*2
        #print("theta works")
        pz=pt/math.tan(theta)
        #print("pz works")
        py=math.sqrt(pt**2/(math.tan(phi)**2+1))
        #print("py works, py: ",py,  "   pt: ",pt)
        px=math.sqrt(pt**2-py**2)
        #print("px works")
        E=math.sqrt(px**2+py**2+pz**2)
        #print("E works")
        #check for correct py,px
        if(math.atan2(px,py)==phi):
            
            return (E,px,py,pz)
        elif ((math.atan2(-px,py)==phi)):
            
            return (E,-px,py,pz)
        elif ((math.atan2(px,-py)==phi)):
            
            return (E,px,-py,pz)
            
        else: # ((math.atan2(-px,-py)==phi)):
            return (E,-px,-py,pz)
    except:
        print("skip this file 3 to 4 error")
        print(pt," ", eta, " ", phi)
        global e3to4
        e3to4+=1
        global writetofile1
        writetofile1=False
        return (0,0,0,0)


                           
def manipulate_vectors(file_create,file_read):
    with open(file_create, "w") as hp:
        wr = csv.writer(hp, dialect='excel')
        wr.writerow(fieldnames1)
        
        with open( file_read, 'r') as file_name:
            csv_reader = csv.reader(file_name, delimiter=',')
            linecount=0
            
            for row in csv_reader:
                list_3vect=[]
                writetofile=True
                writetofile1=True
                if linecount==0:
                    linecount+=1
                else:
                    linecount+=1
                    for i in range(1,len(row),3):
                        vect=(float(row[i]),float(row[i+1]),float(row[i+2]))
                        list_3vect.append(vect)
        
                    list_3vect.sort(reverse=True,key=takeFirst) 
                    random.shuffle(list_3vect)
                    # manipulate vectors
                    #pre processing
###move this section around#########################################
                    finallist=[row[0]]
                    #only use 30 constituents
                    for k in range(0,len(list_3vect)):
                        finallist.append(list_3vect[k][0])
                        finallist.append(list_3vect[k][1])
                        finallist.append(list_3vect[k][2])
                    wr.writerow(finallist)                

manipulate_vectors(file_create1,file_read1)
manipulate_vectors(file_create2,file_read2)
##move this section around#############################################
                    
                    #shift in eta and phi
                    for k in range(0,len(list_3vect)):
                                if list_3vect[k]!=(0,0,0):
                                    vector=(list_3vect[k][0],list_3vect[k][1]-list_3vect[0][1],list_3vect[k][2]-list_3vect[0][2])
                                    list_3vect[k]=vector

                    
                                        
                    ## rotate
                    #create 4 vectors skip the hardest one
                    vect4_list=[] 
                    for k in range(1,len(list_3vect)):
                        if list_3vect[k]!=(0,0,0):
                            vect4_list.append(convert3to4(list_3vect[k]))
                        else:
                            vect4_list.append((0,0,0,0))
                    
                    #make sure no erros while converting
                    if (writetofile1):
                        #define angle of rotation theta and rotate in 4 vector space eq 3.3 3.4 .3.5 ref https://arxiv.org/pdf/1704.02124.pdf
                        theta=math.atan2(vect4_list[0][2],vect4_list[0][3])+math.pi/2
                        for k in range(0,len(vect4_list)):
                            if vect4_list[k]!=(0,0,0,0):
                                py=vect4_list[k][2]*math.cos(theta)-vect4_list[k][3]*math.sin(theta)
                                pz=vect4_list[k][2]*math.sin(theta)+vect4_list[k][3]*math.cos(theta)
                                vector=(vect4_list[k][0],vect4_list[k][1],py,pz)
                                vect4_list[k]=vector            
                        #convert back to 3 vectors
                        for k in range(0,len(vect4_list)):
                            if vect4_list[k]!=(0,0,0,0):
                                list_3vect[k+1]=convert4to3(vect4_list[k])
                #makes ure no erros
                        if (writetofile):
                            finallist=[row[0]]
                            #only use 30 constituents
                            for k in range(0,len(list_3vect)):
                                finallist.append(list_3vect[k][0]/800)
                                finallist.append(list_3vect[k][1])
                                finallist.append(list_3vect[k][2])
                            wr.writerow(finallist)

            
            
            print("errors 4to3: ", e4to3)
            print("errors 3to4: ", e3to4)  
            

    
