# -*- coding: utf-8 -*-
"""
Created on Fri Oct  3 00:24:21 2014

@author: nikitakalinin
 """

from math import *
import random

#side of the square
n=100

#initial number of grains at every points
backphone=3

global field
field = []
global stack1
stack1 = []
global stack2
stack2 = []

#to measure avalanches
global touched
touched = []

for i in range(n):
    f = []
    for j in range(n):
        f.append(backphone)
    field.append(f)

def add(i,j):
    global field
    field[i][j]=field[i][j]+1
    global stack1
    stack1.append([i,j])

def relax():
    global touched
    touched = []
    for i in range(n):
        g=[]
        for j in range(n):
            g.append(0)
        touched.append(g)

    while len(stack1)>0 or len(stack2)>0:
        while len(stack1)>0:        
            [i,j] = stack1.pop()
            if field[i][j] > 3:
                field[i][j] = field[i][j]-4
                touched[i][j]=touched[i][j]+1
                for [k,l] in ([i-1,j],[i,j-1],[i+1,j],[i,j+1]):
                    if boundary(k,l)!=0:                    
                        field[k][l]=field[k][l]+1
                        if boundary(k,l)>0:
                            if field[k][l]>3:
                                stack2.append([k,l])
                                      
        while len(stack2)>0:        
            [i,j] = stack2.pop()
            if field[i][j] > 3:
                field[i][j] = field[i][j]-4
                touched[i][j]=touched[i][j]+1
                for [k,l] in ([i-1,j],[i,j-1],[i+1,j],[i,j+1]):
                    if boundary(k,l)!=0:                    
                        field[k][l]=field[k][l]+1
                        if boundary(k,l)>0:
                            if field[k][l]>3:
                                stack1.append([k,l])
                        

def boundary(i,j):
    return (i)*(j)*(n-i)*(n-j)


 
       

if __name__ == "__main__":
    #number of grains dropped
    k=1000000
    random.seed(1)
    f1 = open('classical_sandpile/'+str(backphone)+'avalanchesfinal'+str(n)+' '+str(k)+'.txt', 'w')
    f2 = open('classical_sandpile/'+str(backphone)+'volumeavalanchesfinal'+str(n)+' '+str(k)+'.txt', 'w')
    f1.write('side='+str(n)+' '+'number of points ='+ str(k)+'sizes of avalanches: ')
    f2.write('side='+str(n)+' '+'number of points ='+ str(k)+'volumes of avalanches: ')
    
    for l in range(k):
        add(random.randint(1,n-1),random.randint(1,n-1))
        relax()
        counts=0
        countv=0
        for i in range(n):
            for j in range(n):
                if touched[i][j]>0:
                    counts=counts+1
                    countv=countv+touched[i][j]
    
        f1.write(str(counts)+',')
        f2.write(str(countv)+',')
        print l
    f1.close()
    f2.close()


    
