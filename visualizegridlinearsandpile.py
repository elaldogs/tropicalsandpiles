# -*- coding: utf-8 -*-
#============================================================================
# Name        : visualizegrid.py
# Author      : Aldo Guzmán-Sáenz
# Version     :
# Copyright   : 
# Description : This program displays a graphical representation of the output
#               of linearsandpile.cpp
#============================================================================
from Tkinter import *
from struct import unpack

n=1
nunstable=1
delta=1
monomials=[]
unstable=[]

def draw():
    canvas.delete(ALL)
    for [i,j] in monomials:
        canvas.create_rectangle(delta*i, delta*j, delta*(i+1), delta*(j+1), fill="green")
    for [i,j] in unstable:
        canvas.create_oval(delta*i+10,delta*j+10,delta*i-8,delta*j-8,fill="blue")
    canvas.update_idletasks()       
    
root = Tk()
frame = Frame(root, bd=1, relief=SUNKEN)
frame.grid_rowconfigure(0, weight=1)
frame.grid_columnconfigure(0, weight=1)
scrollbarx = Scrollbar(frame, orient=HORIZONTAL)
scrollbarx.grid(row=1, column=0, sticky=E+W)
scrollbary = Scrollbar(frame)
scrollbary.grid(row=0, column=1, sticky=N+S)
canvas = Canvas(frame,bd=0, xscrollcommand=scrollbarx.set, yscrollcommand=scrollbary.set)
canvas.grid(row=0, column=0, sticky=N+S+E+W)
scrollbarx.config(command=canvas.xview)
scrollbary.config(command=canvas.yview)
frame.pack(fill=BOTH,expand=0)

#this assumes that the integer size is 4 bytes
if __name__ == "__main__":
    with  open("./grid.dat") as input:
        n = unpack("i",input.read(4))[0]
        nunstable = unpack("i",input.read(4))[0]
        nmonomials = unpack("i",input.read(4))[0]
        monomials = [ [unpack("i",input.read(4))[0],unpack("i",input.read(4))[0]] for i in xrange(nmonomials)]
        unstable = [ [unpack("i",input.read(4))[0],unpack("i",input.read(4))[0]] for i in xrange(nunstable)]
    delta=max([1,800/n])
    canvas.config(width=delta*n, height=delta*n)
    draw()
    root.mainloop()
    

    
