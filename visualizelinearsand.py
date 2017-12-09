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
import sys


n=1
nunstable=1
delta=1
curve=[]
unstable=[]
radiusx=2
radiusy=2

def draw():
    canvas.delete(ALL)
    for [i,j] in curve:
        canvas.create_rectangle(delta*i, delta*j, delta*(i+1), delta*(j+1), fill="green")
    for [i,j] in unstable:
        canvas.create_oval(delta*(i+0.5)+radiusx,delta*(j+0.5)+radiusx,delta*(i+0.5)-radiusy,delta*(j+0.5)-radiusy,fill="blue")
    canvas.update_idletasks()       




root = Tk()
frame = Frame(root, bd=1, relief=SUNKEN)
frame.grid_rowconfigure(0, weight=1)
frame.grid_columnconfigure(0, weight=1)
xscroll = Scrollbar(frame, orient=HORIZONTAL)
xscroll.grid(row=1, column=0, sticky=E+W)
yscroll = Scrollbar(frame)
yscroll.grid(row=0, column=1, sticky=N+S)
canvas = Canvas(frame,bd=0, xscrollcommand=xscroll.set, yscrollcommand=yscroll.set)
canvas.focus_set()


#this assumes that the integer size is 4 bytes
if __name__ == "__main__":
    with  open("./tsandpile/grid.dat","rb") as input:
        n = unpack("i",input.read(4))[0]
        nunstable = unpack("i",input.read(4))[0]
        curvesize = unpack("i",input.read(4))[0]    
        curve = [ [unpack("i",input.read(4))[0],unpack("i",input.read(4))[0]] for i in xrange(curvesize)]
        unstable = [ [unpack("i",input.read(4))[0],unpack("i",input.read(4))[0]] for i in xrange(nunstable)]

    delta=float(700)/n
    canvas.config(width=delta*n, height=delta*n)
    canvas.grid(row=0, column=0, sticky=N+S+E+W)
    xscroll.config(command=canvas.xview)
    yscroll.config(command=canvas.yview)
    frame.pack(fill=BOTH,expand=0)
    draw()
    root.mainloop()
    

    
