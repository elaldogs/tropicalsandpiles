# -*- coding: utf-8 -*-
#============================================================================
# Name        : visualizeactive.py
# Author      : Aldo Guzmán-Sáenz
# Version     :
# Copyright   : 
# Description : This program processes the file active.dat for human reading
#============================================================================
from struct import unpack


#this assumes that the integer size is 4 bytes
if __name__ == "__main__":
    with  open("./active.dat","rb") as input:
        pseqlen = unpack("i",input.read(4))[0]
        print pseqlen
        points = [(unpack("i",input.read(4))[0],unpack("i",input.read(4))[0],unpack("i",input.read(4))[0]) for i in xrange(pseqlen )]
    print points
    

    
