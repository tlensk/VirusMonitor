# -*- coding: utf-8 -*-
"""
Created on Sat Sep  5 15:46:18 2020

@author: tilen
"""

sep = ","
sep2 = "\t"

fInName = "40_corona.csv"
fIn = open(fInName, "r")
lines = fIn.readlines()
fIn.close()

header = lines[0]
t_header = header.split(sep)
new_header = sep2.join(t_header)
lines = lines[1:]

rows = t_header[0]

body = ""

for line in lines:
    t_line = line.split(sep)
    #print(t_line)
    rows = rows+sep2+t_line[0]
    body = body+sep2.join(t_line[1:])

fOut = open("spltrows","w")
fOut.write(rows[:(-1)])
fOut.close()

fOut = open("spltout","w")
fOut.write(body)
fOut.close()

fOut = open("splthead","w")
fOut.write(new_header)
fOut.close()

fOut = open("splttitle","w")
fOut.write("1"+"\n")
fOut.write("1"+"\n")
fOut.write(fInName+"\n")
fOut.close()
