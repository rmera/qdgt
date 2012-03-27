#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#       untitled.py
#       
#       Copyright 2012 Raul Mera <rmera{at}chemDOThelsinkiDOTfi>
#       
#       This program is free software; you can redistribute it and/or modify
#       it under the terms of the GNU General Public License as published by
#       the Free Software Foundation; either version 2 of the License, or
#       (at your option) any later version.
#       
#       This program is distributed in the hope that it will be useful,
#       but WITHOUT ANY WARRANTY; without even the implied warranty of
#       MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#       GNU General Public License for more details.
#       
#       You should have received a copy of the GNU General Public License
#       along with this program; if not, write to the Free Software
#       Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#       MA 02110-1301, USA.
#       
#       

#To the long life of the Ven. Khenpo Phuntsok Tenzin Rinpoche.

#This is a simple ORCA to Turbomole, Turbomole to ORCA translator.
#instead of messing with C, This will fool gromacs into running TM 
#calculations. Anyway, the idea is that you use an only
#slightly modified version of qm_orca.c to run QMMM with turbomole.
#It is of course less efficient than making a new qm_turbomole.c
#for gromacs, including the need for more disk access.
#It makes things much easier to code, which can be good to build interfaces
#for other QM software. 
#
# The extra cost should be negligible compared with the QM energy and 
#gradient calculation. Maybe if the QM system is very small and the
#MM system very big, one could see some performance hazard compared
#to Gromacs/ORCA, for example.



import sys


#takes ORCA coordinates for atoms and point charges,
#transforms them into TM format.
#The user is supposed to prepare everything else
#in the control file, including the "point charges" option 
#in $drvopt. Also,  after preparing the control file,
#the user must move o copy it to a file named control-template
#in the control file.  for regular QMMM reading TM coordinates
#is not necesary so is not implemented now.


atomn2symbol={"1":"h","6":"c","7":"n","8":"o","29":"cu","30":"zn","16":"s","17":"cl","11":"na"} #incomplete

#Angstrom 2 Bohr
def a2b(angstrom):
	return angstrom/0.529177249
#Bohr 2 Angstrom
def b2a(bohr):
	return bohr*0.529177249

class coordorcatranslator():
	def __init__(self):
		self.charges=[]
		self.coords=[]
		pass
	def readorcacoords(self,orcainp):
		orca=open(orcainp,"r")
		reading=False #True when the coordinates begin.
		for i in orca:
			if "*" in i:
				#The coordinates in orca begin and end with "*".
				reading= not(reading)
			elif reading:
				fields=i.split()
				#Symbol x y z.
				self.coords.append([fields[0],float(fields[1]),float(fields[2]),float(fields[3])])
		orca.close()
	def readorcacharges(self,chargesname):
		charges=open(chargesname,"r")
		charges.readline() #We don't need this.
		for i in charges:
			if i == "\n":
				break
			fields=i.split()
			#Charge x y z.
			self.charges.append([float(fields[0]),float(fields[1]),float(fields[2]),float(fields[3])])
	def writetmcoords(self):
		tmcoords=open("coord","w")
		tmcoords.write("$coord\n")
		for i in self.coords:
			tmcoords.write("{0:20.14f}   {1:20.14f}   {2:20.14f}     {3:2s}\n".format(a2b(i[1]),a2b(i[2]),a2b(i[3]),atomn2symbol[i[0]]))
		tmcoords.write("$user-defined-bonds\n$end\n")
		tmcoords.close()
	def writetmcharges(self):
		template=open("control-template","r")
		target=open("control","w")
		for i in template:
			if "$grad" in i: #Copy the point charges after "$grad".
				target.write("$point_charges\n")
				for j in self.charges:
					target.write("{0:6.3f}   {1:6.3f}   {2:6.3f}     {3:6.3f}\n".format(a2b(j[1]),a2b(j[2]),a2b(j[3]),j[0]))
			target.write(i+"\n")

#Same as above but for gradients.
#Gradients should be both in the same units (Hatree/Bohr)
#This appears to also be the case with Turbomole.
class gradorcatranslator:
	def __init__(self):
		self.energy=0
		self.chargegradients=[]
		self.gradients=[] #Gradient is just a big list, close to orca format.
		pass
	def readtmgrads(self):
		tmgrad=open("gradient","r")
		tmgrad.readline()	# Don't need first line.	
		self.energy=float(tmgrad.readline().split()[6])
		for i in tmgrad:
			if "$end" in i:
				break
			fields=i.split()
			#In tubomole you have first the coordinates, which
			#has 4 fields so we easily skip them.
			if len(fields)!=3: 
				continue
			self.gradients.append(float(fields[0].replace("D","E")))
			self.gradients.append(float(fields[1].replace("D","E")))
			self.gradients.append(float(fields[2].replace("D","E")))
		tmgrad.close()
	def readtmchargegrads(self):
		control=open("control","r")
		readint=False
		for i in control:
			if "$point_charge_gradients" in i:
				readint=True
			elif readint==True:
				fields=i.split()
				if len(fields)!=3:
					break #No more gradients.
				self.chargegradients.append([float(fields[0].replace("D","E")),float(fields[1].replace("D","E")),float(fields[2].replace("D","E"))])
		control.close()
	def writeorcagradients(self,filename):
		orcagrad=open(filename,"w")
		for i in range(7):
			orcagrad.write("#Jackie, Sammo, Bolo\n") #we need not-to-be-read lines.
		orcagrad.write(str(self.energy)+"\n")
		for i in range(3):
			orcagrad.write("#Yaoh!\n") #3 more...
		for i in self.gradients:
			orcagrad.write("{0:22.12f}".format(i)+"\n")  #This should do the trick.
		orcagrad.close()
	def writeorcachargegrads(self,filename):
		outgrads=open(filename,"w")
		outgrads.write("#Bruce\n") #Bruce has his own line.
		for i in self.chargegradients:
			outgrads.write("{0:17.12f} {0:17.12f} {0:17.12f}\n".format(i[0],i[1],i[2]))
		outgrads.close()
		
	
if "-O2T" in sys.argv:
	translator=coordorcatranslator()
	translator.readorcacoords(sys.argv[1])
	translator.readorcacharges(sys.argv[2])
	translator.writetmcoords()
	translator.writetmcharges()
elif "-T2O" in sys.argv:
	translator=gradorcatranslator()
	translator.readtmgrads()
	translator.readtmchargegrads()
	translator.writeorcagradients(sys.argv[1])
	translator.writeorcachargegrads(sys.argv[2])



