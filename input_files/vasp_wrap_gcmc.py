#!/usr/bin/env python

# ----------------------------------------------------------------------
# LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
# http://lammps.sandia.gov, Sandia National Laboratories
# Steve Plimpton, sjplimp@sandia.gov
# ----------------------------------------------------------------------

# Syntax: vasp_wrap_gcmc.py file/zmq POSCARfile

# wrapper on VASP to act as server program using CSlib
#   receives message with list of coords from client
#   creates VASP inputs
#   invokes VASP to calculate self-consistent energy of that config
#   reads VASP outputs
#   sends message with energy and coords (if gcmc move is accepted) to client

# NOTES:
# this code archives VASP input/output in separate directories
# ideally, poscar_write() is the only function user would need to modify as per their simulation cell 
# possible attributes to include in future:
# check to ensure basic VASP input files are in place?
# could make syntax for launching VASP more flexible
#   e.g. command-line arg for # of procs
# detect if VASP had an error and return ERROR field, e.g. non-convergence ??

from __future__ import print_function
import sys

version = sys.version_info[0]
if version == 3:
  sys.exit("The CSlib python wrapper does not yet support python 3")
  
import subprocess
import xml.etree.ElementTree as ET
from cslib import CSlib
import logging
logging.basicConfig(filename='log.vasp',level=logging.DEBUG)

# command to call vasp binary

vaspcmd = "ibrun -np 220 /home1/apps/intel18/impi18_0/vasp/5.4.4/bin/vasp_std > vasp.out"

# enums matching FixClientGCMC class in LAMMPS

SETUP,STEP,UPDATE = range(1,3+1)
DIM,PERIODICITY,ORIGIN,BOX,NATOMS,NTYPES,TYPES,COORDS,UNITS,CHARGE = range(1,10+1)
ENERGY,ERROR = range(1,2+1)  

# -------------------------------------
# functions

# error message and exit

def error(txt):
  logging.debug("ERROR:%s",txt)
  sys.exit(1)

# -------------------------------------
# read initial VASP POSCAR file to setup problem
# return natoms,ntypes,box

def vasp_setup(poscar):

  ps = open(poscar,'r').readlines()

  # box size

  words = ps[2].split()
  xbox = float(words[0])
  words = ps[3].split()
  ybox = float(words[1])
  words = ps[4].split()
  zbox = float(words[2])
  box = [xbox,ybox,zbox]

  ntypes = 0
  natoms = 0
  words = ps[6].split()
  for word in words:
    if word == '#': break
    ntypes += 1
    natoms += int(word)
  
  return natoms,ntypes,box
  
# -------------------------------------
# write a new POSCAR file for VASP

def poscar_write(poscar,natoms,ntypes,types,coords,box,vasprun,natoms_types):

  psold = open(poscar,'r').readlines()
  psnew = open("POSCAR",'w')

  # header, including box size
  
  psnew.write(psold[0])
  psnew.write(psold[1])
  psnew.write("%g %g %g\n" % (box[0],box[1],box[2]))
  psnew.write("%g %g %g\n" % (box[3],box[4],box[5]))
  psnew.write("%g %g %g\n" % (box[6],box[7],box[8]))
  psnew.write(psold[5]) # note that this works only if same types of atoms always present in cell, modify as necessary
  psnew.write("%g %g\n" % (natoms_types[0],natoms_types[1])) # modify as necessary

  # per-atom coords
  # grouped by types

  psnew.write("Cartesian\n")

  for itype in range(1,ntypes+1): 
    for i in range(natoms):
      if types[i] != itype: continue
      x = coords[3*i+0]
      y = coords[3*i+1]
      z = coords[3*i+2]
      aline = "  %g %g %g\n" % (x,y,z)
      psnew.write(aline)

  psnew.close()

  # keeping track of POSCAR at every vasprun

  psnew_ts = open("POSCAR_%s" % str(vasprun),'w')
  psnew_ts.write(psold[0])
  psnew_ts.write(psold[1])
  psnew_ts.write("%g %g %g\n" % (box[0],box[1],box[2]))
  psnew_ts.write("%g %g %g\n" % (box[3],box[4],box[5]))
  psnew_ts.write("%g %g %g\n" % (box[6],box[7],box[8]))
  psnew_ts.write(psold[5]) # note that this works only if same types of atoms always present in cell, modify as necessary
  psnew_ts.write("%g %g\n" % (natoms_types[0],natoms_types[1])) # modify as necessary
 
  # per-atom coords
  # grouped by types

  psnew_ts.write("Cartesian\n")

  for itype in range(1,ntypes+1): 
    for i in range(natoms):
      if types[i] != itype: continue
      x = coords[3*i+0]
      y = coords[3*i+1]
      z = coords[3*i+2]
      aline = "  %g %g %g\n" % (x,y,z)
      psnew_ts.write(aline)

  psnew_ts.close()

# -------------------------------------
# read a VASP output vasprun.xml file
# uses ElementTree module
# see https://docs.python.org/2/library/xml.etree.elementtree.html

def vasprun_read():
  tree = ET.parse('vasprun.xml')
  root = tree.getroot()
  
  #fp = open("vasprun.xml","r")
  #root = ET.parse(fp)
  
  scsteps = root.findall('calculation/scstep')
  energy = scsteps[-1].find('energy')
  for child in energy:
    if child.attrib["name"] == "e_0_energy":
      eout = float(child.text)
  return eout

# -------------------------------------
# read CONTCAR to send to LAMMPS
# remap and convert to Cartesian coordinates; https://www.ruppweb.org/Xray/tutorial/Coordinate%20system%20transformation.htm

def read_contcar(natoms, box, origin):  
  ps = open("CONTCAR",'r').readlines()
  h = [box[0], box[4], box[8], box[7], box[6], box[3]] # h = [xprd, yprd, zprd, yz, xz, xy]
  boxlo = [origin[0], origin[1], origin[2]]
  relaxed_coords = []
  for i in range(8,natoms+8): 
    coords_atom = ps[i].split()
    lamda_x = float(coords_atom[0])
    lamda_y = float(coords_atom[1])
    lamda_z = float(coords_atom[2])
    x = lamda_x; y = lamda_y; z = lamda_z;
    if lamda_x < 0:
      x = lamda_x + 1
    if lamda_x >= 1:
      x = lamda_x - 1
    if lamda_y < 0:
      y = lamda_y + 1
    if lamda_y >= 1:
      y = lamda_y - 1
    if lamda_z < 0:
      z = lamda_z + 1
    if lamda_z >= 1:
      z = lamda_z - 1
    x_cart = h[0]*x + h[5]*y + h[4]*z + boxlo[0];
    y_cart = h[1]*y + h[3]*z + boxlo[1];
    z_cart = h[2]*z + boxlo[2];
    relaxed_coords.append(x_cart)
    relaxed_coords.append(y_cart)
    relaxed_coords.append(z_cart)
  return relaxed_coords

# -------------------------------------
# track VASP run number

def tracking_vasp_runs(vasprun):
  vasprun = vasprun + 1 
  logging.debug("\n")
  logging.debug("vasp run number=%s",vasprun) 
  logging.debug("ntypes =%s",ntypes)
  return vasprun

# -------------------------------------
# reorder coordinates from LAMMPS for VASP POTCAR

def reordering_coords(types,coords):
  thisdict={}
  for i in range (len(types)):
    thisdict.update({types[i]:[]})
  for i in range (len(types)):
    thisdict[types[i]].extend(coords[3*i:3*i+3])
  coords=[]
  for key in thisdict:
      coords.extend(thisdict[key])
  return coords

# -------------------------------------
# find number of each type of atom for POSCAR

def finding_natoms_types(types):
  from itertools import groupby
  types.sort()
  natoms_types = [len(list(group)) for key, group in groupby(types)]
  logging.debug("natoms_types = %s",natoms_types)

  natoms = sum(natoms_types)  
  logging.debug("natoms = %s",natoms)
  return types, natoms_types, natoms

# -------------------------------------
# archive VASP run in directories  

def archiving_vasprun(vasprun):
  import os
  # define the name of the directory to be created
  name = str(vasprun)
  dirpath = os.path.join('./', name)
  try:
    os.mkdir(dirpath)
  except OSError:
    logging.debug ("Creation of the directory %s failed" % dirpath)
  else:
    logging.debug ("Successfully created the directory %s " % dirpath)

  import shutil
  src = os.path.join('./', 'POSCAR'); shutil.copy(src, dirpath)
  src = os.path.join('./', 'CHG'); shutil.copy(src, dirpath)
  src = os.path.join('./', 'CHGCAR'); shutil.copy(src, dirpath)
  src = os.path.join('./', 'CONTCAR'); shutil.copy(src, dirpath)
  src = os.path.join('./', 'DOSCAR'); shutil.copy(src, dirpath)
  src = os.path.join('./', 'EIGENVAL'); shutil.copy(src, dirpath)
  src = os.path.join('./', 'IBZKPT'); shutil.copy(src, dirpath)
  src = os.path.join('./', 'INCAR'); shutil.copy(src, dirpath)
  src = os.path.join('./', 'KPOINTS'); shutil.copy(src, dirpath)
  src = os.path.join('./', 'OSZICAR'); shutil.copy(src, dirpath)
  src = os.path.join('./', 'OUTCAR'); shutil.copy(src, dirpath)
  src = os.path.join('./', 'PCDAT'); shutil.copy(src, dirpath)
 # src = os.path.join('./', 'POTCAR'); shutil.copy(src, dirpath)
  src = os.path.join('./', 'REPORT'); shutil.copy(src, dirpath)
  src = os.path.join('./', 'vasp.out'); shutil.copy(src, dirpath)
  src = os.path.join('./', 'vasprun.xml'); shutil.copy(src, dirpath)
  src = os.path.join('./', 'WAVECAR'); shutil.copy(src, dirpath)
  src = os.path.join('./', 'XDATCAR'); shutil.copy(src, dirpath)

# -------------------------------------
# main program

# command-line args

if len(sys.argv) != 3:
  logging.debug("Syntax: python vasp_wrap_gcmc.py file/zmq POSCARfile")
  sys.exit(1)

mode = sys.argv[1]
poscar_template = sys.argv[2]

# instantiating CSLib in server mode in serial usage

if mode == "file": cs = CSlib(1,mode,"tmp.couple",None) 
elif mode == "zmq": cs = CSlib(1,mode,"*:5555",None)
else:
  logging.debug("Syntax: python vasp_wrap_gcmc.py file/zmq POSCARfile")
  sys.exit(1)
'''
from mpi4py import MPI # Python wrapper on MPI
world=MPI.COMM_WORLD
if mode == "file": cs = CSlib(1,mode,"tmp.couple",world)
elif mode == "zmq": cs = CSlib(1,mode,"*:5555",world)
else:
  print("Syntax: python vasp_wrap_gcmc.py file/zmq POSCARfile")
  sys.exit(1)
'''
natoms,ntypes,box = vasp_setup(poscar_template)

# initial message for GCMC protocol
# receives message with msgID = 0 and protocol from client 
# sends response message with msgID = 0 to complete initial handshake

msgID,nfield,fieldID,fieldtype,fieldlen = cs.recv()
if msgID != 0: error("Bad initial client/server handshake")
protocol = cs.unpack_string(1)
if protocol != "gcmc": error("Mismatch in client/server protocol")
cs.send(0,0)

#timestep = -1 
# initializing before server loop begins, to keep track of VASP runs in separate directories

import glob
dir_list = glob.glob('**/')
if len(dir_list) == 0:
  vasprun = -1
else:
  dir_list_int = []
  for i in range(len(dir_list)):
    dir_list_int.append(int(dir_list[i].strip('/')))
    vasprun = max(dir_list_int) - 1 # -1, as during restart, previous step is run again during setup

# endless server loop

while 1:

  # recv message from client
  # msgID = 0 = all-done message

  msgID,nfield,fieldID,fieldtype,fieldlen = cs.recv()
  if msgID < 0: break

  # SETUP receive at beginning of each run
  # required fields: DIM, PERIODICTY, ORIGIN, BOX, 
  #                  NATOMS, NTYPES, TYPES, COORDS,
  # optional fields: others in enum above, but VASP ignores them

  if msgID == SETUP:
    
    origin = []
    box = []
    natoms_recv = ntypes_recv = 0
    types = []
    coords = []
    
    for field in fieldID:
      if field == DIM:
        dim = cs.unpack_int(DIM)
        if dim != 3: error("VASP only performs 3d simulations")
      elif field == PERIODICITY:
        periodicity = cs.unpack(PERIODICITY,1)
        if not periodicity[0] or not periodicity[1] or not periodicity[2]:
          error("VASP wrapper only currently supports fully periodic systems")
      elif field == ORIGIN:
        origin = cs.unpack(ORIGIN,1)
      elif field == BOX:
        box = cs.unpack(BOX,1) #box = [xprd, 0, 0, xy, yprd, 0, xz, yz, zprd]
      elif field == NATOMS:
        natoms_recv = cs.unpack_int(NATOMS)
        natoms = natoms_recv
        #if natoms != natoms_recv:
          #error("VASP wrapper mis-match in number of atom") 
      elif field == NTYPES:
        ntypes_recv = cs.unpack_int(NTYPES)
        ntypes = ntypes_recv
        #if ntypes != ntypes_recv: 
          #error("VASP wrapper mis-match in number of atom types") 
      elif field == TYPES:
        types = cs.unpack(TYPES,1)
      elif field == COORDS:
        coords = cs.unpack(COORDS,1)

    if not origin or not box or not natoms or not ntypes or \
       not types or not coords:
      error("Required VASP wrapper setup field not received");

    # keeping track in log.vasp
    vasprun = tracking_vasp_runs(vasprun)

    # re-ordering coords for POTCAR
    # this ensures ordering of individual atoms remains same
    coords = reordering_coords(types,coords)

    # finding number of atoms per type
    types, natoms_types, natoms = finding_natoms_types(types)   

    # create POSCAR file 
    # with sorted types, ordered coords, updated # atoms per type
    poscar_write(poscar_template,natoms,ntypes,types,coords,box,vasprun,natoms_types)

    # invoke VASP
    
    logging.debug("Launching VASP ...")
    logging.debug("%s",vaspcmd)
    subprocess.check_output(vaspcmd,stderr=subprocess.STDOUT,shell=True)
    
    # process VASP output; throw exception if VASP encounters error to ensure job continuation
    # in the future, could possibly link this to Custodian to navigate errors better

    try:
      energy = vasprun_read()
      relaxed_coords = read_contcar(natoms, box, origin) # read CONTCAR here from relaxed_coords
    except:
      logging.debug("VASP encountered error during minimization, sending arbitrary large energy to LAMMPS")
      energy = 100000000000000000000000000000000
    logging.debug("VASP energy=%s",energy)
    
    cs.send(msgID,1);
    cs.pack_double(ENERGY,energy)

    # archiving VASP input and outputs at each vasprun to directories to keep track
    archiving_vasprun(vasprun)

  # STEP receive at each timestep of run or minimization
  # required fields: COORDS, NATOMS, NTYPES, TYPES 
  # optional fields: ORIGIN, BOX (LAMMPS will send these if simulation box changes)

  elif msgID == STEP:

    coords = []
    natoms = ntypes = 0 
    types = [] 
    
    for field in fieldID:
      if field == COORDS:
        coords = cs.unpack(COORDS,1)
      elif field == ORIGIN:
        origin = cs.unpack(ORIGIN,1)
      elif field == BOX:
        box = cs.unpack(BOX,1) #box = [xprd, 0, 0, xy, yprd, 0, xz, yz, zprd]
      elif field == NATOMS:
        natoms = cs.unpack_int(NATOMS)
      elif field == NTYPES:
        ntypes = cs.unpack_int(NTYPES)
      elif field == TYPES:
        types = cs.unpack(TYPES,1)
    
    if not coords or not natoms or not ntypes or not types: error("Required VASP wrapper step field not received"); 

    # keeping track in log.vasp
    vasprun = tracking_vasp_runs(vasprun)

    # re-ordering coords for POTCAR
    # this ensures ordering of individual atoms remains same
    coords = reordering_coords(types,coords)

    # finding number of atoms per type
    types, natoms_types, natoms = finding_natoms_types(types)

    # create POSCAR file 
    # with sorted types, ordered coords, updated # atoms per type
    
    poscar_write(poscar_template,natoms,ntypes,types,coords,box,vasprun,natoms_types) 

    # invoke VASP
    
    logging.debug("Launching VASP ...")
    logging.debug("%s",vaspcmd)
    subprocess.check_output(vaspcmd,stderr=subprocess.STDOUT,shell=True)
    
    # process VASP output; throw exception if VASP encounters error to ensure job continuation
    # in the future, could possibly link this to Custodian to navigate errors better

    try:
      energy = vasprun_read()
      relaxed_coords = read_contcar(natoms, box, origin) # read CONTCAR here from relaxed_coords
    except:
      logging.debug("VASP encountered error during minimization, sending arbitrary large energy to LAMMPS")
      energy = 100000000000000000000000000000000
    logging.debug("VASP energy=%s",energy)
    
    cs.send(msgID,1);
    cs.pack_double(ENERGY,energy)

    # archiving VASP input and outputs at each vasprun to directories to keep track
    archiving_vasprun(vasprun)

  elif msgID == UPDATE:
    cs.send(UPDATE,2)
    cs.pack(COORDS,4,3*natoms,relaxed_coords)  #all of VASP's processors send a copy of the entire field anyway
    cs.pack(TYPES,1,natoms,types)    

  else: error("VASP wrapper received unrecognized message")

# final reply to client
  
cs.send(0,0)

# clean-up

del cs

