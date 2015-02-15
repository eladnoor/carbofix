#!/usr/bin/python2.5
################################################################################
# A script for getting the inchi file from KEGG and translating it to a giant
# SDF (a multi-molfile) so that it can be used with WebGCM to get the 
# free energy of formation for each compound in KEGG
################################################################################

import kegg
import pybel

kegg = kegg.Kegg()
kegg.prepare_database()

inchi_file = open(kegg.INCHI_FILE, 'r')
sd_file = open('../kegg/compounds.sdf', 'w')
for line in inchi_file.readlines():
    (cid, inchi) = line.strip().split()
    molecule = pybel.readstring('inchi', inchi)
    mol = molecule.write('sdf')
    molfile_lines = mol.split('\n')
    if (len(molfile_lines) < 4):
        print "ERROR: " + cid
        continue
    molfile_lines[1] = " " + cid
    #print mol.calcdesc()
    sd_file.write('\n'.join(molfile_lines))

sd_file.close()

