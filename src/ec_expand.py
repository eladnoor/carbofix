#!/usr/bin/python2.5
################################################################################
# A script for getting the inchi file from KEGG and translating it to a giant
# SDF (a multi-molfile) so that it can be used with WebGCM to get the 
# free energy of formation for each compound in KEGG
################################################################################

import kegg
import pybel
import openbabel as ob
import oasa
import sys
import html_writer
import math

################################################################################
def GetNeighborIndices(obmol, atom_index):
    neighbors = []    
    for i in range(1, obmol.NumAtoms()+1):
        if (i == atom_index):
            continue
        if (obmol.GetBond(i, atom_index) != None):
            neighbors.append(i)
    return neighbors

def AtomicNum2Formula(atomic_num):
    if (atomic_num == 1):
        return 'H'
    if (atomic_num == 6):
        return 'C'
    if (atomic_num == 7):
        return 'N'
    if (atomic_num == 8):
        return 'O'
    if (atomic_num == 15):
        return 'P'
    if (atomic_num == 16):
        return 'S'
    return 'X'

def BondOrder2Formula(order):
    if (order == 1):
        return '-'
    if (order == 2):
        return '='
    if (order == 3):
        return '#'
    else:
        return ''

def Index2Formula(obmol, atom_index):
    atomic_num = obmol.GetAtom(atom_index).GetAtomicNum()
    return AtomicNum2Formula(atomic_num)

def GetNeighbors(obmol, atom_index):
    neighbors = []    
    for i in range(1, obmol.NumAtoms()+1):
        if (i == atom_index):
            continue
        bond = obmol.GetBond(i, atom_index)
        if (bond != None):
            neighbors.append(BondOrder2Formula(bond.GetBondOrder()) + Index2Formula(obmol, i))
    return neighbors

def IsHydrocarbon(obmol, atom_index):
    """
        Checks whether an atom is a carbon surrounded by hydrogen and carbon atoms (no oxygen, nitrogen, etc.)
        Also if there is a double bond to another carbon atom, this is not good enough.
        If it returns True, it means this atom is a candidate for adding an hydroxyl
    """
    a = obmol.GetAtom(atom_index)
    if (not a.IsCarbon()):
        return False
    carbon_neighbors_counter = 0
    for neighbor in GetNeighbors(obmol, atom_index):
        if (neighbor == '-C'):
            carbon_neighbors_counter += 1
        elif (neighbor == '-H'):
            pass
        else:
            return False # neighbors which are not carbon or hydrogen are not allowed
    return (carbon_neighbors_counter < 4) # if there are 4 carbon neighbors, there is no place for an hydroxyl group

def AnalyseInchi(inchi):
    molecule = pybel.readstring('inchi', inchi)
    obmol = molecule.OBMol
    obmol.DeleteHydrogens()
    
    #print obmol.GetFormula() + ":"
    #molecule.draw()
    
    new_mols = []
    for i in range(1, obmol.NumAtoms()+1):
        #print i, Index2Formula(obmol, i), str(GetNeighbors(obmol, i)),
        carbon_atom = obmol.GetAtom(i)
        if (IsHydrocarbon(obmol, i)):
            new_obmol = ob.OBMol(obmol)
            new_oxygen_atom = new_obmol.NewAtom()
            if (new_oxygen_atom == None):
                raise Exception("Cannot add an hydroxyl group to carbon on index %d" % i)

            new_index = new_oxygen_atom.GetIdx()
            # For some unknown reason, OBMol connects new atoms to some random existing ones.
            # Here we take care of that by removing all the existing bonds
            for j in GetNeighborIndices(new_obmol, new_index):
                new_obmol.DeleteBond(new_obmol.GetBond(j, new_index))
            new_oxygen_atom.SetAtomicNum(8) # oxygen
            if (not new_obmol.AddBond(new_index, i, 1)):
                raise Exception("Cannot add a bond to new hydroxyl group")
            
            new_mol = pybel.Molecule(new_obmol)
            new_mols.append(new_mol)
    
    return (molecule, new_mols)

def Mol2Svg(mol, filename):
    """Create a 2D depiction of the molecule.

    Optional parameters:
      filename -- write to file (default is None)

    OASA is used for 2D coordinate generation and depiction. Tkinter and
    Python Imaging Library are required for image display.
    """
    etab = ob.OBElementTable()
    oasa_mol = oasa.molecule()
    for atom in mol.atoms:
        v = oasa_mol.create_vertex()
        v.symbol = etab.GetSymbol(atom.atomicnum)
        v.charge = atom.formalcharge
        oasa_mol.add_vertex(v)

    for bond in ob.OBMolBondIter(mol.OBMol):
        e = oasa_mol.create_edge()
        e.order = bond.GetBO()
        if bond.IsHash():
            e.type = "h"
        elif bond.IsWedge():
            e.type = "w"
        oasa_mol.add_edge(bond.GetBeginAtomIdx() - 1, bond.GetEndAtomIdx() - 1, e)
    # I'm sure there's a more elegant way to do the following, but here goes...
    # let's set the stereochemistry around double bonds
    mol.write("smi") # Perceive UP/DOWNness
    for bond in ob.OBMolBondIter(mol.OBMol):
        ends = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        if bond.GetBO() == 2:
            stereobonds = [[b for b in ob.OBAtomBondIter(mol.OBMol.GetAtom(x)) if b.GetIdx() != bond.GetIdx() and (b.IsUp() or b.IsDown())]
                           for x in ends]
            if stereobonds[0] and stereobonds[1]: # Needs to be defined at either end
                if stereobonds[0][0].IsUp() == stereobonds[1][0].IsUp():
                    # Either both up or both down
                    stereo = oasa.stereochemistry.cis_trans_stereochemistry.SAME_SIDE
                else:
                    stereo = oasa.stereochemistry.cis_trans_stereochemistry.OPPOSITE_SIDE
                atomids = [(b[0].GetBeginAtomIdx(), b[0].GetEndAtomIdx()) for b in stereobonds]
                extremes = []
                for id, end in zip(ends, atomids):
                    if end[0] == id:
                        extremes.append(end[1])
                    else:
                        extremes.append(end[0])
                center = oasa_mol.get_edge_between(oasa_mol.atoms[ends[0] - 1], oasa_mol.atoms[ends[1] - 1])
                st = oasa.stereochemistry.cis_trans_stereochemistry(
                          center = center, value = stereo,
                          references = (oasa_mol.atoms[extremes[0] - 1], oasa_mol.atoms[ends[0] - 1],
                                        oasa_mol.atoms[ends[1] - 1], oasa_mol.atoms[extremes[1] - 1]))
                oasa_mol.add_stereochemistry(st)
    
    oasa_mol.remove_unimportant_hydrogens()
    oasa.coords_generator.calculate_coords(oasa_mol, bond_length=30)
    
    maxx = max([v.x for v in oasa_mol.vertices])
    minx = min([v.x for v in oasa_mol.vertices])
    maxy = max([v.y for v in oasa_mol.vertices])
    miny = min([v.y for v in oasa_mol.vertices])
    maxcoord = max(maxx - minx, maxy - miny)
    fontsize = 16
    bondwidth = 6
    linewidth = 2
    if maxcoord > 270: # 300  - margin * 2
        for v in oasa_mol.vertices:
            v.x *= 270. / maxcoord
            v.y *= 270. / maxcoord
        fontsize *= math.sqrt(270. / maxcoord)
        bondwidth *= math.sqrt(270. / maxcoord)
        linewidth *= math.sqrt(270. / maxcoord)

    canvas = oasa.cairo_out.cairo_out()
    canvas.show_hydrogens_on_hetero = True
    canvas.font_size = fontsize
    canvas.bond_width = bondwidth
    canvas.line_width = linewidth
    canvas.mol_to_cairo(oasa_mol, filename, format="svg")
   
################################################################################

################################################################################
cid2inchi = {}

#cid2inchi['C00091'] = 'InChI=1/C25H40N7O19P3S/c1-25(2,20(38)23(39)28-6-5-14(33)27-7-8-55-16(36)4-3-15(34)35)10-48-54(45,46)51-53(43,44)47-9-13-19(50-52(40,41)42)18(37)24(49-13)32-12-31-17-21(26)29-11-30-22(17)32/h11-13,18-20,24,37-38H,3-10H2,1-2H3,(H,27,33)(H,28,39)(H,34,35)(H,43,44)(H,45,46)(H2,26,29,30)(H2,40,41,42)/t13-,18-,19-,20?,24-/m1/s1'
#cid2inchi['C00100'] = 'InChI=1/C24H40N7O17P3S/c1-4-15(33)52-8-7-26-14(32)5-6-27-22(36)19(35)24(2,3)10-45-51(42,43)48-50(40,41)44-9-13-18(47-49(37,38)39)17(34)23(46-13)31-12-30-16-20(25)28-11-29-21(16)31/h11-13,17-19,23,34-35H,4-10H2,1-3H3,(H,26,32)(H,27,36)(H,40,41)(H,42,43)(H2,25,28,29)(H2,37,38,39)/t13-,17-,18-,19?,23-/m1/s1'
#cid2inchi['C00827'] = 'InChI=1/C24H40N7O18P3S/c1-12(32)23(37)53-7-6-26-14(33)4-5-27-21(36)18(35)24(2,3)9-46-52(43,44)49-51(41,42)45-8-13-17(48-50(38,39)40)16(34)22(47-13)31-11-30-15-19(25)28-10-29-20(15)31/h10-13,16-18,22,32,34-35H,4-9H2,1-3H3,(H,26,33)(H,27,36)(H,41,42)(H,43,44)(H2,25,28,29)(H2,38,39,40)/t12?,13-,16-,17-,18?,22-/m1/s1'
#cid2inchi['C00022'] = 'InChI=1/C3H4O3/c1-2(4)3(5)6/h1H3,(H,5,6)'

#molecule = pybel.readstring('inchi', cid2inchi['C00100'])
#obmol = molecule.OBMol
#pybel.Molecule(obmol).draw()
#a = obmol.NewAtom()
#a.SetAtomicNum(8)
#obmol.AddBond(2, a.GetIdx(), 1)
#pybel.Molecule(obmol).draw()
#sys.exit(0)

kegg = kegg.Kegg()
kegg.prepare_database()

if (len(cid2inchi) == 0):
    inchi_file = open(kegg.INCHI_FILE, 'r')
    sd_file = open('../kegg/compounds.sdf', 'w')
    for line in inchi_file.readlines():
        (cid, inchi) = line.strip().split()
        cid2inchi[cid] = inchi
    sd_file.close()

inchi2cid = {}
inchi_pairs = set()
for (cid, inchi) in cid2inchi.iteritems():
    (subs, prod_list) = AnalyseInchi(inchi)
    for prod in prod_list:
        inchi_pairs.add((subs.write("INCHI"), prod.write("INCHI")))

    if (subs.write("INCHI") in inchi2cid):
        continue
    else:
        inchi2cid[subs.write("INCHI")] = cid

res_file = open('../kegg/reaction_extra.txt', 'w')
html_file = html_writer.HtmlWriter('../res/ec_expand.html')
rcounter = 0
for (s, p) in inchi_pairs:
    print p
    if (p in inchi2cid):
        cid_s = inchi2cid[s]
        cid_p = inchi2cid[p]
        res_file.write("NERW M%05d: " % rcounter + cid_s + " => " + cid_p + "\n")
        html_file.write('<p><img name="Substrate" src="http://www.kegg.com/Fig/compound/%s.gif" border=0>' % cid_s)
        html_file.write(' ... ')
        html_file.write('<img name="Product" src="http://www.kegg.com/Fig/compound/%s.gif" border=0></p>' % cid_p)
        rcounter += 1
        
res_file.close()
