#!/usr/bin/python2.5
#
# Except for the default python packages for python 2.5, this program requires:
# pydot - a python interface for graphviz (which requires the libraries graphviz and networkx)
################################################################################

import os
import sys
import re
import util
import pydot
import csv
import pylab
import gzip

def parse_reaction_formula_side(s):
    """ parse the side formula, e.g. '2 C00001 + C00002 + 3 C00003'
        return the set of CIDs, ignore stoichiometry
    """
    compound_bag = {}
    for member in re.split('\s+\+\s+', s):
        tokens = member.split(None, 1)
        if (len(tokens) == 1):
            amount = 1
            key = member
        else:
            try:
                amount = int(tokens[0])
            except ValueError:
                raise KeggParseException("Non-specific reaction: " + s)
            key = tokens[1]
            
        if (key[0] != 'C'):
            raise KeggNonCompoundException("Compound ID doesn't start with C: " + key)
        try:
            cid = int(key[1:])
            compound_bag[cid] = compound_bag.get(cid, 0) + amount
        except ValueError:
            raise KeggParseException("Non-specific reaction: " + s)
    
    return compound_bag
  
def parse_reaction_formula(formula):
    """ parse a two-sided formula such as: 2 C00001 => C00002 + C00003 
        return the set of substrates, products and the direction of the reaction
    """
    tokens = re.split("([^=^<]+) (<*=>*) ([^=^>]+)", formula)
    if (len(tokens) != 5):
        return None
    left_bag = parse_reaction_formula_side(tokens[1])
    direction = tokens[2] # the direction: <=, => or <=>
    right_bag = parse_reaction_formula_side(tokens[3])
    
    return (left_bag, right_bag, direction)

def formula_to_atombag(formula):
    """
        Given a string representing a chemical formula, returns a bag containing the atoms in the formula
    """
    if (formula == "?" or formula.find("(") != -1 or formula.find(")") != -1):
        return {}

    atom_bag = {}
    for (atom, count) in re.findall("([A-Z][a-z]*)([0-9]*)", formula):
        if (count == ''):
            count = 1
        else:
            count = int(count)
        atom_bag[atom] = atom_bag.get(atom, 0) + count
    
    if ("R" in atom_bag): # this formula is not full ('R' is a wildcard not an atom)
        return {}
    
    return atom_bag

################################################################################
#                               EXCEPTIONS                                     #
################################################################################

class KeggParseException(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)
        
class KeggNonCompoundException(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)

class KeggReactionNotBalancedException(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)
    
class Reaction:
    def __init__(self, name, rid, sparse_reaction, weight=1):
        self.name = name
        self.rid = rid
        self.sparse_reaction = sparse_reaction
        self.weight = weight
        
    def get_cids(self):
        return set(self.sparse_reaction.keys())
    
    def unique_string(self):
        return " + ".join([("%d C%05d") % (coeff, cid) for (cid, coeff) in sorted(self.sparse_reaction.iteritems())])
    
    def is_futile(self):
        return max([abs(x) for x in self.sparse_reaction.values()]) > 0.01
    
    def is_balanced(self, cid2atom_bag):
        atom_diff = {}
        for (cid, coeff) in self.sparse_reaction.iteritems():
            atom_bag = cid2atom_bag.get(cid, {})
            if (atom_bag == None):
                # this reaction cannot be checked since there is an unspecific compound
                return False
            for (atomic_number, atom_count) in atom_bag.iteritems():
                atom_diff[atomic_number] = atom_diff.get(atomic_number, 0) + coeff * atom_count

        # ignore H and O inconsistencies
        atom_diff[1] = 0
        atom_diff[8] = 0
        
        return max([abs(x) for x in atom_diff.values()]) < 0.01
    
class Kegg:
    def __init__(self, log_file=None, max_carbons=None): # CO2, HCO3-
        util._mkdir('../kegg')
        self.COMPOUND_FILE = '../kegg/compound.txt.gz'
        self.REACTION_FILE = '../kegg/reaction.txt.gz'
        self.INCHI_FILE = '../kegg/inchi.txt.gz'
    
        self.edge_color = "cadetblue"
        self.edge_fontcolor = "indigo"
        self.edge_coeff_fontcolor = "darkolivegreen"
        self.node_fontcolor_cofactor = "dodgerblue"
        self.node_fontcolor_environment = "green"
        self.node_fontcolor = "white"
        self.node_fillcolor = "dodgerblue"
        self.font = "verdana"

        # cid2uid is a used for two purposes. One is to have canonical IDs for compounds
        # according to their INCHI labels (i.e. if two CIDs have the same INCHI, all occurences of the second
        # one will be changed to the first appearing CID). If a CID is not in the INCHI file, but is used
        # by one of the reactions, it might cause the reaction to be skipped. If compound formulas are to be
        # added using "database_updates.txt", the SETC line should appear at the top.
        self.cid2uid_map = {}
        self.compound2names_map = {}
        self.compound2atoms_map = {} # each value is a map of atoms to counts
        self.banned_compounds = set() # if one of these CIDs participates in a reaction, skip it
        self.banned_reactions = set() # skip all the of the reactions with these RIDs 
        self.cofactor_set = set()
        self.ignored_set = set()
        
        self.reactions = []
        
        if (log_file != None):
            self.LOG_FILE = log_file
        else:
            self.LOG_FILE = sys.stderr
            
    def prepare_all(self, UPDATE_FILE):
        self.parse_database()    
        self.update_database(UPDATE_FILE)
        self.preprocess_database()
    
    def parse_database(self):
        self.parse_compound_file()
        self.parse_inchi_file()
        self.parse_reaction_file()
   
    def parse_compound_file(self):
        compound_file = gzip.open(self.COMPOUND_FILE, 'r')
        for line in compound_file.readlines():
            if (line.find('%') != -1):
                line = line[0:line.find('%')]
            (compound_id, remainder) = line.rstrip().split(': ')
            cid = int(compound_id[1:])
            (name, formula) = remainder.split(' = ')
    
            atom_bag = formula_to_atombag(formula)
            if (len(atom_bag) > 0):
                self.compound2names_map[cid] = name.strip()
                self.compound2atoms_map[cid] = atom_bag
        compound_file.close()

    def formula_to_sparse(self, formula, ignore_protons=True):
        """
            translates a formula to a sparse-reaction
        """
        (left_bag, right_bag, direction) = parse_reaction_formula(formula)
        
        sparse_reaction = {}
        for (cid, count) in left_bag.iteritems():
            sparse_reaction[cid] = sparse_reaction.get(cid, 0) - count 
        for (cid, count) in right_bag.iteritems():
            sparse_reaction[cid] = sparse_reaction.get(cid, 0) + count 

        return (sparse_reaction, direction)
    
    def parse_inchi_file(self):
        """
            Parses the INCHI file and creates the cid2uid map for removing duplicates (CIDs of the same compound).
        """
        inchi_file = csv.reader(gzip.open(self.INCHI_FILE, 'r'), delimiter='\t')
        inchi2uid_map = {}
        for row in inchi_file:
            if (len(row) != 2):
                continue
            (cid, inchi) = row
            cid = int(cid[1:])
            if (not inchi in inchi2uid_map):
                inchi2uid_map[inchi] = cid
            self.cid2uid_map[cid] = inchi2uid_map[inchi]

    def parse_reaction_file(self, reaction_fname=None):
        """ Parse the formulas in reaction2formula_map and create a new map of reaction_id to pairs of
            substrates/products, taking the direction of each reaction into account.
            Note that two different reactions can have the same substrate-product pair, which means that
            there is redundancy in the database.
        """
        if (reaction_fname == None):
            reaction_fname = self.REACTION_FILE
        
        for cid in self.compound2names_map.keys():
            if (not cid in self.cid2uid_map):
                self.cid2uid_map[cid] = cid
        reaction_file = gzip.open(reaction_fname, 'r')
        for line in reaction_file.readlines():
            line = line.strip()
            if (line.find('%') != -1):
                line = line[0:line.find('%')]
            if (line == ""):
                continue
                
            if (len(line.split(':')) == 2):
                (reaction_id, formula) = line.split(':')
                rid = int(reaction_id.strip()[1:])
            else:
                raise Exception("Syntax Error - " + line)
    
            try:
                (sparse_reaction, direction) = self.formula_to_sparse(formula.strip())
                self.add_reaction(rid, direction, sparse_reaction)
            except KeggParseException:
                continue
            except KeggNonCompoundException:
                continue

    def update_database(self, fname):
        """
            Updates the database of reactions and compounds, using the update_database file.
            Commands are: SETR, NEWR, CARB, TRNS, DELR, SETC, NEWC, DELC, COFR
            
            In stiochio_mode, all reactions are not touched, so only SETC, NEWC, DELC, COFR are used.
        """
        
        update_file = open(fname, 'r')

        for line in update_file.readlines():
            if (line.find('%') != -1):
                line = line[0:line.find('%')]
            line = line.strip()
            if (line == ''):
                continue
            (command, rest) = line.split(' ', 1)
            line = rest.strip()
            
            if (command in ['SETR', 'NEWR', 'CARB', 'TRNS']):
                (reaction_id, formula) = line.split(':')
                rid = int(reaction_id.strip()[1:])
                try:
                    (sparse_reaction, direction) = self.formula_to_sparse(formula.strip())
                except KeggParseException:
                    continue
                self.add_reaction(rid, direction, sparse_reaction, name=reaction_id)
            elif (command == 'DELR'):
                rid = int(line[1:])
                self.banned_reactions.add(rid)
            elif (command in ['SETC', 'NEWC']):
                (cid, remainder) = line.rstrip().split(': ')
                cid = int(cid.strip()[1:])
                if (cid in self.cid2uid_map):
                    cid = self.cid2uid_map[cid]
                else:
                    self.cid2uid_map[cid] = cid # if this CID is not in the INCHI file
                (name, formula) = remainder.split(' = ')
                atom_bag = formula_to_atombag(formula)
                if (len(atom_bag) == 0):
                    raise KeggParseException("SETC must have a valid formula for the compound: " + line)
                self.compound2names_map[cid] = name
                self.compound2atoms_map[cid] = atom_bag
            elif (command == 'DELC'):
                cid = int(line[1:])
                if (cid in self.cid2uid_map):
                    cid = self.cid2uid_map[cid]
                else:
                    self.cid2uid_map[cid] = cid # if this CID is not in the INCHI file
                self.banned_compounds.add(cid)
            elif (command == 'COFR'): # cofactor
                (cid, name) = line.split('=')
                cid = int(cid.strip()[1:])
                if (cid in self.cid2uid_map):
                    cid = self.cid2uid_map[cid]
                else:
                    self.cid2uid_map[cid] = cid # if this CID is not in the INCHI file
                self.compound2names_map[cid] = name.strip()
                self.cofactor_set.add(cid)
            elif (command == 'SKIP'): # ignore this compound
                (cid, name) = line.split('=')
                cid = int(cid.strip()[1:])
                if (cid in self.cid2uid_map):
                    cid = self.cid2uid_map[cid]
                else:
                    self.cid2uid_map[cid] = cid # if this CID is not in the INCHI file
                self.compound2names_map[cid] = name.strip()
                self.ignored_set.add(cid)
            else:
                raise KeggParseException("Unknown command in Database Update file: " + command)
    
        update_file.close()

    def get_all_cids(self):
        cids = set(self.compound2atoms_map.keys()) - self.ignored_set
        return sorted(list(cids))

    def get_compound_name(self, cid):
        return self.compound2names_map.get(cid, None)
      
    def get_compound_carboncount(self, cid):
        # remove the -CoA atoms from the formulas.
        # since we are counting Carbons, but want to ignore co-factors, we have to do this.
        # fortunately, CoA and THF are the only carbon-containing cofactors.
        name = self.get_compound_name(cid)
        if (name == None):
            return None # this means this compound has been removed for some reason

        atom_bag = self.compound2atoms_map.get(cid, {})
        carbon_count = atom_bag.get('C', 0)
        
        if (len(atom_bag) == 0):
            return 0
        else:
            return carbon_count

    def preprocess_database(self):
        """ 
            create a new map of RID to reactants, without the co-factors.
            if the reaction is not balanced, skip it and don't add it to the new map.
        """
        
        unique_reactions = {}
        for r in self.reactions:
            if (r.rid in self.banned_reactions):
                self.LOG_FILE.write("This reaction has been banned: %s\n" % r.name)
            elif (not r.is_futile()):
                self.LOG_FILE.write("This reaction does nothing: %s\n" % r.name)
            elif (not r.is_balanced(self.compound2atoms_map)):
                self.LOG_FILE.write("This reaction isn't balanced: %s\n" % r.name)
            elif (len(self.banned_compounds.intersection(r.get_cids())) > 0):
                self.LOG_FILE.write("This reaction contains a banned compound: %s\n" % r.name)
            else:
                s = r.unique_string()
                if (s not in unique_reactions):
                    unique_reactions[s] = r

        self.reactions = unique_reactions.values()

    def add_reaction(self, rid, direction, sparse_reaction, name=None, weight=1):
        if (name == None):
            name = "R%05d" % rid
        if (direction in ["=>", "<=>"]):
            self.reactions.append(Reaction(name + "_F", rid, sparse_reaction, weight))
        if (direction in ["<=", "<=>"]):
            self.reactions.append(Reaction(name + "_R", rid, Kegg.reverse_sparse_reaction(sparse_reaction), weight))
    
    @staticmethod
    def reverse_sparse_reaction(sparse_reaction):
        backward_reaction = {}
        for (cid, coeff) in sparse_reaction.iteritems():
            backward_reaction[cid] = -coeff
        return backward_reaction
    
    def get_unique_cids_and_reactions(self):
        """
            Gather a set of all the CIDs (unique compound IDs) which are actually used.
            Remove reaction duplicates (i.e. have the same substrates and products,
            and store them in 'unique_reaction_map'.
        """
        cids = set()
        for r in self.reactions:
            cids = cids.union(r.get_cids())
        cids = cids.difference(self.ignored_set)
        cids = cids.difference(self.cofactor_set)
        cids = list(sorted(cids))

        Ncompounds = len(cids)
        Nreactions = len(self.reactions)
        sys.stderr.write("%d reactions with %d unique compounds\n" % (Nreactions, Ncompounds))
    
        # Create the columns, name the reactions (RID) in the stoichiometric matrix
        S = pylab.zeros((Ncompounds, Nreactions))
        rids = []
        f = []
        for j in range(Nreactions):
            r = self.reactions[j]
            rids.append(r.name)
            if (r.weight != 0):
                f.append((j, r.weight))
            
            for (cid, count) in r.sparse_reaction.iteritems():
                if (cid in cids):
                    i = cids.index(cid)
                    S[i, j] = count
                        
        return (f, S, cids, rids)
    
    def create_compound_node(self, Gdot, cid, name):
        node = self.get_node(Gdot, name)
        node.set_label('"%s"' % self.compound2names_map[cid])

        if (cid > 99999):
            node.set_tooltip('"C%05d"' % cid)
            node.set_URL('"http://www.genome.jp/Fig/compound/C%05d.gif"' % cid)
            node.set_fontcolor(self.node_fontcolor_environment) # color for cofactors
            node.set_shape("Mcircle")
            node.set_fontsize("12")
            node.set_fontname(self.font)
        elif (cid in self.cofactor_set):
            node.set_tooltip('"C%05d"' % cid)
            node.set_URL('"http://www.genome.jp/Fig/compound/C%05d.gif"' % cid)
            node.set_fontcolor(self.node_fontcolor_cofactor) # color for cofactors
            node.set_shape("none")
            node.set_fontsize("12")
            node.set_fontname(self.font)
        else:
            node.set_tooltip('"C%05d"' % cid)
            node.set_URL('"http://www.genome.jp/Fig/compound/C%05d.gif"' % cid)
            node.set_shape("box")
            node.set_style("filled")
            node.set_fontcolor(self.node_fontcolor) # color for non-cofcators
            node.set_fillcolor(self.node_fillcolor)
            node.set_fontsize("12")
            node.set_fontname(self.font)

        Gdot.add_node(node)
        return node
    
    def get_node(self, Gdot, name):
        node = Gdot.get_node(name)
        if (node != []):
            return node
        else:
            return pydot.Node(name)

    def create_reaction_nodes(self, Gdot, rid, flux=1):
        node_in = self.get_node(Gdot, "%s in" % rid)
        node_in.set_label("")
        node_in.set_shape("point")
        node_in.set_tooltip('"-> %s"' % rid)
        node_in.set_color(self.edge_color)
        Gdot.add_node(node_in)

        node_out = self.get_node(Gdot, "%s out" % rid)
        node_out.set_label("")
        node_out.set_shape("point")
        node_out.set_tooltip('"%s ->"' % rid)
        node_out.set_color(self.edge_color) # edge connector-point color
        Gdot.add_node(node_out)
        
        self.create_reaction_edge(Gdot, node_in, node_out, rid, flux=flux, arrowhead="none", arrowtail="none")

        return (node_in, node_out)

    def create_reaction_edge(self, Gdot, node_from, node_to, rid, flux=1, arrowhead="none", arrowtail="none"):
        """
            Create an edge for a reaction
        """
        if (node_from == None or node_to == None):
            return None
        edge = pydot.Edge(node_from, node_to)
        edge.set_color(self.edge_color) # edge line color
        edge.set_arrowhead(arrowhead)
        edge.set_arrowtail(arrowtail)
        edge.set_label('"%s x%.2f"' % (rid, flux))
        edge.set_fontcolor(self.edge_fontcolor) # edge label color
        edge.set_fontname(self.font)
        edge.set_fontsize("10")
        #edge.set_URL('"http://www.genome.jp/Fig/reaction/R%05d.gif"' % rid)
        Gdot.add_edge(edge)
        return edge
        
    def create_small_edge(self, Gdot, node_from, node_to, coeff=1, arrowhead="none", arrowtail="none"):
        """
            Create an edge that connects a compound to the 'point' node of a reaction (in or out)
        """
        edge = pydot.Edge(node_from, node_to)
        if (coeff != 1):
            edge.set_label('"%g"' % coeff)
        edge.set_color(self.edge_color)
        edge.set_fontcolor(self.edge_coeff_fontcolor)
        edge.set_arrowhead(arrowhead)
        edge.set_arrowtail(arrowtail)
        Gdot.add_edge(edge)
        return edge
    
    def draw_module(self, mid):
        (S, rids, cids) = self.get_module(mid)
        return self.draw_pathway(S, rids, cids)
            
    def draw_pathway(self, S, fluxes, cids, rids):
        Gdot = pydot.Dot()
        (Nr, Nc) = S.shape
        
        c_nodes = []
        for c in range(Nc):
            node_map = {} # a mapping of all the reactions that this compounds is participating in
            if cids[c] in self.ignored_set:
                for r in pylab.find(S[:,c] != 0): # this is an ignored compound, it shouldn't have any nodes
                    node_map[r] = None
            elif cids[c] in self.cofactor_set:
                for r in pylab.find(S[:,c] != 0): # this is a co-factor, create a new node for each reaction
                    node_map[r] = self.create_compound_node(Gdot, cids[c], str(cids[c]) + "_" + str(rids[r]))
            else:
                node = self.create_compound_node(Gdot, cids[c], str(cids[c]))
                for r in pylab.find(S[:,c] != 0): # point the node_map to the same node for every reaction
                    node_map[r] = node
            c_nodes.append(node_map)
       
        for r in range(Nr):
            rid = rids[r]
            if (rid == "R100000_F"):
                continue
            s_indices = pylab.find(S[r,:] < 0)
            p_indices = pylab.find(S[r,:] > 0)
            if (len(s_indices) == 1 and len(p_indices) == 1):
                c_s = s_indices[0]
                c_p = p_indices[0]
                if (S[r,c_s] == -1 and S[r,c_p] == 1):
                    self.create_reaction_edge(Gdot, c_nodes[c_s][r], c_nodes[c_p][r], rid=rid, flux=fluxes[r], arrowhead="open", arrowtail="none")
                    continue
            
            # this is not a simple 1-to-1 reaction
            (in_node, out_node) = self.create_reaction_nodes(Gdot, rid, flux=fluxes[r])
            for c in s_indices:
                self.create_small_edge(Gdot, c_nodes[c][r], in_node, coeff=-S[r,c], arrowhead="none")
            for c in p_indices:
                self.create_small_edge(Gdot, out_node, c_nodes[c][r], coeff=S[r,c], arrowhead="open")
        
        return Gdot