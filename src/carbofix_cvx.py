#!/usr/bin/python2.5
###############################################################################
# Except for the default python packages for python 2.5, this program requires:
# pydot - a python interface for graphviz (which requires the libraries graphviz
# and networkx)
# glpk - Gnu Linear Programming Kit
# OpenBabel - chemical format conversions
################################################################################

import os
import sys
import pydot
import getopt
import urllib
import types
from html_writer import HtmlWriter
import util
import re
import math
from copy import deepcopy
from bag import bag
from kegg import Kegg
import glpk
import ctypes
import csv
import sqlite3
import pylab
import cvxopt
import cvxopt.lapack
import cvxmod

################################################################################
#                               CONSTANTS & DEFAULTS                           #
################################################################################

def cvxopt_norm1(A, b):
    """
        Solve A*x = b; minimize |x|
        
        where |x| is the NORM1 of x

        A and b, should both be Numpy matrices.
    """
    # First, A might be singular, therefore we use SVD to transform this problem
    # to one with a non-singular matrix

    # If A is MxN, and or rank R:
    # SVD returns: U - MxM unitary matrix, s - R singular values (and zeros), V_H - NxN unitary matrix
    (Nc, Nr) = A.size
    Ns = min(Nc,Nr)
    S_temp = cvxopt.matrix(0, (Ns, 1), tc='d')
    U = cvxopt.spmatrix([], [], [], (Nc, Nr))
    Vt = cvxopt.spmatrix([], [], [], (Nc, Nr))
    sys.stderr.write("Running SVD on full stoichiometric matrix ... ")
    cvxopt.lapack.gesvd(cvxopt.matrix(A), S=S_temp, U=U, Vt=Vt)
    sys.stderr.write("[DONE]\n")
    
    # We need convert s to a MxN matrix (S), so that we get: A = U*S*V_H
    S = cvxopt.spmatrix(S_temp, range(Ns), range(Ns), (Nc, Nr))

    # We've transformed the constraint space from Ax = b to U*(S*V_H*x) = b, using
    # SVD. We use U as our constraint matrix, and we solve for:
    #   z = S*VH*x

    # Optimization variable
    z2 = cvxmod.optvar(name="singular values * vh * flux_vector", rows=Nr, cols=1)
    U2 = cvxopt.matrix(U)
    b2 = cvxopt.matrix(b)
    # Optimization problem
    p = problem(minimize(cvxmod.atoms.norm1(z)), constr = [ U2*z2 == b2 ])
    p.solve()

    z = pylab.matrix(value(z2))
    x = (S*V_H).I * z

    return x

def linprog(f, A, b, lb, ub, verbose=False):
    """
        The constraints are:
            A*x = b
            lb <= x <= ub
        
        The optimization rule is:
            minimize f*x
        
        In matrix A:
            rows are reactions (indexed by 'r') - total of Nr
            columns are compounds (indexed by 'c') - total of Nc
    """
    Nr = A.shape[0]
    Nc = A.shape[1]

    lp = glpk.glp_create_prob()
    glpk.glp_set_prob_name(lp, "bottlenecks")
    glpk.glp_set_obj_dir(lp, glpk.GLP_MIN)
    glpk.glp_add_rows(lp, Nr)
    
    for r in range(Nr):
        glpk.glp_set_row_name(lp, r+1, "Reaction %d" % r)
        glpk.glp_set_row_bnds(lp, r+1, glpk.GLP_UP, 0.0, b[r]) # 0.0 is ignored since the lower bound is -infty

    # Create the columns, name the reactions (RID) and add the values to a temporary sparse matrix
    glpk.glp_add_cols(lp, Nc)
    for c in range(Nc):
        glpk.glp_set_col_name(lp, c+1, "Compound %d" % c)
        glpk.glp_set_col_bnds(lp, c+1, glpk.GLP_DB, lb[c], ub[c])

    # Copy the sparse matrix to C-type structures that GLPK can later use
    size = Nr*Nc 
    ia = (ctypes.c_int    * (1+size))() # array of row indices in sparse matrix
    ja = (ctypes.c_int    * (1+size))() # array of col indices in sparse matrix
    ar = (ctypes.c_double * (1+size))() # array of values in sparse matrix
    for r in range(Nr):
        for c in range(Nc):
            ia[1+r*Nc+c] = ctypes.c_int(r+1)
            ja[1+r*Nc+c] = ctypes.c_int(c+1)
            ar[1+r*Nc+c] = ctypes.c_double(A[r,c])

    glpk.glp_load_matrix(lp, size, ia, ja, ar)
    
    for c in range(Nc):
        glpk.glp_set_obj_coef(lp, c+1, f[c])

    parm = glpk.glp_smcp()
    glpk.glp_init_smcp(parm)
    if (verbose):
        parm.msg_lev = glpk.GLP_MSG_ON
    else:
        parm.msg_lev = glpk.GLP_MSG_OFF
    retval = glpk.glp_simplex(lp, parm)
    if (retval != 0):
        return None
    
    f_min = glpk.glp_get_obj_val(lp)

    solution = pylab.zeros((Nc, 1))
    for c in range(Nc):
        solution[c,0] = glpk.glp_get_col_prim(lp, c+1)
    
    m = pylab.dot(A,solution) - b
    for i in range(Nr):
        if (m[i,0] > 0):
            return None
    return solution

class Carbofix:
    def __init__(self, default_cost=None, update_file='../rec/database_updates.txt'):
        util._mkdir('../res')
        self.LOG_FILE = open('../res/carbocyc.log', 'w')
        self.UPDATE_FILE = update_file
        self.ATMOSPHERE = 'C90002'
        self.BIOMASS = 'C90001'
        self.UPPER_BOUND = 100.0
        self.EPSILON_FLUX = 0.000001
        self.rid2cost = self.read_cost_from_database()
        self.S = None # Universal Stoichiometric Matrix
        if (default_cost == None):
            self.default_cost = util.median(self.rid2cost.values())
        else:
            self.default_cost = default_cost
        
    def read_cost_from_database(self):
        try:
            comm = sqlite3.connect('../../enzyme_db/rec/enzymes.sqlite')
            cursor = comm.cursor()
            cursor.execute('SELECT rid, side, MAX(tn) FROM merged_final WHERE tn IS NOT NULL GROUP BY rid, side;')
            rid2cost = {}
            for row in cursor:
                (rid, side, cost) = row
                key = (int(rid), int(side))
                rid2cost[key] = float(cost)
            comm.close()
            return rid2cost
        except sqlite3.OperationalError:
            return {}

    def __del__(self):
        self.LOG_FILE.close()
        
    def get_cost(self, rid, direction): # Use the specific activity, not k_cat !!!
        if (direction == "<="):
            side = 1
        elif (direction == "=>"):
            side = -1
        else:
            raise Exception("unknown direction - " + direction)

        if (rid[0] != 'R'):
            return 0

        rid = int(rid[1:])
        if ((rid, side) in self.rid2cost):
            return 1.0 / self.rid2cost[(rid, side)]
        elif (self.default_cost != None):
            return 1.0 / self.default_cost
        else:
            return 1.0
    
    def build_problem(self, sources, targets, use_cost=False):
        reaction_list = self.KEGG.get_all_reactions()
    
        # Gather a set of all the CIDs (unique compound IDs) which are actually used.
        # Remove reaction duplicates (i.e. have the same substrates and products,
        # and store them in 'unique_reaction_map'.
        all_cid_set = set(sources)
        unique_reaction_map = {}
        for (subs, prod, rid, direction) in reaction_list:
            all_cid_set = all_cid_set.union(subs.to_set())
            all_cid_set = all_cid_set.union(prod.to_set())
            
            if (rid[0] != 'R'):
                cost = 0.0
            elif (use_cost):
                cost = self.get_cost(rid, direction)
            else:
                cost = 1.0

            unique_str = self.KEGG.cidbag_to_keggstr(subs) + " => " + self.KEGG.cidbag_to_keggstr(prod)
            if ((not unique_str in unique_reaction_map) or unique_reaction_map[unique_str][4] > cost):
                unique_reaction_map[unique_str] = (subs, prod, rid, direction, cost)
        all_cid_set.add(self.BIOMASS)
        
        # Create an index for each CID, and a map to convert CIDs to their index (row) in the matrix
        self.cid_list = list(sorted(all_cid_set))
        self.cids2index = {}
        for i in range(len(self.cid_list)):
            self.cids2index[self.cid_list[i]] = i
    
        # Translate each reaction to a sparse vector of stoichiometric coefficients,
        # according to the cids2index map.
        self.sparse_columns = []
        self.rid_list = []
        cost_weights = {}
        counter = 0
        for (subs, prod, rid, direction, cost) in unique_reaction_map.itervalues():
            svector = []
            
            if (len(subs.intersection(prod)) > 0):
                continue
            
            for (cid, count) in subs.iteritems():
                svector.append((self.cids2index[cid], -count))
            for (cid, count) in prod.iteritems():
                svector.append((self.cids2index[cid], count))
    
            self.rid_list.append(rid + ":" + direction)
            self.sparse_columns.append(svector)

            cost_weights[counter] = cost
            counter += 1
                
        # Add unbalanced reactions for input and output reactions
        for cid in sources:
            self.sparse_columns.append([(self.cids2index[cid], 1)]) # allow unbalanced influx of each carbon source
            self.rid_list.append('IN_' + cid + ':=>')
            self.KEGG.add_reaction('IN_' + cid, '=>', bag([self.ATMOSPHERE]), bag([cid]))
        for cid in targets:
            self.sparse_columns.append([(self.cids2index[cid], -1), (self.cids2index[self.BIOMASS], 1)]) # allow any target to convert to BIOMASS
            self.rid_list.append('EX_' + cid + ':=>')
            self.KEGG.add_reaction('EX_' + cid, '=>', bag([cid]), bag([self.BIOMASS]))
        
        sys.stderr.write("%d reactions with %d unique compounds\n" % (len(self.rid_list), len(self.cid_list)))
    
        Nr = len(self.cid_list) # each row is a compound
        Nc = len(self.rid_list) # each column is a reaction

        # convert the 'sparce_columns' to a full stoichimetric matrix S
        x = []
        I = []
        J = []
        for c in range(len(self.sparse_columns)):
            for (r, value) in self.sparse_columns[c]:
                x.append(value)
                I.append(r)
                J.append(c)
        self.S = cvxopt.spmatrix(x, I, J, (Nr, Nc))
                
        self.b = cvxopt.matrix(0.0, (Nr, 1))
        i_biomass = self.cid_list.index(self.BIOMASS)
        self.b[i_biomass] = 1.0 # target (BIOMASS) has an export rate of 1.0 (AU)
    
    def initialize(self, sources, targets, use_cost=False):
        self.KEGG = Kegg(self.LOG_FILE, sources=sources)
        self.KEGG.prepare_all(self.UPDATE_FILE)
    
        self.KEGG.compound2names_map[self.BIOMASS] = 'biomass'
        self.KEGG.compound2names_map[self.ATMOSPHERE] = 'air'
        self.KEGG.compound2atoms_map[self.BIOMASS] = bag()
        self.KEGG.compound2atoms_map[self.ATMOSPHERE] = bag()
        
        self.KEGG.preprocess_database()
        self.build_problem(sources=sources, targets=targets, use_cost=use_cost)

    def find_path(self, experiment_name, sources, targets, use_cost=False, milp_factor=0):
        self.initialize(sources=sources, targets=targets, use_cost=use_cost)

        # Find a solution with a minimal total flux
        (solution, min_flux) = self.solve()
        Gdot = self.solution_to_Gdot(solution, use_cost)
        self.write_solution(Gdot, experiment_name, 0)
        return

        if (solution == None):
            sys.stderr.write("Couldn't find any cycle!\n")
        else:
            sys.stderr.write("Found a minimal cycle a total flux of %g:\n" % min_flux)
            solution_counter = 0

            # Use iterative MILP to find the best suboptimal solutions
            flux = min_flux
            while (True):
                print solution
                C = self.concentrate(solution)
                print C
                Gdot = self.solution_to_Gdot(solution, use_cost)
                self.write_solution(Gdot, experiment_name, solution_counter)
                solution_counter += 1
                if (flux > min_flux * milp_factor):
                    break
                
                self.add_milp_constraint(solution)
                (solution, flux) = self.solve()
                if (solution == None):
                    break
    
    def solve(self):
        """
            Use the CVXOPT solver to solve the problem.
            The return value is the vector of fluxes.
        """

        sys.stderr.write("Writing the (%dx%d) universal stoichiometric matrix to file and the 'b' ... " % (self.S.size[0], self.S.size[1]))
        csv_file = csv.writer(open('../res/USM.csv', 'w'))
        for i in range(len(self.S)):
            csv_file.writerow((self.S.I[i], self.S.J[i], self.S.V[i]))
        del csv_file
        
        csv_file = csv.writer(open('../res/b.csv', 'w'))
        for i in range(len(self.b)):
            if (self.b[i,0] != 0):
                csv_file.writerow((i, self.b[i,0]))
        del csv_file
        
        sys.stderr.write("[DONE]\n")
        
        x = cvxopt_norm1(self.S, self.b)

        solution = {}
        for j in range(x.shape[0]):
            if (x[j] > self.EPSILON_FLUX):
                solution[j] = x[j]
            
        return (solution, sum(x))
    
    def add_general_milp_constraints(self, lp):
        """
            add the binary indicator variables
        """
        Nc = glpk.glp_get_num_cols(lp)
        Nr = glpk.glp_get_num_rows(lp)
        
        glpk.glp_add_cols(lp, Nc)
        for j in range(Nc):
            glpk.glp_set_col_kind(lp, Nc+j+1, glpk.GLP_BV) # BV - Binary Value

        # add the constraints that cause each indicator to be 0 if its corresponding
        # flux variable is 0 and 1 if the flux is positive.
        glpk.glp_add_rows(lp, Nc)
        for j in range(Nc):
            self.set_mat_row(lp, Nr+j, [j, Nc+j], [1, -self.UPPER_BOUND])
            glpk.glp_set_row_bnds(lp, Nr+j+1, glpk.GLP_UP, -self.UPPER_BOUND, 0.0) # the lower bound is ignored

    def add_milp_constraint(self, lp, solution):
        """
            add the constraints that prevent previous solutions from being chosen again
        """
        Nr = glpk.glp_get_num_rows(lp)
        glpk.glp_add_rows(lp, 1)
        ind = solution.keys()
        val = [1.0] * len(ind)
        self.set_mat_row(lp, Nr, ind, val)
        glpk.glp_set_row_bnds(lp, Nr+1, glpk.GLP_UP, 0.0, len(ind)-1.0) # the lower bound is ignored

    def solution_to_Gdot(self, solution, use_cost):
        path = []
        weights = []
        path_str = []
        for j in solution.keys():
            path.append(self.rid_list[j])
            weights.append(solution[j])
            (rid, direction) = self.rid_list[j].split(':')
            cost = self.get_cost(rid, direction)
            if (use_cost and cost != None and cost != 0):
                path_str.append("\"x%g [%.3f]\"" % (solution[j], cost))
            else:
                path_str.append("\"x%g\"" % (solution[j]))
    
        Gdot = self.KEGG.build_graph(path, weights, path_str)
        
        atmosphere_node = self.KEGG.get_node(Gdot, self.ATMOSPHERE)
        atmosphere_node.set_color('green')
        atmosphere_node.set_shape('Mcircle')
    
        biomass_node = self.KEGG.get_node(Gdot, self.BIOMASS)
        biomass_node.set_color('magenta')
        biomass_node.set_shape('Mcircle')
        
        for id in self.KEGG.cofactor_set:
            cofactor_node = self.KEGG.get_node(Gdot, id, create_if_missing=False)
            if (cofactor_node != None):
                cofactor_node.set_color('black')
                cofactor_node.set_shape('none')           
        
        return Gdot
    
    def concentrate(self, solution):
        smatrix = []
        cid2index = {}
        cid_counter = 0
        dG0_f = []
        cid2dG0_f = {} # @@@ use the GIBBS program to fill in this dictionary

        reactions = solution.keys()
        for i in range(len(reactions)):
            for (row_index, stoich_coeff) in self.sparse_columns[reactions[i]]:
                cid = self.cid_list[row_index]
                if (cid not in cid2index):
                    cid2index[cid] = cid_counter
                    cid_counter += 1
                    dG0_f.append(cid2dG0_f.get(cid, 0.0))
                j = cid2index[cid]
                smatrix.append((i, j, stoich_coeff))

        S = pylab.zeros((len(reactions), len(cid2index)))
        for (i, j, coeff) in smatrix:
            S[i, j] = coeff
        
        dG0_f = pylab.matrix(dG0_f).T
        R = 8.31 # gas constant (J/K mol)
        T = 300 # temperature (K)
        f = pylab.ones((S.shape[1], 1))
        dG0_r = pylab.dot(S, dG0_f)
        ub = pylab.ones((S.shape[1], 1)) * 1e-2
        lb = pylab.ones((S.shape[1], 1)) * 1e-6
        
        C = linprog(f, S, dG0_r/(-R*T), lb, ub)
        return C
    
    def write_solution(self, Gdot, experiment_name, solution_index):
        svg_filename = '../res/%s_%02d.svg' % (experiment_name, solution_index)
        sys.stderr.write("Writing to SVG file '%s' ... " % svg_filename)
        Gdot.write_svg(svg_filename, prog='dot')
        sys.stderr.write("[DONE]\n")

################################################################################
#                               MAIN                                           #
################################################################################

def main():
    cf = Carbofix()
    
    # LEGEND:
    
    # C00118 - D-glyceraldehyde 3-phosphate
    # C00197 - 3-phospho-glycerate
    # C00048 - glyoxylate
    # C00036 - oxaloacetate
    # C00024 - acetyl-CoA
    # C00469 - ethanol
    # C00132 - methanol
    # C01438 - methane
    # C06142 - butanol
    
    # R00024 - Rubisco (ribulose-1,5P => 2 glycerate-3P)
    # R00134 - Formate oxidoreductase (CO2 => Formate)
    # R00345 - PEP carboxylase (oxaloacetate => phosphoenolpyruvate)
    # R00344 - pyruvate carboxylase (pyruvate => oxaloacetate)
    # R00742 - Acetyl-CoA carboxylase (acetyl-CoA => malonyl-CoA)
    # R01859 - Propionyl-CoA carboxylase (Propanoyl-CoA => (S)-Methylmalonyl-CoA)
    # R01197 - 2-Ketoglutarate Synthase (succinyl-CoA => 2-ketoglutarate)
    # R00709 - Isocitrate Dehydrogenase (2-ketoglutarate => isocitrate)
    # R00210 - pyruvate:NADP+ 2-oxidoreductase (acetyl-CoA => pyruvate)
    # R00353 - Malonyl-CoA:pyruvate carboxyltransferase (Malonyl-CoA + Pyruvate <=> Acetyl-CoA + Oxaloacetate)


    cf.find_path('CO2 to GLYOXYLATE', sources=['C00011', 'C00288'], targets=['C00048'], use_cost=False, milp_factor=0)
    
    #s = eval('{10535: 2.0, 10537: 1.0, 5132: 1.0, 8950: 1.0, 1751: 1.0, 2648: 1.0}')
    #cf.initialize(sources=['C00011', 'C00288'], targets=['C00048'])
    #print cf.concentrate(s)
    
    #cf.find_path('Glucose to CO2', sources=['C00031'], targets=['C00011'], use_cost=False, milp_factor=0)

####################################################################################################

if (__name__ == '__main__'):
    main()
