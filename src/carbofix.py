import sys
import numpy as np
from scipy.sparse import find
import util
from kegg import Kegg
import csv
import gzip
from pulp import LpProblem, LpMaximize, LpMinimize, LpVariable, LpAffineExpression,\
                 solvers, LpContinuous, LpBinary, LpStatusOptimal, lpSum, LpStatus
                 
################################################################################
#                               CONSTANTS & DEFAULTS                           #
################################################################################

class Carbofix:
    def __init__(self, default_cost=None, update_file='../rec/database_updates.txt'):
        util._mkdir('../res')
        self.LOG_FILE = open('../res/carbocyc.log', 'w')
        self.UPDATE_FILE = update_file
        self.ATMOSPHERE = 100000
        self.BIOMASS = 100001
        self.UPPER_BOUND = 100.0

    def __del__(self):
        self.LOG_FILE.close()
        
    def find_path(self, experiment_name, sources, targets, milp_factor=0):
        kegg = Kegg(self.LOG_FILE)
        kegg.prepare_all(self.UPDATE_FILE)
        self.add_import_and_export(kegg, sources, targets)
        f, S, cids, rids = self.build_lp(kegg)

        # Find a solution with a minimal total flux
        solution = self.linprog(f, S, cids, rids, verbose=True)
        #(solution, min_flux) = self.solve(lp, verbose=True, milp=False)
        if solution is None:
            sys.stderr.write("Couldn't find any cycle\n")
        else:
            flux = sum([flux for (r, flux) in solution])
            
            #self.add_general_milp_constraints(lp)

            # Use iterative MILP to find the best suboptimal solutions
            min_flux = flux
            previous_solutions = []
            with open("../res/%s.txt" % experiment_name, "w") as solution_file:
                while True:
                    self.write_solution(solution_file, len(previous_solutions), S,
                                        solution, cids, rids)
                                        
                    sol_S, sol_fluxes, sol_cids, sol_rids = \
                        self.get_solution_matrix(S, solution, cids, rids)
                        
                    sys.stderr.write("Found a minimal cycle with %d unique reactions and a total flux of %g\n" % (len(solution), sum(sol_fluxes)))
                    
                    Gdot = kegg.draw_pathway(sol_S, sol_fluxes, sol_cids, sol_rids)
                    Gdot.write_svg('../res/%s_%03d.svg' %
                                   (experiment_name, len(previous_solutions)),
                                   prog='dot')
                    if flux > min_flux * milp_factor:
                        return
                    else:
                        # create the MILP problem to constrain the previous solutions not
                        # to reappear again.
                        previous_solutions.append(solution)
                        solution = self.linprog(f, S, cids, rids,
                                                previous_solutions=previous_solutions,
                                                verbose=True)
                        if solution is None:
                            sys.stderr.write("Couldn't find any more cycles\n")
                            return
                        flux = sum([_f for (r, _f) in solution])

    def add_import_and_export(self, kegg, sources, targets):
        kegg.compound2names_map[self.BIOMASS] = 'biomass'
        kegg.compound2names_map[self.ATMOSPHERE] = 'air'
        kegg.compound2atoms_map[self.BIOMASS] = {}
        kegg.compound2atoms_map[self.ATMOSPHERE] = {}
        
        # Add input and output reactions
        kegg.add_reaction(100000, '=>', {self.ATMOSPHERE : 1}, name='IMPORT', weight=0) # atmosphere is "free"
        for cid in sources:
            # input any of the sources for the atmosphere
            kegg.add_reaction(100000 + cid, '=>', {self.ATMOSPHERE : -1, cid : 1}, name='IMPORT_C%05d' % cid, weight=0)
        
        for cid in targets:
            kegg.add_reaction(200000 + cid, '=>', {self.BIOMASS : 1, cid : -1}, name='EXPORT_C%05d' % cid, weight=0)

    def build_lp(self, kegg, debug_mode=False):
        f, S, cids, rids = kegg.get_unique_cids_and_reactions()

        # convert the KEGG rids to unique ones by adding a counter        
        unique_rids = ['%s_%05d' % (rid, r) for (r, rid) in enumerate(rids)]

        if debug_mode:
            Ncompounds = S.shape[0]
            csv_output = csv.writer(gzip.open('../res/universal_matrix.csv.gz', 'w'))
            for c in xrange(Ncompounds):
                csv_output.writerow(['%g' % x for x in S[c,:]])

            rid_output = open('../res/reactions.txt', 'w')
            rid_output.write('\n'.join(rids))
            
            cid_output = open('../res/compounds.txt', 'w')
            cid_output.write('\n'.join(['C%05d' % c for c in cids]))
        
        return f, S, cids, unique_rids

    def linprog(self, f, S, cids, rids, previous_solutions=[], verbose=False):
        """
            S is a NCxNR matrix where the rows represent compounds and the columns represent reactions.
            f is the goal vector, and is NRx1
            b is the linear constraint vector, and is NCx1
            
            cids is the list of CIDs from KEGG (NC long)
            reactions is a list of pairs of RID and direction (NR long)
        """
        Ncompounds, Nreactions = S.shape
        prob = LpProblem('CarboFix', sense=LpMinimize)
        prob.solver = solvers.CPLEX(msg=verbose)
        
        # reaction fluxes are the continuous variables
        flux_dict = LpVariable.dicts('v', rids, lowBound=0, upBound=self.UPPER_BOUND)
        
        # each row in the stoichiometric matrix becomes an affine constraint
        # (a row in S*v = 0)
        for c, cid in enumerate(cids):
            rids_changing_c, coeffs = find(S[c,:])[1:]
            fluxes_changing_c = map(lambda r: flux_dict[rids[r]], rids_changing_c)
            S_times_v = LpAffineExpression(zip(fluxes_changing_c, coeffs), name='S_times_v')
            if cid == self.BIOMASS:
                prob.addConstraint(S_times_v == 1, 'mass_balance_%s' % cid)
            else:
                prob.addConstraint(S_times_v == 0, 'mass_balance_%s' % cid)

        # add boolean indicator variables for each reaction, and minimize their sum
        gamma_dict = LpVariable.dicts('gamma', rids, lowBound=0, upBound=1, cat='Integer')
        
        # Make each gamma_r into a flux indicator
        # (i.e. so that if v_r > 0 then gamma_i must be equal to 1).
        # We use the following constraint:
        #          v_r <= M*gamma_i
        for rid in rids:
            prob.addConstraint(flux_dict[rid] <= self.UPPER_BOUND * gamma_dict[rid],
                               'gamma_bound_%s' % rid)

        # map the weights of the goal function to the new boolean indicators
        gamma_list = map(gamma_dict.get, rids)
        objective = LpAffineExpression(zip(gamma_list, f.flat),
                                       name='f_times_gamma')
        prob.setObjective(objective)

        if previous_solutions != []:
            # for each previous solution, add constraints on the indicator variables,
            # so that that solution will not repeat (i.e. the sum over the previous reaction set must be
            # less than the size of that set).
            for i, prevsol in enumerate(previous_solutions):
                prev_gammas = map(lambda (r, flux): gamma_dict[rids[r]], prevsol)
                prob.addConstraint(lpSum(prev_gammas) <= len(prev_gammas) - 1,
                                   'eliminate_solution_%03d' % i) 
    
        prob.writeLP('../res/carbofix_%03d.lp' % len(previous_solutions))
        prob.solve()
        if prob.status != LpStatusOptimal:
            return
        
        # make a sparse representation of the solution flux vector
        solution = []
        for r, rid in enumerate(rids):
            if gamma_dict[rid].varValue == 1:
                solution.append((r, flux_dict[rid].varValue))
        
        return solution
        
    def write_solution(self, handle, counter, S, solution, cids, rids):

        def write_compound_and_coeff(cid, coeff):
            if coeff == 1:
                return "C%05d" % cid
            else:
                return "%d C%05d" % (coeff, cid)

        def write_reaction(handle, prefix, S, r, cids, rids):
            all_cids = [cids[c] for c in find(S[:,r])[1]]
            if (self.ATMOSPHERE in all_cids or self.BIOMASS in all_cids):
                return None
            
            left = []
            for c in find(S[:,r] < 0)[0]:
                left.append(write_compound_and_coeff(cids[c], -S[c,r]))
            right = []
            for c in find(S[:,r] > 0)[0]:
                right.append(write_compound_and_coeff(cids[c], S[c,r]))
            handle.write(prefix + rids[r] + "  " + " + ".join(left) + " => " + " + ".join(right) + "\n")
        
        handle.write('ENTRY       M-PATHWAY_%03d\n' % counter)
        handle.write('SKIP        FALSE\n')
        handle.write('NAME        M-PATHWAY_%03d\n' % counter)
        handle.write('TYPE        MARGIN\n')
        handle.write('CONDITIONS  pH=7.0,I=0.0,T=300\n')
        handle.write('C_MID       0.0001\n')
        for i in range(len(solution)):
            r, flux = solution[i]
            if (i == 0):
                write_reaction(handle, 'REACTION    ', S, r, cids, rids)
            else:
                write_reaction(handle, '            ', S, r, cids, rids)
        handle.write('///\n')
        handle.flush()
    
    def get_solution_matrix(self, S, solution, cids, rids):
        sol_indices, sol_fluxes = zip(*solution)
        sol_rids = map(lambda r: rids[r], sol_indices)
        S_interm = S[:, sol_indices].T
        nonzero_compounds = find(np.sum(S_interm != 0, 0))[1]
        sol_cids = [cids[c] for c in nonzero_compounds]
        sol_S = S_interm[:, nonzero_compounds]
        return sol_S, sol_fluxes, sol_cids, sol_rids

################################################################################
#                               MAIN                                           #
################################################################################

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

def main():
    cf = Carbofix()
    
    cf.find_path('CO2 to GLYOXYLATE', sources=[11, 288], targets=[48], milp_factor=20)
    #cf.find_path('CO2 to GA3P', sources=[9, 11, 288], targets=[118], milp_factor=20)
    
####################################################################################################

if (__name__ == '__main__'):
    main()
