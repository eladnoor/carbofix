import sys
import pylab
import util
from kegg import Kegg
import cplex, linprog
import csv
import gzip

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
        (f, S, b, cids, rids) = self.build_lp(kegg)

        # Find a solution with a minimal total flux
        solution = self.linprog(f, S, b, cids, rids, verbose=True)
        #(solution, min_flux) = self.solve(lp, verbose=True, milp=False)
        if (solution == None):
            sys.stderr.write("Couldn't find any cycle!\n")
        else:
            flux = sum([flux for (r, flux) in solution])
            
            #self.add_general_milp_constraints(lp)

            # Use iterative MILP to find the best suboptimal solutions
            min_flux = flux
            previous_solutions = []
            solution_file = open("../res/%s.txt" % experiment_name, "w")
            while (True):
                #print solution
                self.write_solution(solution_file, len(previous_solutions), S, solution, cids, rids)
                (sol_S, sol_fluxes, sol_cids, sol_rids) = self.get_solution_matrix(S, solution, cids, rids)
                sys.stderr.write("Found a minimal cycle with %d unique reactions and a total flux of %g\n" % (len(solution), sum(sol_fluxes)))
                
                #C = self.concentrate(sol_S, sol_cids)
                #sys.stderr.write(str(C))
                Gdot = kegg.draw_pathway(sol_S, sol_fluxes, sol_cids, sol_rids)
                Gdot.write_svg('../res/%s_%03d.svg' % (experiment_name, len(previous_solutions)), prog='dot')
                if (flux > min_flux * milp_factor):
                    break
                else:
                    # create the MILP problem to constrain the previous solutions not
                    # to reappear again.
                    previous_solutions.append(solution)
                    solution = self.linprog(f, S, b, cids, rids, previous_solutions=previous_solutions, verbose=True)
                    if (solution == None):
                        sys.stderr.write("Couldn't find any cycle!\n")
                        return
                    flux = sum([flux for (r, flux) in solution])

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

    def build_lp(self, kegg):
        (f, S, cids, rids) = kegg.get_unique_cids_and_reactions()
        
        if False: # change to True in order to export Stoichiometric data
            Ncompounds = S.shape[0]
            csv_output = csv.writer(gzip.open('../res/universal_matrix.csv.gz', 'w'))
            for c in xrange(Ncompounds):
                csv_output.writerow(['%g' % x for x in S[c,:]])

            rid_output = open('../res/reactions.txt', 'w')
            rid_output.write('\n'.join(rids))
            
            cid_output = open('../res/compounds.txt', 'w')
            cid_output.write('\n'.join(['C%05d' % c for c in cids]))
        
        # target (BIOMASS) has an export rate of 1.0 (AU), non-target compounds are at steady state (b = 0)
        b = [(cids.index(self.BIOMASS), 1.0)]
        return (f, S, b, cids, rids)

    def linprog(self, f, S, b, cids, rids, previous_solutions=[], verbose=False):
        """
            S is a NCxNR matrix where the rows represent compounds and the columns represent reactions.
            f is the goal vector, and is NRx1
            b is the linear constraint vector, and is NCx1
            
            cids is the list of CIDs from KEGG (NC long)
            reactions is a list of pairs of RID and direction (NR long)
        """
        (Ncompounds, Nreactions) = S.shape
        cpl = cplex.Cplex()
        if (verbose):
            cpl.set_log_stream(sys.stderr)
            cpl.set_results_stream(sys.stderr)
            cpl.set_warning_stream(sys.stderr)
        else:
            cpl.set_log_stream(None)
            cpl.set_results_stream(None)
            cpl.set_warning_stream(None)
        
        # reaction fluxes are the continuous variables
        cpl.set_problem_name('LP')
        for r in range(Nreactions):
            cpl.variables.add(lb=[0], ub=[self.UPPER_BOUND], names=[rids[r]])
        
        # add a linear constraint on the fluxes for each compound (mass balance)
        ind_constraints = 0
        for c in range(Ncompounds):
            cpl.linear_constraints.add(senses='E', names=["C%05d" % cids[c]])
            for r in pylab.find(S[c,:] != 0):
                cpl.linear_constraints.set_coefficients(ind_constraints, r, S[c,r])
            ind_constraints += 1

        for r in range(Nreactions):
            # add boolean indicator variables for each reaction, and minimize their sum
            cpl.variables.add(types='B', names=[rids[r] + "_gamma"])

            # add a constrain for each integer variable, so that they will be indicators of the flux variables
            # using v_i - M*gamma_i <= 0
            cpl.linear_constraints.add(senses='L', names=[rids[r] + "_bound"])
            cpl.linear_constraints.set_coefficients(ind_constraints, r, 1.0)
            cpl.linear_constraints.set_coefficients(ind_constraints, Nreactions+r, -self.UPPER_BOUND)
            b.append((ind_constraints, 0.0))
            ind_constraints += 1

        # map the weights of the goal function to the new boolean indicators
        f_bool = []
        for (r, weight) in f:
            f_bool.append((r+Nreactions, weight))

        if (previous_solutions != []):
            # for each previous solution, add constraints on the indicator variables,
            # so that that solution will not repeat (i.e. the sum over the previous reaction set must be
            # less than the size of that set).
            for i in range(len(previous_solutions)):
                cpl.linear_constraints.add(senses='L', names=["solution_%03d" % i])
                for (r, flux) in previous_solutions[i]:
                    cpl.linear_constraints.set_coefficients(ind_constraints, Nreactions+r, 1)
                b.append((ind_constraints, len(previous_solutions[i]) - 1))
                ind_constraints += 1
    
        cpl.objective.set_linear(f_bool)
        cpl.linear_constraints.set_rhs(b)
    
        cpl.write('../res/carbofix_%03d.lp' % len(previous_solutions), filetype='lp')
        cpl.solve()
        if (cpl.solution.get_status() == cplex.callbacks.SolveCallback.status.MIP_optimal):
            solution = pylab.array(cpl.solution.get_values())
            r_active = pylab.find(solution[Nreactions:] == 1)
            return [(r, solution[r]) for r in r_active] # a sparse vector of the solution fluxes
        elif (cpl.solution.get_status() == cplex.callbacks.SolveCallback.status.optimal):
            solution = pylab.array(cpl.solution.get_values())
            r_active = pylab.find(solution > 0)
            return [(r, solution[r]) for r in r_active] # a sparse vector of the solution fluxes
        else:
            return None
        
    def write_solution(self, handle, counter, S, solution, cids, rids):

        def write_compound_and_coeff(cid, coeff):
            if (coeff == 1):
                return "C%05d" % cid
            else:
                return "%d C%05d" % (coeff, cid)

        def write_reaction(handle, prefix, S, r, cids, rids):
            all_cids = [cids[c] for c in pylab.find(S[:,r] != 0)]
            if (self.ATMOSPHERE in all_cids or self.BIOMASS in all_cids):
                return None
            
            s = rids[r] + "  "
            left = []
            for c in pylab.find(S[:,r] < 0):
                left.append(write_compound_and_coeff(cids[c], -S[c,r]))
            right = []
            for c in pylab.find(S[:,r] > 0):
                right.append(write_compound_and_coeff(cids[c], S[c,r]))
            handle.write(prefix + rids[r] + "  " + " + ".join(left) + " => " + " + ".join(right) + "\n")
        
        handle.write('ENTRY       M-PATHWAY_%03d\n' % counter)
        handle.write('SKIP        FALSE\n')
        handle.write('NAME        M-PATHWAY_%03d\n' % counter)
        handle.write('TYPE        MARGIN\n')
        handle.write('CONDITIONS  pH=7.0,I=0.0,T=300\n')
        handle.write('C_MID       0.0001\n')
        for i in range(len(solution)):
            (r, flux) = solution[i]
            if (i == 0):
                write_reaction(handle, 'REACTION    ', S, r, cids, rids)
            else:
                write_reaction(handle, '            ', S, r, cids, rids)
        handle.write('///\n')
        handle.flush()
    
    def get_solution_matrix(self, S, solution, cids, rids):
        sol_indices = [r for (r, flux) in solution]
        sol_fluxes = [flux for (r, flux) in solution]
        sol_rids = [rids[r] for (r, flux) in solution]
        S_interm = S[:, sol_indices].T
        nonzero_compounds = pylab.find(pylab.sum(S_interm != 0, 0))
        sol_cids = [cids[c] for c in nonzero_compounds]
        sol_S = S_interm[:, nonzero_compounds]
        return (sol_S, sol_fluxes, sol_cids, sol_rids)

    def concentrate(self, S, cids):
        cid2dG0_f = {} # @@@ use the GIBBS program to fill in this dictionary

        dG0_f = pylab.matrix([cid2dG0_f.get(cid, 0) for cid in cids]).T
        R = 8.31 # gas constant (J/K mol)
        T = 300 # temperature (K)
        f = pylab.ones((S.shape[1], 1))
        dG0_r = pylab.dot(S, dG0_f)
        ub = pylab.ones((S.shape[1], 1)) * 1e-0
        lb = pylab.ones((S.shape[1], 1)) * 1e-10
        
        C = linprog.linprog(f=f, A=S, b=dG0_r/(-R*T), lb=lb, ub=ub)
        return C

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
    
    #cf.find_path('CO2 to GLYOXYLATE', sources=[11, 288], targets=[48], milp_factor=20)
    cf.find_path('CO2 to GA3P', sources=[9, 11, 288], targets=[118], milp_factor=20)
    
####################################################################################################

if (__name__ == '__main__'):
    main()
