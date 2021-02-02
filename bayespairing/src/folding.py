import RNA
import math

def callback_structure_energy(s, data):
    if s:
        data['structures'].append(s)
        data['energies'].append(data['fc'].eval_structure(s))

class Fold:
    def __init__(self, s):
        self.md = RNA.md()
        self.md.uniq_ML = 1
        self.fc = RNA.fold_compound(s, self.md)
        self.kT = (RNA.K0+self.md.temperature)*RNA.GASCONST

    def non_redundant_sampling(self, nb, delta=None):
        """
        Boltzmann non-redundant sampling
        Return the results in an iterator of structure, energy

        If delta (kcal/mol) is given, return top nb structure with energy within delta from MFE
        """

        if delta:
            delta = int(round(delta,0))
            sub = self.fc.subopt(delta*100)
            if nb < len(sub):
                return map(lambda t: (t.structure,t.energy), sub[:nb])
            else:
                return map(lambda t: (t.structure,t.energy), sub)
        else:
            data = {
                'fc': self.fc,
                'structures': [],
                'energies': []
            }
            
            #ROMAN ADDED: rescale
            ss_mfe, mfe = self.fc.mfe()
            #print(mfe)
            self.fc.exp_params_rescale(mfe)
            self.fc.pf()
            
            self.fc.pbacktrack(nb, callback_structure_energy, data, RNA.PBACKTRACK_NON_REDUNDANT)
            return zip(data['structures'], data['energies'])

    def constraint_folding(self, c=""):
        """
        Folding with hard constraint and enforce constraint
        Return minimum free energy (mfe) structure, mfe, and free ensemble energy (fee)
        """
        if c:
            self.fc.hc_init()
            self.fc.hc_add_from_db(c, RNA.CONSTRAINT_DB_ENFORCE_BP | RNA.CONSTRAINT_DB_DEFAULT)
        ss_mfe, mfe = self.fc.mfe()
        _, fee = self.fc.pf()
        

        # Reset hard constraint
        #self.fc.hc_init()
        return ss_mfe, mfe, fee

    def compute_Boltzmann(self, e):
        """
        Return the value of Boltzmann function for a given energy
        B(e) = exp(-e/kT)
        """
        return math.exp(-e/self.kT)