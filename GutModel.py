from pysb import *
from pysb.simulator import ScipyOdeSimulator
import numpy as np
import matplotlib.pyplot as plt

Model()

n_levels = 10
Monomer('Bacteroides', ['energy'], {'energy': ['_%d' % i for i in range(n_levels)]})
Monomer('Clostridium', ['energy'], {'energy': ['_%d' % i for i in range(n_levels)]})
Monomer('Bifidobacterium', ['energy'], {'energy': ['_%d' % i for i in range(n_levels)]})
Monomer('Desulfobrivio', ['energy'], {'energy': ['_%d' % i for i in range(n_levels)]})

Monomer('Inulin')
Monomer('Glucose')
Monomer('Lactose')
Monomer('Fructo')
Monomer('ChondSulf')

Parameter('k_Bact_Inulin', 1)

[Rule('BactEatInulin_%d_%d' % (i,i+1), Bacteroides(energy='_%d' % i) + Inulin() >> Bacteroides(energy= '_%d' % (i+1)), k_Bact_Inulin)
 for i in range(n_levels - 1)]

print(model.rules)


