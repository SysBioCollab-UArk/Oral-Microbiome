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
Monomer('Lactate')

Parameter('k_Bact_Inulin', 1)

[Rule('BactEatInulin_%d_%d' % (i,i+1), Bacteroides(energy='_%d' % i) + Inulin() >> Bacteroides(energy= '_%d' % (i+1)), k_Bact_Inulin)
 for i in range(n_levels - 1)]

Parameter('k_Clost_Inulin', 1)

[Rule('ClostEatInulin_%d_%d' % (i,i+1), Clostridium(energy='_%d' % i) + Inulin() >> Clostridium(energy= '_%d' % (i+1)), k_Clost_Inulin)
 for i in range(n_levels - 1)]

Parameter('k_Desulfo_ChondSulf', 1)

[Rule('DesulfoEatChond_%d_%d' % (i,i+1), Desulfobrivio(energy='_%d' % i) + ChondSulf() >> Desulfobrivio(energy= '_%d' % (i+1)), k_Desulfo_ChondSulf)
 for i in range(n_levels - 1)]

Parameter('k_Desulfo_Lactate', 1)

[Rule('DesulfoEatLactate_%d_%d' % (i,i+1), Desulfobrivio(energy='_%d' % i) + Lactate() >> Desulfobrivio(energy= '_%d' % (i+1)), k_Desulfo_Lactate)
 for i in range(n_levels - 1)]

Parameter('k_Bact_Glucose')

[Rule('BactEatGlucose_%d_%d' % (i,i+1), Bacteroides(energy='_%d' % i) + Glucose() >> Bacteroides(energy= '_%d' % (i+1)), k_Bact_Glucose)
 for i in range(n_levels - 1)]

Parameter('k_Bact_Lactose')

[Rule('BactEatLactose_%d_%d' % (i,i+1), Bacteroides(energy='_%d' % i) + Lactose() >> Bacteroides(energy= '_%d' % (i+1)), k_Bact_Lactose)
 for i in range(n_levels - 1)]

Parameter('k_Bact_Fructo')

[Rule('BactEatFructo_%d_%d' % (i,i+1), Bacteroides(energy='_%d' % i) + Fructo() >> Bacteroides(energy= '_%d' % (i+1)), k_Bact_Fructo)
 for i in range(n_levels - 1)]

Parameter('k_Clost_Glucose')

[Rule('ClostEatGlucose_%d_%d' % (i,i+1), Clostridium(energy='_%d' % i) + Glucose() >> Clostridium(energy= '_%d' % (i+1)), k_Clost_Glucose)
 for i in range(n_levels - 1)]

Parameter('k_Clost_Lactose')

[Rule('ClostEatLactose_%d_%d' % (i,i+1), Clostridium(energy='_%d' % i) + Lactose() >> Clostridium(energy= '_%d' % (i+1)), k_Clost_Lactose)
 for i in range(n_levels - 1)]

Parameter('k_Clost_Fructo')

[Rule('ClostEatFructo_%d_%d' % (i,i+1), Clostridium(energy='_%d' % i) + Fructo() >> Clostridium(energy= '_%d' % (i+1)), k_Clost_Fructo)
 for i in range(n_levels - 1)]

Parameter('k_Bifido_Glucose')

[Rule('BifidoEatGlucose_%d_%d' % (i,i+1), Bifidobacterium(energy='_%d' % i) + Glucose() >> Bifidobacterium(energy= '_%d' % (i+1)), k_Bifido_Glucose)
 for i in range(n_levels - 1)]

Parameter('k_Bifido_Lactose')

[Rule('BifidoEatLactose_%d_%d' % (i,i+1), Bifidobacterium(energy='_%d' % i) + Lactose() >> Bifidobacterium(energy= '_%d' % (i+1)), k_Bifido_Lactose)
 for i in range(n_levels - 1)]

Parameter('k_Bifido_Fructo')

[Rule('BifidoEatFructo_%d_%d' % (i,i+1), Bifidobacterium(energy='_%d' % i) + Fructo() >> Bifidobacterium(energy= '_%d' % (i+1)), k_Bifido_Fructo)
 for i in range(n_levels - 1)]

Parameter('k_Bifido_Lactate')

Rule('BifidoMakeLactate', Bifidobacterium() >> Bifidobacterium() + Lactate(), k_Bifido_Lactate)

#energy loss rules

Parameter('k_energy_loss', 14.4/26.3)

[Rule('Bact_Loss_%d_%d' % (i,i-1), Bacteroides(energy='_%d' % i) >> Bacteroides(energy= '_%d' % (i-1)), k_energy_loss)
 for i in range(1, n_levels)]

[Rule('Clost_Loss_%d_%d' % (i,i-1), Clostridium(energy='_%d' % i) >> Clostridium(energy= '_%d' % (i-1)), k_energy_loss)
 for i in range(1, n_levels)]

[Rule('Desulfo_Loss_%d_%d' % (i,i-1), Desulfobrivio(energy='_%d' % i) >> Desulfobrivio(energy= '_%d' % (i-1)), k_energy_loss)
 for i in range(1, n_levels)]

[Rule('Bifido_Loss_%d_%d' % (i,i-1), Bifidobacterium(energy='_%d' % i) >> Bifidobacterium(energy= '_%d' % (i-1)), k_energy_loss)
 for i in range(1, n_levels)]

#monomer division
#TODO: fix rules to conserve energy (some kind of if statament) + look for death + doubling times + consumption rates of metabolites
div_threshold = int(n_levels/2)

td_Bact = 1 #find that number
td_Clost = 1
td_Desulfo = 1
td_Bifido = 1
Parameter('k_Bact_division', np.log(2)/td_Bact)
Parameter('k_Clost_division', np.log(2)/td_Clost)
Parameter('k_Desulfo_division', np.log(2)/td_Desulfo)
Parameter('k_Bifido_division', np.log(2)/td_Bifido)

[Rule('Bact_divides_%d' % i, Bacteroides(energy='_%d' % i) >> Bacteroides(energy='_%d' % (i/2)) + Bacteroides(energy='_%d' % (i/2)), k_Bact_division)
 for i in range(div_threshold, n_levels)]

[Rule('Clost_divides_%d' % i, Clostridium(energy='_%d' % i) >> Clostridium(energy='_%d' % (i/2)) + Clostridium(energy='_%d' % (i/2)), k_Clost_division)
 for i in range(div_threshold, n_levels)]

[Rule('Desulfo_divides_%d' % i, Desulfobrivio(energy='_%d' % i) >> Desulfobrivio(energy='_%d' % (i/2)) + Desulfobrivio(energy='_%d' % (i/2)), k_Desulfo_division)
 for i in range(div_threshold, n_levels)]

[Rule('Bifido_divides_%d' % i, Bifidobacterium(energy='_%d' % i) >> Bifidobacterium(energy='_%d' % (i/2)) + Bifidobacterium(energy='_%d' % (i/2)), k_Bifido_division)
 for i in range(div_threshold, n_levels)]

print(model.rules)





