from pysb import *
from pysb.simulator import ScipyOdeSimulator
import numpy as np
import matplotlib.pyplot as plt

Model()
t_step = 26.3 #seconds
n_levels = 10

#bacteria
Monomer('Bacteroides', ['energy', 'stuck'], {'energy': ['_%d' % i for i in range(n_levels)], 'stuck': ['u', 's', 'p']})
Monomer('Clostridium', ['energy', 'stuck'], {'energy': ['_%d' % i for i in range(n_levels)], 'stuck': ['u', 's', 'p']})
Monomer('Bifidobacterium', ['energy', 'stuck'], {'energy': ['_%d' % i for i in range(n_levels)], 'stuck': ['u', 's', 'p']})
Monomer('Desulfobrivio', ['energy', 'stuck'], {'energy': ['_%d' % i for i in range(n_levels)], 'stuck': ['u', 's', 'p']})

#metabolites
Monomer('Inulin')
Monomer('Glucose')
Monomer('Lactose')
Monomer('Fructo')
Monomer('ChondSulf')
Monomer('Lactate')

Observable('Bact_tot', Bacteroides())
Observable('Clost_tot', Clostridium())
Observable('Bifido_tot', Bifidobacterium())
Observable('Desulfo_tot', Desulfobrivio())

#metabolite production rules

Parameter('k_Inulin_prod', 10/t_step)
Parameter('k_Glucose_prod', 30/t_step)
Parameter('k_Lactose_prod', 15/t_step)
Parameter('k_Fructo_prod', 25/t_step)
Parameter('k_ChondSulf_prod', 0.1/t_step)
Parameter('k_Lactate_prod', 0/t_step)

Rule('Inulin_prod', None >> Inulin(), k_Inulin_prod)
Rule('Glucose_prod', None >> Glucose(), k_Glucose_prod)
Rule('Lactose_prod', None >> Lactose(), k_Lactose_prod)
Rule('Fructo_prod', None >> Fructo(), k_Fructo_prod)
Rule('ChondSulf_prod', None >> ChondSulf(), k_ChondSulf_prod)
Rule('Lactate_prod', None >> Lactate(), k_Lactate_prod)

#bacterial consumption rules
hung_threshold = int(0.8*n_levels - 1)
Parameter('k_Bact_Inulin_Basal', 1)
Parameter('k_Bact_Inulin_Hungry', 10)

[Rule('BactEatInulin_Basal_%d_%d' % (i,i+1), Bacteroides(energy='_%d' % i) + Inulin() >> Bacteroides(energy= '_%d' % (i+1)), k_Bact_Inulin_Basal)
 for i in range(hung_threshold, n_levels -1)]

[Rule('BactEatInulin_Hungry_%d_%d' % (i,i+1), Bacteroides(energy='_%d' % i) + Inulin() >> Bacteroides(energy= '_%d' % (i+1)), k_Bact_Inulin_Hungry)
 for i in range(0, hung_threshold)]

#Todo: finish adding basal and hungry rules below.
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

Parameter('k_energy_loss', 14.4/t_step)

[Rule('Bact_Loss_%d_%d' % (i,i-1), Bacteroides(energy='_%d' % i) >> Bacteroides(energy= '_%d' % (i-1)), k_energy_loss)
 for i in range(1, n_levels)]

[Rule('Clost_Loss_%d_%d' % (i,i-1), Clostridium(energy='_%d' % i) >> Clostridium(energy= '_%d' % (i-1)), k_energy_loss)
 for i in range(1, n_levels)]

[Rule('Desulfo_Loss_%d_%d' % (i,i-1), Desulfobrivio(energy='_%d' % i) >> Desulfobrivio(energy= '_%d' % (i-1)), k_energy_loss)
 for i in range(1, n_levels)]

[Rule('Bifido_Loss_%d_%d' % (i,i-1), Bifidobacterium(energy='_%d' % i) >> Bifidobacterium(energy= '_%d' % (i-1)), k_energy_loss)
 for i in range(1, n_levels)]

#TODO: fix rules to conserve energy (some kind of if statament) + look for death + doubling times + consumption rates of metabolites
#death: no mention of death, but any agent that moves past the final horizantal element is removed from the simulation
#death: bacteria move from element to element as a function of velocity of the flow field (displacement per timestep)
#death: probabilistic function controls likelyhood of getting stuck in the membrane. Probability of escaping is 10% (per some time step?) and 5% to become perm embedded
#since our program does not have a spacial element we will have to add a death element instead
#DT: all genera have same doubling times and they are all constant: 330 timesteps I believe
#MC: probability of consuming metabolite is function of metabolite concentrate and number of other hungry (energy<8) bacteria on the element
#MC: the function is not shown but I believe a good starting point would be (# of metabolites /# of hungry bacteria) per time step

div_threshold = int(n_levels/2)
td_Bact = 330*t_step
td_Clost = 330*t_step
td_Desulfo = 330*t_step
td_Bifido = 330*t_step
Parameter('k_Bact_division', np.log(2)/td_Bact)
Parameter('k_Clost_division', np.log(2)/td_Clost)
Parameter('k_Desulfo_division', np.log(2)/td_Desulfo)
Parameter('k_Bifido_division', np.log(2)/td_Bifido)

#bacteria division rules

[Rule('Bact_divides_%d' % i, Bacteroides(energy='_%d' % i) >> Bacteroides(energy='_%d' % (i/2)) + Bacteroides(energy='_%d' % ((i + 1)/2)), k_Bact_division)
 for i in range(div_threshold, n_levels)]

[Rule('Clost_divides_%d' % i, Clostridium(energy='_%d' % i) >> Clostridium(energy='_%d' % (i/2)) + Clostridium(energy='_%d' % ((i + 1)/2)), k_Clost_division)
 for i in range(div_threshold, n_levels)]

[Rule('Desulfo_divides_%d' % i, Desulfobrivio(energy='_%d' % i) >> Desulfobrivio(energy='_%d' % (i/2)) + Desulfobrivio(energy='_%d' % ((i + 1)/2)), k_Desulfo_division)
 for i in range(div_threshold, n_levels)]

[Rule('Bifido_divides_%d' % i, Bifidobacterium(energy='_%d' % i) >> Bifidobacterium(energy='_%d' % (i/2)) + Bifidobacterium(energy='_%d' % ((i + 1)/2)), k_Bifido_division)
 for i in range(div_threshold, n_levels)]

#death rules
Parameter('k_Death', 1/t_step)
Rule('Bact_death', Bacteroides(energy = '_0') >> None, k_Death)
Rule('Clost_death', Clostridium(energy = '_0') >> None, k_Death)
Rule('Desulfo_death', Desulfobrivio(energy = '_0') >> None, k_Death)
Rule('Bifido_death', Bifidobacterium(energy = '_0') >> None, k_Death)
#leave rate
#rate of sticking, unsticking, and time to leave system (k) - find this

Parameter('k_Bact_stuck', 1/(t_step*100))
Parameter('k_Clost_stuck', 1/(t_step*100))
Parameter('k_Desulfo_stuck', 1/(t_step*100))
Parameter('k_Bifido_stuck', 1/(t_step*100))
Parameter('k_Bact_unstuck', 1/(t_step*10))
Parameter('k_Clost_unstuck', 1/(t_step*10))
Parameter('k_Desulfo_unstuck', 1/(t_step*10))
Parameter('k_Bifido_unstuck', 1/(t_step*10))
Parameter('k_Bact_permstuck', 1/(t_step*20))
Parameter('k_Clost_permstuck', 1/(t_step*20))
Parameter('k_Desulfo_permstuck', 1/(t_step*20))
Parameter('k_Bifido_permstuck', 1/(t_step*20))
#Parameter('k_Bact_removed',1/(33*2.25*10000*t_step) - k_Bact_stuck + k_Bact_stuck*k_Bact_unstuck)
#Parameter('k_Clost_removed', 1/(33*2.25*10000*t_step) - k_Clost_stuck + k_Clost_stuck*k_Clost_unstuck)
#Parameter('k_Desulfo_removed', 1/(33*2.25*10000*t_step) - k_Desulfo_stuck + k_Desulfo_stuck*k_Delsulfo_unstuck)
#Parameter('k_Bifido_removed', 1/(33*2.25*10000*t_step) - k_Bifido_stuck + k_Bifido_stuck*k_Bifido_unstuck)

#todo: s,p, u rules
#removal rules
#Rule('Bact_removal', Bacteroides() >> None, k_Bact_removed)
#Rule('Clost_removal', Clostridium() >> None, k_Clost_removed)
#Rule('Desulfo_removal', Desulfobrivio() >> None, k_Desulfo_removed)
#Rule('Bifido_removal', Bifidobacterium() >> None, k_Bifido_removed)

print(model.rules)




