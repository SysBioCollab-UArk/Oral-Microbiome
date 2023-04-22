from pysb import *
from pysb.simulator import ScipyOdeSimulator
import numpy as np
import matplotlib.pyplot as plt
from ReadGutlogo import read_Gutlogo
from sympy import Piecewise

Model()
n_levels = 10
t_step = 26.3 #seconds

#todo: add carrying capacity rules*
#yes, stuck bacteria can reproduce, and offspring bacteria is unstuck
#parameters to find: hungry consumption rate, basal consumption rate, lactate production, stuck rates
#bacteria
Monomer('Bacteroides', ['energy', 'stuck'], {'energy': ['_%d' % i for i in range(n_levels)], 'stuck': ['u', 's', 'p']})
Monomer('Clostridium', ['energy', 'stuck'], {'energy': ['_%d' % i for i in range(n_levels)], 'stuck': ['u', 's', 'p']})
Monomer('Bifidobacterium', ['energy', 'stuck'], {'energy': ['_%d' % i for i in range(n_levels)], 'stuck': ['u', 's', 'p']})
Monomer('Desulfobrivio', ['energy', 'stuck'], {'energy': ['_%d' % i for i in range(n_levels)], 'stuck': ['u', 's', 'p']})

Parameter('Bact_0', 5490)
Parameter('Clost_0', 921)
Parameter('Bifido_0', 23562)
Parameter('Desulfo_0', 70)

Initial(Bacteroides(energy = '_%d' % (n_levels - 1), stuck ='u'), Bact_0)
Initial(Clostridium(energy = '_%d' % (n_levels - 1), stuck ='u'), Clost_0)
Initial(Bifidobacterium(energy = '_%d' % (n_levels - 1), stuck ='u'), Bifido_0)
Initial(Desulfobrivio(energy = '_%d' % (n_levels - 1), stuck ='u'), Desulfo_0)


#metabolites
Monomer('Inulin')
Monomer('Glucose')
Monomer('Lactose')
Monomer('Fructo')
Monomer('ChondSulf')
Monomer('Lactate')

Parameter('Inulin_0', 100)
Parameter('Glucose_0', 100)
Parameter('Lactose_0', 100)
Parameter('Fructo_0', 100)
Parameter('ChondSulf_0', 100)
Parameter('Lactate_0', 100)

Initial(Inulin(), Inulin_0)
Initial(Glucose(), Glucose_0)
Initial(Lactose(), Lactose_0)
Initial(Fructo(), Fructo_0)
Initial(ChondSulf(), ChondSulf_0)
Initial(Lactate(), Lactate_0)


#initial metabolites

hung_threshold = int(0.8*n_levels - 1)

Observable('Bact_tot', Bacteroides())
Observable('Clost_tot', Clostridium())
Observable('Bifido_tot', Bifidobacterium())
Observable('Desulfo_tot', Desulfobrivio())
Observable('Pop_tot', Bacteroides() + Clostridium() + Bifidobacterium() + Desulfobrivio())
Observable('Metab_tot', Inulin() + Glucose() + Lactose() + Fructo() + ChondSulf() + Lactate())
Observable('Bact_E100', Bacteroides(energy = '_%d' % (n_levels - 1)))
Observable('Bact_E90', Bacteroides(energy = '_%d' % (n_levels - 2)))
Observable('Bact_E80', Bacteroides(energy = '_%d' % (n_levels - 3)))
Observable('Bact_E70', Bacteroides(energy = '_%d' % (n_levels - 4)))

obs_to_plot = ['Clost_tot', 'Bifido_tot', 'Desulfo_tot', 'Bact_tot']

Hungry_bact = [Bacteroides(energy = '_%d' % i) for i in range(hung_threshold + 1)] + \
              [Clostridium(energy = '_%d' % i) for i in range(hung_threshold + 1)] + \
              [Bifidobacterium(energy = '_%d' % i) for i in range(hung_threshold + 1)] + \
              [Desulfobrivio(energy = '_%d' % i) for i in range(hung_threshold + 1)]
Hungry_obs = Hungry_bact[0]
for i in range(1, len(Hungry_bact)):
 Hungry_obs += Hungry_bact[i]
Observable('Hungry_tot', Hungry_obs)

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
#  k_hungry = 1/(N*M), N = total number of hungry bacteria, M = total number of metabolites
# (Nb/N * 1/Nb) * (Mi/M * 1/Mi) = prob of selecting any given bacteroide with any given inulin
#todo finish 25s and 50s, run 1000 ticks and download csv for netlogo
n_Bact_Inulin = 250  #energy increase in gutlogo code
Parameter('k_Bact_Inulin_Basal', 0)
#Parameter('k_Bact_Inulin_Hungry', 10)
Expression('k_Bact_Inulin_Hungry', Piecewise((0, (Metab_tot < 1) | (Hungry_tot < 1)),
                                             (n_Bact_Inulin/(Metab_tot*Hungry_tot)/t_step, True)))

[Rule('BactEatInulin_Basal_%d_%d' % (i,i+1), Bacteroides(energy='_%d' % i) + Inulin() >> Bacteroides(energy= '_%d' % (i+1)), k_Bact_Inulin_Basal)
 for i in range(hung_threshold, n_levels -1)]
[Rule('BactEatInulin_Hungry_%d_%d' % (i,i+1), Bacteroides(energy='_%d' % i) + Inulin() >> Bacteroides(energy= '_%d' % (i+1)), k_Bact_Inulin_Hungry)
 for i in range(0, hung_threshold)]

n_Clost_Inulin = 250
Parameter('k_Clost_Inulin_Basal', 0)
#Parameter('k_Clost_Inulin_Hungry', 10)
Expression('k_Clost_Inulin_Hungry', Piecewise((0, (Metab_tot < 1) | (Hungry_tot < 1)),
                                             (n_Clost_Inulin/(Metab_tot*Hungry_tot)/t_step, True)))
[Rule('ClostEatInulin_Basal_%d_%d' % (i,i+1), Clostridium(energy='_%d' % i) + Inulin() >> Clostridium(energy= '_%d' % (i+1)), k_Clost_Inulin_Basal)
 for i in range(hung_threshold, n_levels - 1)]
[Rule('ClostEatInulin_Hungry_%d_%d' % (i,i+1), Clostridium(energy='_%d' % i) + Inulin() >> Clostridium(energy= '_%d' % (i+1)), k_Clost_Inulin_Hungry)
 for i in range(0, hung_threshold)]

#todo we multiplied by 10, why??
#todo increase energy by chunks, decrease rate constant, change n_levels to 100, add species for new energies above n_levels
n_Desulfo_Choldsulf = 250
Parameter('k_Desulfo_ChondSulf_Basal', 0)
#Parameter('k_Desulfo_ChondSulf_Hungry', 10)
Expression('k_Desulfo_ChondSulf_Hungry', Piecewise((0, (Metab_tot < 1) | (Hungry_tot < 1)),
                                             (n_Desulfo_Choldsulf/(Metab_tot*Hungry_tot)/t_step, True)))
[Rule('DesulfoEatChond_Basal_%d_%d' % (i,i+1), Desulfobrivio(energy='_%d' % i) + ChondSulf() >> Desulfobrivio(energy= '_%d' % (i+1)), k_Desulfo_ChondSulf_Basal)
 for i in range(hung_threshold, n_levels - 1)]
[Rule('DesulfoEatChond_Hungry_%d_%d' % (i,i+1), Desulfobrivio(energy='_%d' % i) + ChondSulf() >> Desulfobrivio(energy= '_%d' % (i+1)), k_Desulfo_ChondSulf_Hungry)
 for i in range(0, hung_threshold)]

n_Desulfo_Lactate = 500
Parameter('k_Desulfo_Lactate_Basal', 0)
#Parameter('k_Desulfo_Lactate_Hungry', 10)
Expression('k_Desulfo_Lactate_Hungry', Piecewise((0, (Metab_tot < 1) | (Hungry_tot < 1)),
                                             (n_Desulfo_Lactate/(Metab_tot*Hungry_tot)/t_step, True)))
[Rule('DesulfoEatLactate_Basal_%d_%d' % (i,i+1), Desulfobrivio(energy='_%d' % i) + Lactate() >> Desulfobrivio(energy= '_%d' % (i+1)), k_Desulfo_Lactate_Basal)
 for i in range(hung_threshold, n_levels - 1)]
[Rule('DesulfoEatLactate_Hungry_%d_%d' % (i,i+1), Desulfobrivio(energy='_%d' % i) + Lactate() >> Desulfobrivio(energy= '_%d' % (i+1)), k_Desulfo_Lactate_Hungry)
 for i in range(0, hung_threshold)]

n_Bact_Glucose = 500
Parameter('k_Bact_Glucose_Basal', 0)
#Parameter('k_Bact_Glucose_Hungry', 10)
Expression('k_Bact_Glucose_Hungry', Piecewise((0, (Metab_tot < 1) | (Hungry_tot < 1)),
                                             (n_Bact_Glucose/(Metab_tot*Hungry_tot)/t_step, True)))
[Rule('BactEatGlucose_Basal_%d_%d' % (i,i+1), Bacteroides(energy='_%d' % i) + Glucose() >> Bacteroides(energy= '_%d' % (i+1)), k_Bact_Glucose_Basal)
 for i in range(hung_threshold, n_levels - 1)]
[Rule('BactEatGlucose_Hungry_%d_%d' % (i,i+1), Bacteroides(energy='_%d' % i) + Glucose() >> Bacteroides(energy= '_%d' % (i+1)), k_Bact_Glucose_Hungry)
 for i in range(0, hung_threshold)]

n_Bact_Lactose = 250
Parameter('k_Bact_Lactose_Basal', 0)
#Parameter('k_Bact_Lactose_Hungry', 10)
Expression('k_Bact_Lactose_Hungry', Piecewise((0, (Metab_tot < 1) | (Hungry_tot < 1)),
                                             (n_Bact_Lactose/(Metab_tot*Hungry_tot)/t_step, True)))
[Rule('BactEatLactose_Basal_%d_%d' % (i,i+1), Bacteroides(energy='_%d' % i) + Lactose() >> Bacteroides(energy= '_%d' % (i+1)), k_Bact_Lactose_Basal)
 for i in range(hung_threshold, n_levels - 1)]
[Rule('BactEatLactose_Hungry_%d_%d' % (i,i+1), Bacteroides(energy='_%d' % i) + Lactose() >> Bacteroides(energy= '_%d' % (i+1)), k_Bact_Lactose_Hungry)
 for i in range(0, hung_threshold)]

n_Bact_Fructo = 250
Parameter('k_Bact_Fructo_Basal', 0)
#Parameter('k_Bact_Fructo_Hungry', 10)
Expression('k_Bact_Fructo_Hungry', Piecewise((0, (Metab_tot < 1) | (Hungry_tot < 1)),
                                             (n_Bact_Fructo/(Metab_tot*Hungry_tot)/t_step, True)))
[Rule('BactEatFructo_Basal_%d_%d' % (i,i+1), Bacteroides(energy='_%d' % i) + Fructo() >> Bacteroides(energy= '_%d' % (i+1)), k_Bact_Fructo_Basal)
 for i in range(hung_threshold, n_levels - 1)]
[Rule('BactEatFructo_Hungry_%d_%d' % (i,i+1), Bacteroides(energy='_%d' % i) + Fructo() >> Bacteroides(energy= '_%d' % (i+1)), k_Bact_Fructo_Hungry)
 for i in range(0, hung_threshold)]

n_Clost_Glucose = 500
Parameter('k_Clost_Glucose_Basal', 0)
#Parameter('k_Clost_Glucose_Hungry', 10)
Expression('k_Clost_Glucose_Hungry', Piecewise((0, (Metab_tot < 1) | (Hungry_tot < 1)),
                                             (n_Clost_Glucose/(Metab_tot*Hungry_tot)/t_step, True)))
[Rule('ClostEatGlucose_Basal_%d_%d' % (i,i+1), Clostridium(energy='_%d' % i) + Glucose() >> Clostridium(energy= '_%d' % (i+1)), k_Clost_Glucose_Basal)
 for i in range(hung_threshold, n_levels - 1)]
[Rule('ClostEatGlucose_Hungry_%d_%d' % (i,i+1), Clostridium(energy='_%d' % i) + Glucose() >> Clostridium(energy= '_%d' % (i+1)), k_Clost_Glucose_Hungry)
 for i in range(0, hung_threshold)]

n_Clost_Lactose = 250
Parameter('k_Clost_Lactose_Basal', 0)
#Parameter('k_Clost_Lactose_Hungry', 10)
Expression('k_Clost_Lactose_Hungry', Piecewise((0, (Metab_tot < 1) | (Hungry_tot < 1)),
                                             (n_Clost_Lactose/(Metab_tot*Hungry_tot)/t_step, True)))
[Rule('ClostEatLactose_Basal_%d_%d' % (i,i+1), Clostridium(energy='_%d' % i) + Lactose() >> Clostridium(energy= '_%d' % (i+1)), k_Clost_Lactose_Basal)
 for i in range(hung_threshold, n_levels - 1)]
[Rule('ClostEatLactose_Hungry_%d_%d' % (i,i+1), Clostridium(energy='_%d' % i) + Lactose() >> Clostridium(energy= '_%d' % (i+1)), k_Clost_Lactose_Hungry)
 for i in range(0, hung_threshold)]

n_Clost_Fructo = 250
Parameter('k_Clost_Fructo_Basal', 0)
#Parameter('k_Clost_Fructo_Hungry', 10)
Expression('k_Clost_Fructo_Hungry', Piecewise((0, (Metab_tot < 1) | (Hungry_tot < 1)),
                                             (n_Clost_Fructo/(Metab_tot*Hungry_tot)/t_step, True)))
[Rule('ClostEatFructo_Basal_%d_%d' % (i,i+1), Clostridium(energy='_%d' % i) + Fructo() >> Clostridium(energy= '_%d' % (i+1)), k_Clost_Fructo_Basal)
 for i in range(hung_threshold, n_levels - 1)]
[Rule('ClostEatFructo_Hungry_%d_%d' % (i,i+1), Clostridium(energy='_%d' % i) + Fructo() >> Clostridium(energy= '_%d' % (i+1)), k_Clost_Fructo_Hungry)
 for i in range(0, hung_threshold)]

n_Bifido_Glucose = 250
Parameter('k_Bifido_Glucose_Basal', 0)
# Parameter('k_Bifido_Glucose_Hungry', 10)
Expression('k_Bifido_Glucose_Hungry', Piecewise((0, (Metab_tot < 1) | (Hungry_tot < 1)),
                                             (n_Bifido_Glucose/(Metab_tot*Hungry_tot)/t_step, True)))
[Rule('BifidoEatGlucose_Basal_%d_%d' % (i,i+1), Bifidobacterium(energy='_%d' % i) + Glucose() >> Bifidobacterium(energy= '_%d' % (i+1)), k_Bifido_Glucose_Basal)
 for i in range(hung_threshold, n_levels - 1)]
[Rule('BifidoEatGlucose_Hungry_%d_%d' % (i,i+1), Bifidobacterium(energy='_%d' % i) + Glucose() >> Bifidobacterium(energy= '_%d' % (i+1)), k_Bifido_Glucose_Hungry)
 for i in range(0, hung_threshold)]

n_Bifido_Lactose = 500
Parameter('k_Bifido_Lactose_Basal', 0)
# Parameter('k_Bifido_Lactose_Hungry', 10)
Expression('k_Bifido_Lactose_Hungry', Piecewise((0, (Metab_tot < 1) | (Hungry_tot < 1)),
                                             (n_Bifido_Lactose/(Metab_tot*Hungry_tot)/t_step, True)))
[Rule('BifidoEatLactose_Basal_%d_%d' % (i,i+1), Bifidobacterium(energy='_%d' % i) + Lactose() >> Bifidobacterium(energy= '_%d' % (i+1)), k_Bifido_Lactose_Basal)
 for i in range(hung_threshold, n_levels - 1)]
[Rule('BifidoEatLactose_Hungry_%d_%d' % (i,i+1), Bifidobacterium(energy='_%d' % i) + Lactose() >> Bifidobacterium(energy= '_%d' % (i+1)), k_Bifido_Lactose_Hungry)
 for i in range(0, hung_threshold)]

n_Bifido_Fructo = 500
Parameter('k_Bifido_Fructo_Basal', 0)
# Parameter('k_Bifido_Fructo_Hungry', 10)
Expression('k_Bifido_Fructo_Hungry', Piecewise((0, (Metab_tot < 1) | (Hungry_tot < 1)),
                                             (n_Bifido_Fructo/(Metab_tot*Hungry_tot)/t_step, True)))
[Rule('BifidoEatFructo_Basal_%d_%d' % (i,i+1), Bifidobacterium(energy='_%d' % i) + Fructo() >> Bifidobacterium(energy= '_%d' % (i+1)), k_Bifido_Fructo_Basal)
 for i in range(hung_threshold, n_levels - 1)]
[Rule('BifidoEatFructo_Hungry_%d_%d' % (i,i+1), Bifidobacterium(energy='_%d' % i) + Fructo() >> Bifidobacterium(energy= '_%d' % (i+1)), k_Bifido_Fructo_Hungry)
 for i in range(0, hung_threshold)]

#lactate production rule
Parameter('k_Bifido_Lactate', 0.005/t_step)

Rule('BifidoMakeLactate', Bifidobacterium() >> Bifidobacterium() + Lactate(), k_Bifido_Lactate)

#energy loss rules
Parameter('k_energy_loss', 1/(1440/n_levels * t_step))

[Rule('Bact_Loss_%d_%d' % (i,i-1), Bacteroides(energy='_%d' % i) >> Bacteroides(energy= '_%d' % (i-1)), k_energy_loss)
 for i in range(1, n_levels)]

[Rule('Clost_Loss_%d_%d' % (i,i-1), Clostridium(energy='_%d' % i) >> Clostridium(energy= '_%d' % (i-1)), k_energy_loss)
 for i in range(1, n_levels)]

[Rule('Desulfo_Loss_%d_%d' % (i,i-1), Desulfobrivio(energy='_%d' % i) >> Desulfobrivio(energy= '_%d' % (i-1)), k_energy_loss)
 for i in range(1, n_levels)]

[Rule('Bifido_Loss_%d_%d' % (i,i-1), Bifidobacterium(energy='_%d' % i) >> Bifidobacterium(energy= '_%d' % (i-1)), k_energy_loss)
 for i in range(1, n_levels)]

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

[Rule('Bact_divides_%d' % i, Bacteroides(energy='_%d' % i) >> Bacteroides(energy='_%d' % (i/2)) + Bacteroides(energy='_%d' % ((i + 1)/2), stuck = 'u'), k_Bact_division)
 for i in range(div_threshold, n_levels)]

[Rule('Clost_divides_%d' % i, Clostridium(energy='_%d' % i) >> Clostridium(energy='_%d' % (i/2)) + Clostridium(energy='_%d' % ((i + 1)/2), stuck = 'u'), k_Clost_division)
 for i in range(div_threshold, n_levels)]

[Rule('Desulfo_divides_%d' % i, Desulfobrivio(energy='_%d' % i) >> Desulfobrivio(energy='_%d' % (i/2)) + Desulfobrivio(energy='_%d' % ((i + 1)/2), stuck = 'u'), k_Desulfo_division)
 for i in range(div_threshold, n_levels)]

[Rule('Bifido_divides_%d' % i, Bifidobacterium(energy='_%d' % i) >> Bifidobacterium(energy='_%d' % (i/2)) + Bifidobacterium(energy='_%d' % ((i + 1)/2), stuck = 'u'), k_Bifido_division)
 for i in range(div_threshold, n_levels)]

#death rules
Parameter('k_Death', 1/t_step)
Rule('Bact_death', Bacteroides(energy = '_0') >> None, k_Death)
Rule('Clost_death', Clostridium(energy = '_0') >> None, k_Death)
Rule('Desulfo_death', Desulfobrivio(energy = '_0') >> None, k_Death)
Rule('Bifido_death', Bifidobacterium(energy = '_0') >> None, k_Death)
#leave rate
#rate of sticking, unsticking, and time to leave system (k) - find this

#Parameter('k_Bact_stuck', 1/(t_step*100))
#Parameter('k_Clost_stuck', 1/(t_step*100))
#Parameter('k_Desulfo_stuck', 1/(t_step*100))
#Parameter('k_Bifido_stuck', 1/(t_step*100))
#Parameter('k_Bact_unstuck', 1/(t_step*10))
#Parameter('k_Clost_unstuck', 1/(t_step*10))
#Parameter('k_Desulfo_unstuck', 1/(t_step*10))
#Parameter('k_Bifido_unstuck', 1/(t_step*10))
#Parameter('k_Bact_permstuck', 1/(t_step*20))
#Parameter('k_Clost_permstuck', 1/(t_step*20))
#Parameter('k_Desulfo_permstuck', 1/(t_step*20))
#Parameter('k_Bifido_permstuck', 1/(t_step*20))
#Parameter('k_Bact_removed',1/(33*2.25*10000*t_step) - k_Bact_stuck + k_Bact_stuck*k_Bact_unstuck)
#Parameter('k_Clost_removed', 1/(33*2.25*10000*t_step) - k_Clost_stuck + k_Clost_stuck*k_Clost_unstuck)
#Parameter('k_Desulfo_removed', 1/(33*2.25*10000*t_step) - k_Desulfo_stuck + k_Desulfo_stuck*k_Delsulfo_unstuck)
#Parameter('k_Bifido_removed', 1/(33*2.25*10000*t_step) - k_Bifido_stuck + k_Bifido_stuck*k_Bifido_unstuck)

#removal rules

flow_rate = 2.22 #cm/min
Parameter('k_Bact_unstuck_removed', flow_rate/(7.98 * 100)/60) #/s (7.98 cm is the length of an element; there are 100 total elements; 60 seconds per minute)
Parameter('k_Clost_unstuck_removed', flow_rate/(7.98 * 100)/60) #/s (7.98 cm is the length of an element; there are 100 total elements; 60 seconds per minute)
Parameter('k_Desulfo_unstuck_removed', flow_rate/(7.98 * 100)/60) #/s (7.98 cm is the length of an element; there are 100 total elements; 60 seconds per minute)
Parameter('k_Bifido_unstuck_removed', flow_rate/(7.98 * 100)/60) #/s (7.98 cm is the length of an element; there are 100 total elements; 60 seconds per minute)

Rule('Bact_unstuck_removal', Bacteroides(stuck = 'u') >> None, k_Bact_unstuck_removed)
Rule('Clost_unstuck_removal', Clostridium(stuck = 'u') >> None, k_Clost_unstuck_removed)
Rule('Desulfo_unstuck_removal', Desulfobrivio(stuck = 'u') >> None, k_Desulfo_unstuck_removed)
Rule('Bifido_unstuck_removal', Bifidobacterium(stuck = 'u') >> None, k_Bifido_unstuck_removed)


#stuck rules
maxStuckChance = 0.5 #probability
midStuckConc = 10 #concentration
lowStuckBound = 0.02 #probability

Expression('k_Bact_stuck',Piecewise(
 (0, maxStuckChance * midStuckConc / (midStuckConc + Pop_tot) / t_step < lowStuckBound),
 (maxStuckChance * midStuckConc / (midStuckConc + Pop_tot) / t_step, True)))
Parameter('k_Bact_unstuck', 0.1/t_step)
Parameter('k_Bact_permstuck', 0.05/t_step)

Rule('Bact_stuck', Bacteroides(stuck='u') | Bacteroides(stuck= 's'), k_Bact_stuck, k_Bact_unstuck)
Rule('Bact_permstuck', Bacteroides(stuck='s') >> Bacteroides(stuck= 'p'), k_Bact_permstuck)

Expression('k_Clost_stuck',Piecewise(
 (0, maxStuckChance * midStuckConc / (midStuckConc + Pop_tot) / t_step < lowStuckBound),
 (maxStuckChance * midStuckConc / (midStuckConc + Pop_tot) / t_step, True)))
Parameter('k_Clost_unstuck', 0.1/t_step)
Parameter('k_Clost_permstuck', 0.05/t_step)

Rule('Clost_stuck', Clostridium(stuck='u') | Clostridium(stuck= 's'), k_Clost_stuck, k_Clost_unstuck)
Rule('Clost_permstuck', Clostridium(stuck='s') >> Clostridium(stuck= 'p'), k_Clost_permstuck)

Expression('k_Desulfo_stuck',Piecewise(
 (0, maxStuckChance * midStuckConc / (midStuckConc + Pop_tot) / t_step < lowStuckBound),
 (maxStuckChance * midStuckConc / (midStuckConc + Pop_tot) / t_step, True)))
Parameter('k_Desulfo_unstuck', 0.1/t_step)
Parameter('k_Desulfo_permstuck', 0.05/t_step)

Rule('Desulfo_stuck', Desulfobrivio(stuck='u') | Desulfobrivio(stuck= 's'), k_Desulfo_stuck, k_Desulfo_unstuck)
Rule('Desulfo_permstuck', Desulfobrivio(stuck='s') >> Desulfobrivio(stuck= 'p'), k_Desulfo_permstuck)

Expression('k_Bifido_stuck',Piecewise(
 (0, maxStuckChance * midStuckConc / (midStuckConc + Pop_tot) / t_step < lowStuckBound),
 (maxStuckChance * midStuckConc / (midStuckConc + Pop_tot) / t_step, True)))
Parameter('k_Bifido_unstuck', 0.1/t_step)
Parameter('k_Bifido_permstuck', 0.05/t_step)

Rule('Bifido_stuck', Bifidobacterium(stuck='u') | Bifidobacterium(stuck= 's'), k_Bifido_stuck, k_Bifido_unstuck)
Rule('Bifido_permstuck', Bifidobacterium(stuck='s') >> Bifidobacterium(stuck= 'p'), k_Bifido_permstuck)

#todo inflow rules
#bacterial creation rules
Parameter('k_Bact_creation', 0)
Rule('Bact_creation', None >> Bacteroides(energy = '_%d' % (n_levels - 1), stuck = 'u'), k_Bact_creation)

Parameter('k_Clost_creation', 0)
Rule('Clost_creation', None >> Clostridium(energy = '_%d' % (n_levels - 1), stuck = 'u'), k_Clost_creation)

Parameter('k_Desulfo_creation', 0)
Rule('Desulfo_creation', None >> Desulfobrivio(energy = '_%d' % (n_levels - 1), stuck = 'u'), k_Desulfo_creation)

Parameter('k_Bifido_creation', 0)
Rule('Bifido_creation', None >> Bifidobacterium(energy = '_%d' % (n_levels - 1), stuck = 'u'), k_Bifido_creation)

print(model.rules)

#simulations
n_steps = 5000
t_span = np.linspace(0, n_steps*t_step, int(n_steps)*1+1)
sim = ScipyOdeSimulator(model, t_span, verbose = True)
result = sim.run()

for obs in obs_to_plot:
 plt.plot(range(n_steps + 1), result.observables[obs], label = obs)
 plt.xlabel('time (number of time steps)')
 plt.ylabel('# of cells')
 #plt.yscale('log', base = 10)
 plt.legend(loc=0)

print(result.observables['Bact_tot'])
print(n_steps*t_step)
print(t_span)
print(result.expressions['k_Bact_Inulin_Hungry'])
print(result.observables['Hungry_tot'])
print(result.observables['Metab_tot'])
print(result.observables['Pop_tot'])

plt.tight_layout()

plt.figure()
x, pops, label = read_Gutlogo('populations-4.csv')
for y, l in zip(pops, label):
    plt.plot(x, y, label=l)

plt.xlabel('Time (number of time steps)')
plt.ylabel('Population (number of agents)')
plt.title('Kinetics Model of Gut Microbiome')
plt.legend(loc=0)

plt.figure()
for obs in ['Bact_E100', 'Bact_E90', 'Bact_E80', 'Bact_E70']:
    plt.plot(range(n_steps + 1), result.observables[obs], label = obs)
plt.xlabel('time (number of time steps)')
plt.ylabel('# of cells')
plt.legend(loc=0)

plt.show()

#todo: upload csv to project (oral micro folder top left), create a new file to read the data in and plot it (np.genfromtxt), google read data from a file in python
#todo: run some longer sims, push to github
#"Time","Closts","Bifidos","Desulfos","Bacteroides"
