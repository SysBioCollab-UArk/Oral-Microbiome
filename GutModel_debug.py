from pysb import *
from pysb.simulator import ScipyOdeSimulator
import numpy as np
import matplotlib.pyplot as plt
from ReadGutlogo import read_Gutlogo
from sympy import Piecewise

# test comment
Model()
n_levels = 10
t_step = 60 # 26.3 seconds
hung_threshold = int(0.8*n_levels - 1)
deltaE_50 = int(np.ceil(n_levels/2))
deltaE_25 = int(np.ceil(deltaE_50/2))
print(deltaE_50)
print(deltaE_25)

# todo: add carrying capacity rules*
# yes, stuck bacteria can reproduce, and offspring bacteria is unstuck
# parameters to find: hungry consumption rate, basal consumption rate, lactate production, stuck rates

# bacteria
Monomer('Bacteroides', ['energy', 'stuck'],
        {'energy': ['_%d' % i for i in range(0, n_levels + deltaE_50)], 'stuck': ['u', 's', 'p']})

print(model.monomers)

Parameter('Bact_0', 5490)

Initial(Bacteroides(energy='_%d' % (n_levels - 1), stuck='u'), Bact_0)

# metabolites
Monomer('Inulin')
Monomer('Glucose')
Monomer('Lactose')
Monomer('Fructo')
Monomer('ChondSulf')
Monomer('Lactate')

Parameter('Inulin_0', 0)
Parameter('Glucose_0', 0)
Parameter('Lactose_0', 0)
Parameter('Fructo_0', 0)
Parameter('ChondSulf_0', 0)
Parameter('Lactate_0', 0)

Initial(Inulin(), Inulin_0)
Initial(Glucose(), Glucose_0)
Initial(Lactose(), Lactose_0)
Initial(Fructo(), Fructo_0)
Initial(ChondSulf(), ChondSulf_0)
Initial(Lactate(), Lactate_0)

# observables

Observable('Bact_tot', Bacteroides())
Observable('Pop_tot', Bacteroides())
Observable('Metab_tot', Inulin() + Glucose() + Lactose() + Fructo() + ChondSulf() + Lactate())

Observable('Inulin_tot', Inulin())
Observable('Glucose_tot', Glucose())
Observable('Lactose_tot', Lactose())
Observable('Fructo_tot', Fructo())
Observable('ChondSulf_tot', ChondSulf())
Observable('Lactate_tot', Lactate())

obs_Bact_0_100 = [Observable('Bact_E%d' % i, Bacteroides(energy='_%d' % i)) for i in range(n_levels)]
obs_Bact_gt_100 = [Observable('Bact_E%d' % i, Bacteroides(energy='_%d' % i)) for i in range(n_levels, n_levels +
                                                                                            deltaE_50)]
# print(model.observables)
# quit()

obs_to_plot = ['Bact_tot']

Hungry_bact = [Bacteroides(energy='_%d' % i) for i in range(hung_threshold + 1)]

Hungry_obs = Hungry_bact[0]
for i in range(1, len(Hungry_bact)):
    Hungry_obs += Hungry_bact[i]
Observable('Hungry_tot', Hungry_obs)

# metabolite production rules
Parameter('k_Inulin_prod', 10/t_step) #10/t_step)
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

# bacterial consumption rules
#  k_hungry = 1/(N*M), N = total number of hungry bacteria, M = total number of metabolites
# (Nb/N * 1/Nb) * (Mi/M * 1/Mi) = prob of selecting any given bacteroide with any given inulin

n_Bact_Inulin = 1  # energy increase in gutlogo code
# Parameter('k_Bact_Inulin_Basal', 0)
# Parameter('k_Bact_Inulin_Hungry', 10)
Expression('k_Bact_Inulin_Hungry', Piecewise((0, (Metab_tot < 1) | (Hungry_tot < 1)),
                                             (n_Bact_Inulin/(Metab_tot*Hungry_tot)/t_step, True)))
# [Rule('BactEatInulin_Basal_%d_%d' % (i, i+deltaE_25),
#       Bacteroides(energy='_%d' % i) + Inulin() >> Bacteroides(energy='_%d' % (i+deltaE_25)), k_Bact_Inulin_Basal)
#  for i in range(hung_threshold, n_levels - 1)]
[Rule('BactEatInulin_Hungry_%d_%d' % (i, i+deltaE_25),
      Bacteroides(energy='_%d' % i) + Inulin() >> Bacteroides(energy='_%d' % (i+deltaE_25)),
      k_Bact_Inulin_Hungry)
 for i in range(hung_threshold + 1)]

n_Bact_Glucose = 1
# Parameter('k_Bact_Glucose_Basal', 0)
# Parameter('k_Bact_Glucose_Hungry', 10)
Expression('k_Bact_Glucose_Hungry', Piecewise((0, (Metab_tot < 1) | (Hungry_tot < 1)),
                                              (n_Bact_Glucose/(Metab_tot*Hungry_tot)/t_step, True)))
# [Rule('BactEatGlucose_Basal_%d_%d' % (i,i+1), Bacteroides(energy='_%d' % i) + Glucose() >> Bacteroides(energy= '_%d' % (i+deltaE_50)), k_Bact_Glucose_Basal)
#  for i in range(hung_threshold, n_levels - 1)]
[Rule('BactEatGlucose_Hungry_%d_%d' % (i, i+deltaE_50),
      Bacteroides(energy='_%d' % i) + Glucose() >> Bacteroides(energy='_%d' % (i+deltaE_50)),
      k_Bact_Glucose_Hungry)
 for i in range(hung_threshold + 1)]

n_Bact_Lactose = 1
# Parameter('k_Bact_Lactose_Basal', 0)
# Parameter('k_Bact_Lactose_Hungry', 10)
Expression('k_Bact_Lactose_Hungry', Piecewise((0, (Metab_tot < 1) | (Hungry_tot < 1)),
                                              (n_Bact_Lactose/(Metab_tot*Hungry_tot)/t_step, True)))
# [Rule('BactEatLactose_Basal_%d_%d' % (i,i+1), Bacteroides(energy='_%d' % i) + Lactose() >> Bacteroides(energy= '_%d' % (i+deltaE_25)), k_Bact_Lactose_Basal)
#  for i in range(hung_threshold, n_levels - 1)]
[Rule('BactEatLactose_Hungry_%d_%d' % (i, i+deltaE_25),
      Bacteroides(energy='_%d' % i) + Lactose() >> Bacteroides(energy= '_%d' % (i+deltaE_25)),
      k_Bact_Lactose_Hungry)
 for i in range(hung_threshold + 1)]

n_Bact_Fructo = 1
# Parameter('k_Bact_Fructo_Basal', 0)
# Parameter('k_Bact_Fructo_Hungry', 10)
Expression('k_Bact_Fructo_Hungry', Piecewise((0, (Metab_tot < 1) | (Hungry_tot < 1)),
                                             (n_Bact_Fructo/(Metab_tot*Hungry_tot)/t_step, True)))
# [Rule('BactEatFructo_Basal_%d_%d' % (i,i+1), Bacteroides(energy='_%d' % i) + Fructo() >> Bacteroides(energy= '_%d' % (i+deltaE_25)), k_Bact_Fructo_Basal)
#  for i in range(hung_threshold, n_levels - 1)]
[Rule('BactEatFructo_Hungry_%d_%d' % (i, i+deltaE_25),
      Bacteroides(energy='_%d' % i) + Fructo() >> Bacteroides(energy='_%d' % (i+deltaE_25)),
      k_Bact_Fructo_Hungry)
 for i in range(hung_threshold + 1)]

# energy loss rules
Parameter('k_energy_loss', 1/(1440/n_levels * t_step))

[Rule('Bact_Loss_%d_%d' % (i, i-1),
      Bacteroides(energy='_%d' % i) >> Bacteroides(energy='_%d' % (i-1)), k_energy_loss)
 for i in range(1, n_levels + deltaE_50)]

# print(model.rules)
# quit()

# death: no mention of death, but any agent that moves past the final horizontal element is removed from the simulation
# death: bacteria move from element to element as a function of velocity of the flow field (displacement per timestep)
# death: probabilistic function controls likelyhood of getting stuck in the membrane. Probability of escaping is 10% (per some time step?) and 5% to become perm embedded
# since our program does not have a spacial element we will have to add a death element instead
# DT: all genera have same doubling times and they are all constant: 330 timesteps I believe
# MC: probability of consuming metabolite is function of metabolite concentrate and number of other hungry (energy<8) bacteria on the element
# MC: the function is not shown but I believe a good starting point would be (# of metabolites /# of hungry bacteria) per time step

div_threshold = int(n_levels/2)
td_Bact = 330*t_step

Parameter('k_Bact_division', (np.log(2)/td_Bact))
#todo; account for number of cells during division like in gutlogo code below
#todo; play around in gutlogo code on netogo; isolate the cause of growth at 350 timesteps
#todo; double check that their stuck rules match their paper (stuck --> stuck + unstuck)
"""
to reproduceBact
;; reproduce the chosen turtle
  if energy > 50 and count turtles-here < 1000[ ;;turtles-here check to model space limit
    let tmp (energy / 2 )
    set energy (tmp) ;; parent's energy is halved
    hatch 1 [
      rt random-float 360
      set energy tmp ;; child gets half of parent's energy
      set isSeed false
      set isStuck false
	    set age 0
    ]
  ]
end
"""
# bacteria division rules
[Rule('Bact_divides_%d' % i, Bacteroides(energy='_%d' % i) >>
      Bacteroides(energy='_%d' % (i/2)) + Bacteroides(energy='_%d' % ((i + 1)/2), stuck='u'), k_Bact_division)
 for i in range(div_threshold, n_levels + deltaE_50)]

# death rules
Parameter('k_Death', 1/t_step)
Rule('Bact_death', Bacteroides(energy='_0') >> None, k_Death)

# removal rules

flow_rate = 2.22  # cm/min
# 7.98 cm is the length of an element; there are 100 total elements; 60 seconds per minute
Parameter('k_Bact_unstuck_removed', flow_rate / (7.98*100) / 60)  # /s

Rule('Bact_unstuck_removal', Bacteroides(stuck='u') >> None, k_Bact_unstuck_removed)

# stuck rules
maxStuckChance = 0.5  # probability
midStuckConc = 10  # concentration
lowStuckBound = 0.02  # probability

Expression('k_Bact_stuck', Piecewise(
 (0, maxStuckChance * midStuckConc / (midStuckConc + Pop_tot) / t_step < lowStuckBound),
 (maxStuckChance * midStuckConc / (midStuckConc + Pop_tot) / t_step, True)))
Parameter('k_Bact_unstuck', 0.1/t_step)
Parameter('k_Bact_permstuck', 0.05/t_step)

Rule('Bact_stuck', Bacteroides(stuck='u') | Bacteroides(stuck='s'), k_Bact_stuck, k_Bact_unstuck)
Rule('Bact_permstuck', Bacteroides(stuck='s') >> Bacteroides(stuck='p'), k_Bact_permstuck)


# TODO: turn on inflow rules
# bacterial creation rules
Parameter('k_Bact_creation', 0)
Rule('Bact_creation', None >> Bacteroides(energy='_%d' % (n_levels - 1), stuck='u'), k_Bact_creation)

print(model.rules)
print()
for ic in model.initial_conditions:
    print(ic)
# quit()

# simulations
n_steps = 1000 # 5000
t_span = np.linspace(0, n_steps*t_step, int(n_steps)*1+1)
sim = ScipyOdeSimulator(model, t_span, verbose=True)
result = sim.run()

for obs in obs_to_plot:
    plt.plot(range(n_steps + 1), result.observables[obs], label=obs)
    plt.xlabel('time (number of time steps)')
    plt.ylabel('# of cells')
    # plt.yscale('log', base = 10)
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

# plot metabolites
plt.figure()
for obs in [Inulin_tot, Glucose_tot, Lactose_tot, Fructo_tot, ChondSulf_tot, Lactate_tot]:
    plt.plot(range(n_steps + 1), result.observables[obs.name], label=obs.name)
plt.xlabel('time (number of time steps)')
plt.ylabel('# of molecules')
plt.legend(loc=0)

# plot bacterial species at all energy levels
for obs_0_100 in [obs_Bact_0_100]:
    plt.figure()
    for obs in obs_0_100:
        plt.plot(range(n_steps + 1), result.observables[obs.name], label=obs.name)
    plt.xlabel('time (number of time steps)')
    plt.ylabel('# of cells')
    # plt.yscale('log', base=10)
    plt.legend(loc=0)

for obs_gt_100 in [obs_Bact_gt_100]:
    plt.figure()
    for obs in obs_gt_100:
        plt.plot(range(n_steps + 1), result.observables[obs.name], '--', label=obs.name)
    plt.xlabel('time (number of time steps)')
    plt.ylabel('# of cells')
    plt.yscale('log', base=10)
    plt.legend(loc=0)

plt.show()
