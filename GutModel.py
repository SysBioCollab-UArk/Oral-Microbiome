from pysb import *
from pysb.simulator import ScipyOdeSimulator
import numpy as np
import matplotlib.pyplot as plt
from sympy import Piecewise
from ReadGutlogo import read_Gutlogo

Model()

# PARAMETER SETTINGS FROM GUTLOGO

# Absorption represents the percentage of metabolites that pass through the gut into the bloodstream.
Parameter('absorption', 0)
# A 'reserveFraction' of greater than zero will activate the reserve metabolite module of the simulation. This module
# inhibits the bacteria from accessing a portion of the metabolites. The portion of metabolite available is a function
# of position down the gut. Therefore, a colony further down the gut would have access to a larger portion of the
# remaining bacteria than a colony closer to the start.
Parameter('reservefraction', 0)

# initial bacterial numbers (#)
Parameter("initnumclosts", 921)
Parameter("initnumdesulfos", 70)
Parameter("initnumbifidos", 23562)
Parameter("initnumbacteroides", 5490)

# inflow of bacteria: # of bacteria that enter the system every 'tickinflow' steps
Parameter("inconcbifidos", 0)
Parameter("inconcclosts", 0)
Parameter("inconcbacteroides", 0)
Parameter("inconcdesulfos", 0)

# Bacterial doubling times (steps)
Parameter("bacteroiddoub", 330)
Parameter("clostdoub", 330)
Parameter("bifidodoub", 330)
Parameter("desulfodoub", 330)

# inflow of metabolites per timestep (#/step)
Parameter("inflowfo", 25)
Parameter("inflowcs", 0.1)
Parameter("inflowglucose", 30)
Parameter("inflowlactose", 15)
Parameter("inflowinulin", 10)
Parameter("inflowlactate", 0)

# number of lactate molecules produced by bifido (#/step)
Parameter("bifido_lactate_production", 0.005)

# parameters of sticking and unsticking of bacteria in the gut lining
Parameter("seedpercent", 5.0)  # probability/step
Parameter("seedchance", 5.0)  # probability/step
Parameter("unstuckchance", 10)  # probability/step
Parameter("lowstuckbound", 2)  # probability/step
Parameter("maxstuckchance", 50)  # probability/step
Parameter("midstuckconc", 10)  # # of cells

# Fluid flow rate throughout the gut (patches/step)
Parameter("flowdist", 0.28)
# How often bacteria enter the gut (steps)
Parameter("tickinflow", 480)

# Define PySB model elements

# TODO: Not sure why we're using seconds as the time unit. Each step = 1 minute. Could we use minutes as the time unit?
t_step = 60  # seconds

# TODO: add carrying capacity rules - STILL NEED TO DO THIS (06/08/24)
# yes, stuck bacteria can reproduce, and offspring bacteria is unstuck
# parameters to find: hungry consumption rate, basal consumption rate, lactate production, stuck rates

# BACTERIA
n_levels = 10  # TODO: eventually want to increase this to 100 levels
deltaE_50 = int(np.ceil(n_levels / 2))  # = 5 if n_levels = 10
deltaE_25 = int(np.ceil(deltaE_50 / 2))  # = 3 if n_levels = 10
Monomer('Bacteroides', ['energy', 'stuck'],
        {'energy': ['_%d' % i for i in range(0, n_levels + deltaE_50)], 'stuck': ['u', 's', 'p']})
Monomer('Clostridium', ['energy', 'stuck'],
        {'energy': ['_%d' % i for i in range(0, n_levels + deltaE_50)], 'stuck': ['u', 's', 'p']})
Monomer('Bifidobacterium', ['energy', 'stuck'],
        {'energy': ['_%d' % i for i in range(0, n_levels + deltaE_50)], 'stuck': ['u', 's', 'p']})
Monomer('Desulfobrivio', ['energy', 'stuck'],
        {'energy': ['_%d' % i for i in range(0, n_levels + deltaE_50)], 'stuck': ['u', 's', 'p']})

# total # of initial bacteria (unstuck and stuck)

Expression('Bact_unstuck_0', initnumbacteroides * (1-seedpercent / 100))
Expression('Clost_unstuck_0', initnumclosts * (1-seedpercent / 100))
Expression('Bifido_unstuck_0', initnumbifidos * (1-seedpercent / 100))
Expression('Desulfo_unstuck_0', initnumdesulfos * (1-seedpercent / 100))

Initial(Bacteroides(energy='_%d' % (n_levels - 1), stuck='u'), Bact_unstuck_0)
Initial(Clostridium(energy='_%d' % (n_levels - 1), stuck='u'), Clost_unstuck_0)
Initial(Bifidobacterium(energy='_%d' % (n_levels - 1), stuck='u'), Bifido_unstuck_0)
Initial(Desulfobrivio(energy='_%d' % (n_levels - 1), stuck='u'), Desulfo_unstuck_0)

Expression('Bact_permstuck_0', initnumbacteroides - Bact_unstuck_0)
Expression('Clost_permstuck_0', initnumclosts - Clost_unstuck_0)
Expression('Bifido_permstuck_0', initnumbifidos - Bifido_unstuck_0)
Expression('Desulfo_permstuck_0', initnumdesulfos - Desulfo_unstuck_0)

Initial(Bacteroides(energy='_%d' % (n_levels - 1), stuck='p'), Bact_permstuck_0)
Initial(Clostridium(energy='_%d' % (n_levels - 1), stuck='p'), Clost_permstuck_0)
Initial(Bifidobacterium(energy='_%d' % (n_levels - 1), stuck='p'), Bifido_permstuck_0)
Initial(Desulfobrivio(energy='_%d' % (n_levels - 1), stuck='p'), Desulfo_permstuck_0)

Observable('Bact_tot', Bacteroides())
Observable('Clost_tot', Clostridium())
Observable('Bifido_tot', Bifidobacterium())
Observable('Desulfo_tot', Desulfobrivio())
Observable('Pop_tot', Bacteroides() + Clostridium() + Bifidobacterium() + Desulfobrivio())

# Temporary #########
# Observable('Bact_E150', Bacteroides(energy='_%d' % (n_levels + deltaE_50 - 1)))
# Observable('Bact_E125', Bacteroides(energy='_%d' % (n_levels + deltaE_25 - 1)))
# Observable('Bact_E100', Bacteroides(energy='_%d' % (n_levels - 1)))
# Observable('Bact_E90', Bacteroides(energy='_%d' % (n_levels - 2)))
# Observable('Bact_E80', Bacteroides(energy='_%d' % (n_levels - 3)))
# Observable('Bact_E70', Bacteroides(energy='_%d' % (n_levels - 4)))
# Observable('Bact_E0', Bacteroides(energy='_0'))

obs_Bact_0_100 = [Observable('Bact_E%d' % i, Bacteroides(energy='_%d' % i)) for i in range(n_levels)]
obs_Bact_gt_100 = [Observable('Bact_E%d' % i, Bacteroides(energy='_%d' % i)) for i in range(n_levels, n_levels + deltaE_50)]
obs_Clost_0_100 = [Observable('Clost_E%d' % i, Clostridium(energy='_%d' % i)) for i in range(n_levels)]
obs_Clost_gt_100 = [Observable('Clost_E%d' % i, Clostridium(energy='_%d' % i)) for i in range(n_levels, n_levels + deltaE_50)]
obs_Bifido_0_100 = [Observable('Bifido_E%d' % i, Bifidobacterium(energy='_%d' % i)) for i in range(n_levels)]
obs_Bifido_gt_100 = [Observable('Bifido_E%d' % i, Bifidobacterium(energy='_%d' % i)) for i in range(n_levels, n_levels + deltaE_50)]
obs_Desulfo_0_100 = [Observable('Desulfo_E%d' % i, Desulfobrivio(energy='_%d' % i)) for i in range(n_levels)]
obs_Desulfo_gt_100 = [Observable('Desulfo_E%d' % i, Desulfobrivio(energy='_%d' % i)) for i in range(n_levels, n_levels + deltaE_50)]
#########

# METABOLITES

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

Observable('Inulin_tot', Inulin())
Observable('Glucose_tot', Glucose())
Observable('Lactose_tot', Lactose())
Observable('Fructo_tot', Fructo())
Observable('ChondSulf_tot', ChondSulf())
Observable('Lactate_tot', Lactate())
Observable('Metab_tot', Inulin() + Glucose() + Lactose() + Fructo() + ChondSulf() + Lactate())

# RULES

# metabolite production rules
Expression('k_inflowinulin', inflowinulin/t_step)  # 10/t_step)
Expression('k_inflowglucose', inflowglucose/t_step) #30
Expression('k_inflowlactose', inflowlactose/t_step) #15
Expression('k_inflowfo', inflowfo/t_step) #25
Expression('k_inflowcs', inflowcs/t_step) #0.1
Expression('k_inflowlactate', inflowlactate/t_step) #0

Rule('Inulin_prod', None >> Inulin(), k_inflowinulin)
Rule('Glucose_prod', None >> Glucose(), k_inflowglucose)
Rule('Lactose_prod', None >> Lactose(), k_inflowlactose)
Rule('Fructo_prod', None >> Fructo(), k_inflowfo)
Rule('ChondSulf_prod', None >> ChondSulf(), k_inflowcs)
Rule('Lactate_prod', None >> Lactate(), k_inflowlactate)

# BACTERIAL EATING RULES
hung_threshold = int(0.8 * n_levels - 1)
# k_hungry = 1/(N*M), N = total number of hungry bacteria, M = total number of metabolites
# (Nb/N * 1/Nb) * (Mi/M * 1/Mi) = prob of selecting any given bacteroide with any given inulin

Hungry_bact = [Bacteroides(energy='_%d' % i) for i in range(hung_threshold + 1)] + \
              [Clostridium(energy='_%d' % i) for i in range(hung_threshold + 1)] + \
              [Bifidobacterium(energy='_%d' % i) for i in range(hung_threshold + 1)] + \
              [Desulfobrivio(energy='_%d' % i) for i in range(hung_threshold + 1)]
Hungry_obs = Hungry_bact[0]
for i in range(1, len(Hungry_bact)):
    Hungry_obs += Hungry_bact[i]

Observable('Hungry_tot', Hungry_obs)

deltaE_Bact_Inulin = deltaE_25  # energy increase in GutLogo code
Expression('k_Bact_Inulin_Hungry', Piecewise((0, (Metab_tot < 1) | (Hungry_tot < 1)),
                                             (1/(Metab_tot*Hungry_tot)/t_step, True)))
[Rule('BactEatInulin_Hungry_%d_%d' % (i, i+deltaE_Bact_Inulin),
      Bacteroides(energy='_%d' % i) + Inulin() >> Bacteroides(energy='_%d' % (i+deltaE_Bact_Inulin)),
      k_Bact_Inulin_Hungry) for i in range(hung_threshold + 1)]

deltaE_Clost_Inulin = deltaE_25
Expression('k_Clost_Inulin_Hungry', Piecewise((0, (Metab_tot < 1) | (Hungry_tot < 1)),
                                              (1/(Metab_tot*Hungry_tot)/t_step, True)))
[Rule('ClostEatInulin_Hungry_%d_%d' % (i, i+deltaE_Clost_Inulin),
      Clostridium(energy='_%d' % i) + Inulin() >> Clostridium(energy= '_%d' % (i+deltaE_Clost_Inulin)),
      k_Clost_Inulin_Hungry) for i in range(hung_threshold + 1)]

deltaE_Desulfo_Choldsulf = deltaE_25
Expression('k_Desulfo_ChondSulf_Hungry', Piecewise((0, (Metab_tot < 1) | (Hungry_tot < 1)),
                                                   (1/(Metab_tot*Hungry_tot)/t_step, True)))
[Rule('DesulfoEatChond_Hungry_%d_%d' % (i, i+deltaE_Desulfo_Choldsulf),
      Desulfobrivio(energy='_%d' % i) + ChondSulf() >> Desulfobrivio(energy='_%d' % (i+deltaE_Desulfo_Choldsulf)),
      k_Desulfo_ChondSulf_Hungry) for i in range(hung_threshold + 1)]

deltaE_Desulfo_Lactate = deltaE_50
Expression('k_Desulfo_Lactate_Hungry', Piecewise((0, (Metab_tot < 1) | (Hungry_tot < 1)),
                                                 (1/(Metab_tot*Hungry_tot)/t_step, True)))
[Rule('DesulfoEatLactate_Hungry_%d_%d' % (i, i+deltaE_Desulfo_Lactate),
      Desulfobrivio(energy='_%d' % i) + Lactate() >> Desulfobrivio(energy='_%d' % (i+deltaE_Desulfo_Lactate)),
      k_Desulfo_Lactate_Hungry) for i in range(hung_threshold + 1)]

deltaE_Bact_Glucose = deltaE_50
Expression('k_Bact_Glucose_Hungry', Piecewise((0, (Metab_tot < 1) | (Hungry_tot < 1)),
                                              (1/(Metab_tot*Hungry_tot)/t_step, True)))
[Rule('BactEatGlucose_Hungry_%d_%d' % (i, i+deltaE_Bact_Glucose),
      Bacteroides(energy='_%d' % i) + Glucose() >> Bacteroides(energy='_%d' % (i+deltaE_Bact_Glucose)),
      k_Bact_Glucose_Hungry) for i in range(hung_threshold + 1)]

deltaE_Bact_Lactose = deltaE_25
Expression('k_Bact_Lactose_Hungry', Piecewise((0, (Metab_tot < 1) | (Hungry_tot < 1)),
                                              (1/(Metab_tot*Hungry_tot)/t_step, True)))
[Rule('BactEatLactose_Hungry_%d_%d' % (i, i+deltaE_Bact_Lactose),
      Bacteroides(energy='_%d' % i) + Lactose() >> Bacteroides(energy= '_%d' % (i+deltaE_Bact_Lactose)),
      k_Bact_Lactose_Hungry) for i in range(hung_threshold + 1)]

deltaE_Bact_Fructo = deltaE_25
Expression('k_Bact_Fructo_Hungry', Piecewise((0, (Metab_tot < 1) | (Hungry_tot < 1)),
                                             (1/(Metab_tot*Hungry_tot)/t_step, True)))
[Rule('BactEatFructo_Hungry_%d_%d' % (i, i+deltaE_Bact_Fructo),
      Bacteroides(energy='_%d' % i) + Fructo() >> Bacteroides(energy='_%d' % (i+deltaE_Bact_Fructo)),
      k_Bact_Fructo_Hungry) for i in range(hung_threshold + 1)]

deltaE_Clost_Glucose = deltaE_50
Expression('k_Clost_Glucose_Hungry', Piecewise((0, (Metab_tot < 1) | (Hungry_tot < 1)),
                                               (1/(Metab_tot*Hungry_tot)/t_step, True)))
[Rule('ClostEatGlucose_Hungry_%d_%d' % (i, i+deltaE_Clost_Glucose),
      Clostridium(energy='_%d' % i) + Glucose() >> Clostridium(energy='_%d' % (i+deltaE_Clost_Glucose)),
      k_Clost_Glucose_Hungry) for i in range(hung_threshold + 1)]

deltaE_Clost_Lactose = deltaE_25
Expression('k_Clost_Lactose_Hungry', Piecewise((0, (Metab_tot < 1) | (Hungry_tot < 1)),
                                               (1/(Metab_tot*Hungry_tot)/t_step, True)))
[Rule('ClostEatLactose_Hungry_%d_%d' % (i, i+deltaE_Clost_Lactose),
      Clostridium(energy='_%d' % i) + Lactose() >> Clostridium(energy='_%d' % (i+deltaE_Clost_Lactose)),
      k_Clost_Lactose_Hungry) for i in range(hung_threshold + 1)]

deltaE_Clost_Fructo = deltaE_25
Expression('k_Clost_Fructo_Hungry', Piecewise((0, (Metab_tot < 1) | (Hungry_tot < 1)),
                                              (1/(Metab_tot*Hungry_tot)/t_step, True)))
[Rule('ClostEatFructo_Hungry_%d_%d' % (i, i+deltaE_Clost_Fructo),
      Clostridium(energy='_%d' % i) + Fructo() >> Clostridium(energy='_%d' % (i+deltaE_Clost_Fructo)),
      k_Clost_Fructo_Hungry) for i in range(hung_threshold + 1)]

deltaE_Bifido_Glucose = deltaE_25
Expression('k_Bifido_Glucose_Hungry', Piecewise((0, (Metab_tot < 1) | (Hungry_tot < 1)),
                                                (1/(Metab_tot*Hungry_tot)/t_step, True)))
[Rule('BifidoEatGlucose_Hungry_%d_%d' % (i, i+deltaE_Bifido_Glucose),
      Bifidobacterium(energy='_%d' % i) + Glucose() >> Bifidobacterium(energy='_%d' % (i+deltaE_Bifido_Glucose)),
      k_Bifido_Glucose_Hungry) for i in range(hung_threshold + 1)]

deltaE_Bifido_Lactose = deltaE_50
Expression('k_Bifido_Lactose_Hungry', Piecewise((0, (Metab_tot < 1) | (Hungry_tot < 1)),
                                                (1/(Metab_tot*Hungry_tot)/t_step, True)))
[Rule('BifidoEatLactose_Hungry_%d_%d' % (i, i+deltaE_Bifido_Lactose),
      Bifidobacterium(energy='_%d' % i) + Lactose() >> Bifidobacterium(energy='_%d' % (i+deltaE_Bifido_Lactose)),
      k_Bifido_Lactose_Hungry) for i in range(hung_threshold + 1)]

deltaE_Bifido_Fructo = deltaE_50
Expression('k_Bifido_Fructo_Hungry', Piecewise((0, (Metab_tot < 1) | (Hungry_tot < 1)),
                                               (1/(Metab_tot*Hungry_tot)/t_step, True)))
[Rule('BifidoEatFructo_Hungry_%d_%d' % (i, i+deltaE_Bifido_Fructo),
      Bifidobacterium(energy='_%d' % i) + Fructo() >> Bifidobacterium(energy='_%d' % (i+deltaE_Bifido_Fructo)),
      k_Bifido_Fructo_Hungry) for i in range(hung_threshold + 1)]

# lactate production rule
Expression('k_bifido_lactate_production', bifido_lactate_production/t_step)
Rule('BifidoMakesLactate', Bifidobacterium() >> Bifidobacterium() + Lactate(), k_bifido_lactate_production)

# energy loss rules
Parameter('k_energy_loss', 1/(1440/n_levels * t_step))  # 1440 minutes = 24 hours, time cells can survive w/o eating

[Rule('Bact_Loss_%d_%d' % (i, i-1),
      Bacteroides(energy='_%d' % i) >> Bacteroides(energy='_%d' % (i-1)), k_energy_loss)
 for i in range(1, n_levels + deltaE_50)]

[Rule('Clost_Loss_%d_%d' % (i, i-1),
      Clostridium(energy='_%d' % i) >> Clostridium(energy='_%d' % (i-1)), k_energy_loss)
 for i in range(1, n_levels + deltaE_50)]

[Rule('Desulfo_Loss_%d_%d' % (i, i-1),
      Desulfobrivio(energy='_%d' % i) >> Desulfobrivio(energy='_%d' % (i-1)), k_energy_loss)
 for i in range(1, n_levels + deltaE_50)]

[Rule('Bifido_Loss_%d_%d' % (i, i-1),
      Bifidobacterium(energy='_%d' % i) >> Bifidobacterium(energy='_%d' % (i-1)), k_energy_loss)
 for i in range(1, n_levels + deltaE_50)]

# print(model.rules)
# quit()

# death: no mention of death, but any agent that moves past the final horizontal element is removed from the simulation
# death: bacteria move from element to element as a function of velocity of the flow field (displacement per timestep)
# death: probabilistic function controls likelihood of getting stuck in the membrane. Probability of escaping is 10%
#        (per some time step?) and 5% to become perm embedded
# since our program does not have a spacial element we will have to add a death element instead
# DT: all genera have same doubling times and they are all constant: 330 timesteps I believe
# MC: probability of consuming metabolite is function of metabolite concentrate and number of other hungry (energy<8)
#     bacteria on the element
# MC: the function is not shown but I believe a good starting point would be (# of metabolites /# of hungry bacteria)
#     per time step

div_threshold = int(n_levels/2)
# bacteroiddoub = 330*t_step  # doubling time
# clostdoub = 330*t_step
# desulfodoub = 330*t_step
# bifidodoub = 330*t_step
Expression('k_Bact_division', np.log(2)/(bacteroiddoub * t_step))  # division rate: X = Xo * exp(kt), X/Xo = 2, t = td
Expression('k_Clost_division', np.log(2)/(clostdoub * t_step))
Expression('k_Desulfo_division', np.log(2)/(desulfodoub * t_step))
Expression('k_Bifido_division', np.log(2)/(bifidodoub * t_step))
# todo; account for number of cells during division like in gutlogo code below
# todo; play around in gutlogo code on netogo; isolate the cause of growth at 350 timesteps
# todo; double check that their stuck rules match their paper (stuck --> stuck + unstuck)
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
      Bacteroides(energy='_%d' % (i/2)) + Bacteroides(energy='_%d' % ((i+1)/2), stuck='u'), k_Bact_division)
 for i in range(div_threshold, n_levels + deltaE_50)]

[Rule('Clost_divides_%d' % i, Clostridium(energy='_%d' % i) >>
      Clostridium(energy='_%d' % (i/2)) + Clostridium(energy='_%d' % ((i+1)/2), stuck='u'), k_Clost_division)
 for i in range(div_threshold, n_levels + deltaE_50)]

[Rule('Desulfo_divides_%d' % i, Desulfobrivio(energy='_%d' % i) >>
      Desulfobrivio(energy='_%d' % (i/2)) + Desulfobrivio(energy='_%d' % ((i+1)/2), stuck='u'), k_Desulfo_division)
 for i in range(div_threshold, n_levels + deltaE_50)]

[Rule('Bifido_divides_%d' % i, Bifidobacterium(energy='_%d' % i) >>
      Bifidobacterium(energy='_%d' % (i/2)) + Bifidobacterium(energy='_%d' % ((i+1)/2), stuck='u'), k_Bifido_division)
 for i in range(div_threshold, n_levels + deltaE_50)]

# death by energy
Parameter('k_Death_Energy', 1e6)  # death should be instantaneous when energy = 0
Rule('Bact_death_energy', Bacteroides(energy='_0') >> None, k_Death_Energy)
Rule('Clost_death_energy', Clostridium(energy='_0') >> None, k_Death_Energy)
Rule('Desulfo_death_energy', Desulfobrivio(energy='_0') >> None, k_Death_Energy)
Rule('Bifido_death_energy', Bifidobacterium(energy='_0') >> None, k_Death_Energy)

# death by age
# TODO: Don't think cells can die by old age in GutLogo. ** Double check this **
Parameter('k_Death_Age', 1/(1000 * t_step)) # age is a counter in GutLogo code going from 0 to 1000
Rule('Bact_death_age', Bacteroides() >> None, k_Death_Age)
Rule('Clost_death_age', Clostridium() >> None, k_Death_Age)
Rule('Desulfo_death_age', Desulfobrivio() >> None, k_Death_Age)
Rule('Bifido_death_age', Bifidobacterium() >> None, k_Death_Age)

# removal rules (by flow)

# flowdist in patches/min #2.22
# 7.98 cm is the length of an element; there are 100 total elements; 60 seconds per minute
Expression('k_Bact_unstuck_removed', flowdist / 100 / 60)  # /s
Expression('k_Clost_unstuck_removed', flowdist / 100 / 60)  # /s
Expression('k_Desulfo_unstuck_removed', flowdist / 100 / 60)  # /s
Expression('k_Bifido_unstuck_removed', flowdist / 100 / 60)  # /s

Rule('Bact_unstuck_removal', Bacteroides(stuck='u') >> None, k_Bact_unstuck_removed)
Rule('Clost_unstuck_removal', Clostridium(stuck='u') >> None, k_Clost_unstuck_removed)
Rule('Desulfo_unstuck_removal', Desulfobrivio(stuck='u') >> None, k_Desulfo_unstuck_removed)
Rule('Bifido_unstuck_removal', Bifidobacterium(stuck='u') >> None, k_Bifido_unstuck_removed)

# stuck rules
Expression('k_stuck',
           Piecewise(
               (0, maxstuckchance * (1 - Pop_tot / (midstuckconc + Pop_tot)) < lowstuckbound),
               (maxstuckchance / 100 * (1 - Pop_tot / (midstuckconc + Pop_tot)) / t_step, True)
           )
)
Expression('k_unstuckchance', unstuckchance / 100 / t_step)  # probability per step #0.1
Expression('k_seedchance', seedchance / 100 / t_step)  # probability per step #0.05

Rule('Bact_stuck', Bacteroides(stuck='u') | Bacteroides(stuck='s'), k_stuck, k_unstuckchance)
Rule('Bact_permstuck', Bacteroides(stuck='s') >> Bacteroides(stuck='p'), k_seedchance)

Rule('Clost_stuck', Clostridium(stuck='u') | Clostridium(stuck='s'), k_stuck, k_unstuckchance)
Rule('Clost_permstuck', Clostridium(stuck='s') >> Clostridium(stuck='p'), k_seedchance)

Rule('Desulfo_stuck', Desulfobrivio(stuck='u') | Desulfobrivio(stuck='s'), k_stuck, k_unstuckchance)
Rule('Desulfo_permstuck', Desulfobrivio(stuck='s') >> Desulfobrivio(stuck='p'), k_seedchance)

Rule('Bifido_stuck', Bifidobacterium(stuck='u') | Bifidobacterium(stuck='s'), k_stuck, k_unstuckchance)
Rule('Bifido_permstuck', Bifidobacterium(stuck='s') >> Bifidobacterium(stuck='p'), k_seedchance)

# BACTERIA INPUT RULES (COMING INTO THE GUT)

# TODO: NOTE - This section has been cleaned up. Follow this example to clean up the rest of the code. In particular,
#  we don't want to use the GutLogo parameters directly in the rules unless they are the rate parameter exactly. If
#  they are only part of the rate parameter, crete an expression that can be turned into a parameter object later on

# Bacteria inflow rates (#/s)
Expression("k_inconcbacteroides", inconcbacteroides / (tickinflow * t_step))
Expression("k_inconcclosts", inconcclosts / (tickinflow * t_step))
Expression("k_inconcdesulfos", inconcdesulfos / (tickinflow * t_step))
Expression("k_inconcbifidos", inconcbifidos / (tickinflow * t_step))

Rule('Bact_inflow', None >> Bacteroides(energy='_%d' % (n_levels - 1), stuck='u'), k_inconcbacteroides)
Rule('Clost_inflow', None >> Clostridium(energy='_%d' % (n_levels - 1), stuck='u'), k_inconcclosts)
Rule('Desulfo_inflow', None >> Desulfobrivio(energy='_%d' % (n_levels - 1), stuck='u'), k_inconcdesulfos)
Rule('Bifido_inflow', None >> Bifidobacterium(energy='_%d' % (n_levels - 1), stuck='u'), k_inconcbifidos)

# print(model.parameters_rules())
# print()
# print(model.expressions)
# print()
# print(model.rules)
# print()
# for ic in model.initial_conditions:
#     print(ic)
# quit()
if __name__ == '__main__':
    print(model.parameters)
    #quit()
    # simulations
    n_steps = 1000  # 5000
    t_span = np.linspace(0, n_steps*t_step, int(n_steps)*1+1)
    sim = ScipyOdeSimulator(model, t_span, verbose=True)
    result = sim.run()

    obs_to_plot = ['Bact_tot', 'Clost_tot', 'Bifido_tot', 'Desulfo_tot']

    for obs in obs_to_plot:
        plt.plot(range(n_steps + 1), result.observables[obs], label=obs)
        plt.xlabel('time (number of time steps)')
        plt.ylabel('# of cells')
        # plt.yscale('log', base = 10)
    plt.legend(loc=0)
    plt.tight_layout()

    # print(result.observables['Bact_tot'])
    # print(n_steps*t_step)
    # print(t_span)
    # print(result.expressions['k_Bact_Inulin_Hungry'])
    # print(result.observables['Hungry_tot'])
    # print(result.observables['Metab_tot'])
    # print(result.observables['Pop_tot'])

    # plot GutLogo simulations
    plt.figure()
    time, cell_counts, settings = read_Gutlogo('GutLogoPopulations_long.csv')
    for species in cell_counts.keys():
        plt.plot(time, cell_counts[species], label=species)
    plt.xlabel('Time step (tick)')
    plt.ylabel('Population (number of agents)')
    # plt.title('Kinetics Model of Gut Microbiome')
    plt.legend(loc=0)
    plt.tight_layout()

    plt.figure()
    time, cell_counts, settings = read_Gutlogo('GutLogoPopulations_short.csv')
    for species in cell_counts.keys():
        plt.plot(time, cell_counts[species], label=species)
    plt.xlabel('Time (number of time steps)')
    plt.ylabel('Population (number of agents)')
    plt.title('Kinetics Model of Gut Microbiome')
    plt.legend(loc=0)
    plt.tight_layout()

    '''
    # plot metabolites
    plt.figure()
    for obs in [Inulin_tot, Glucose_tot, Lactose_tot, Fructo_tot, ChondSulf_tot, Lactate_tot]:
        plt.plot(range(n_steps + 1), result.observables[obs.name], label=obs.name)
    plt.xlabel('time (number of time steps)')
    plt.ylabel('# of molecules')
    plt.legend(loc=0)
    
    # plot bacterial species at all energy levels
    for obs_0_100 in [obs_Bact_0_100, obs_Clost_0_100, obs_Bifido_0_100, obs_Desulfo_0_100]:
        plt.figure()
        for obs in obs_0_100:
            plt.plot(range(n_steps + 1), result.observables[obs.name], label=obs.name)
        plt.xlabel('time (number of time steps)')
        plt.ylabel('# of cells')
        # plt.yscale('log', base=10)
        plt.legend(loc=0)
    
    for obs_gt_100 in [obs_Bact_gt_100, obs_Clost_gt_100, obs_Bifido_gt_100, obs_Desulfo_gt_100]:
        plt.figure()
        for obs in obs_gt_100:
            plt.plot(range(n_steps + 1), result.observables[obs.name], '--', label=obs.name)
        plt.xlabel('time (number of time steps)')
        plt.ylabel('# of cells')
        plt.yscale('log', base=10)
        plt.legend(loc=0)
    
    # plot desulfo E6 and E12 to compare
    plt.figure()
    plt.plot(range(n_steps + 1), result.observables['Desulfo_E6'], label='Desulfo_E6')
    plt.plot(range(n_steps + 1), result.observables['Desulfo_E12'], label='Desulfo_E12')
    plt.xlabel('time (number of time steps)')
    plt.ylabel('# of cells')
    plt.yscale('log', base=2)
    plt.ylim(bottom=2**(-11))
    plt.xlim(left=3000)
    plt.legend(loc=0)
    '''
    plt.show()

    #todo there seems to be a carrying capacity that does not line up, why?; they all die of age later; "turtles-here" what is this?
    """
    **age**: Positive integer value representing the number of ticks the bacteria has been alive for. Used to determine if
     a bacteria would reproduce on the current tick. Seed colony bacteria are given a random age from 0 to 1000.
     
     Offspring bacteria can NEVER be stuck; age starts at 0
     
     bacteria must be a certain age to reproduce
     
     incoming bacteria have a randomized age form 0 to 1000
     
     if (age mod clostDoub = 0 and age != 0)[
          reproduceBact
        ]
     set age (age + 1)
     
     bact - 8500
    """

    # todo; model the unstuck, permstuck, stuck thing--exponential function? not constant?
    # todo; look at similar literature from citations regarding basic biology

    # todo; run gutlogo again under simplistic conditions
    # (turn stuck variables off; use one bacteria type) run that-save the output. Can we reproduce that with our model, by
    # reintroducing one variable at a time

