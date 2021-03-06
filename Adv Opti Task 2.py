# Specify term (5 or 7 or 0 for all) and half (1 or 2) to solve for
solve_term = 0
solve_half = 1

T = 92 # the last time-index, corresponding to the Fri 1300-1330 timeslot
EOD_timings = [24, 48, 58, 82, 92] # last time-index of each day
HASS_block_timings = [14,15,16,17,18,19,25,26,27,72,73,74,75,76,77,89,90,91] # time-indices blocked for HASS
venue_overlap = {6:[4,5], 32:[30,31]} # e.g. 1.415/1.416: [1.415, 1.416]

#############################################################################
import gurobipy as gp
from gurobipy import GRB
import pandas as pd
import os
import warnings
import re
import time

### Read data to dataframe and create multidict #############################
#os.chdir("set your wd here")
ttdata_full = pd.read_csv('TT Data v6.csv')
ttdata_otherpillars = pd.read_csv('TT_NonESD_Data v2.csv')

assert type(solve_term) == int and type(solve_half) == int
with warnings.catch_warnings():
	warnings.simplefilter("ignore")

	ttdata_half = ttdata_full[ttdata_full["term_half"]==solve_half].reset_index() # both terms
	ttdata_T5 = ttdata_half[ttdata_half["term"]=="ESD T5"].reset_index() # T5 only
	ttdata_T7 = ttdata_half[ttdata_half["term"]!="ESD T5"].reset_index() # T7 & TAE only

	ttdata_other_half = ttdata_otherpillars[ttdata_otherpillars["term_half"]==solve_half]

if solve_term == 0:
	df = ttdata_half
elif solve_term == 5:
	df = ttdata_T5
else:
	df = ttdata_T7

# multidict for ESD subjects
jobs, proc, instructor, venue, subject, prec, class_num, ft_core, ft_av, ft_ba, ft_fn, ft_sc, ft_ui = gp.multidict({
	df["job_ix"][i]: [
		df["proc_time"][i],
		df["instructor"][i],
		df["venue_ix"][i],
		df["subject_ix"][i],
		df["precedence"][i],
		df["class_num"][i],
		df["track_core"][i],
		df["track_avi"][i],
		df["track_ba"][i],
		df["track_fin"][i],
		df["track_scl"][i],
		df["track_uis"][i],
	] for i in range(len(df))
})

# multidict for non-ESD subjects
blocker, non_ESD_venue, block_start, block_end = gp.multidict({
	ttdata_other_half["blocker"][i]: [
		ttdata_other_half["venue_ix"][i],
		ttdata_other_half["release_date"][i],
		ttdata_other_half["deadline"][i]
	] for i in range(len(ttdata_other_half))
})

print("\nConverted csv data into multidict")

### Define variables & objective function ###################################
completion_time = {(j,v,t): t + proc[j] for j in jobs for v in set(venue.values()) for t in range(1,T+1)}

print("Problem space: %d ESD jobs x %d venues x %d time indices = %d\n"
% (len(jobs), len(set(venue.values())), T, len(jobs)*len(set(venue.values()))*T))

print("Creating Gurobi model:")
m = gp.Model("Adv Opti Task 2 by Kutosotase 2021")
X = m.addVars(completion_time.keys(), vtype=GRB.BINARY, name="X")

time0 = time.time() # to check solution time

# Variables to ensure classes of a subject are spread out across days
EOD_splits = [(j,k,t) for j in jobs for k in jobs for t in range(4)]
E = m.addVars(EOD_splits, vtype=GRB.BINARY, name="E")

# Variables to prevent clashes between focus track subjects
if solve_term  == 0:
	ft_list = [ft_core, ft_av, ft_ba, ft_fn, ft_sc, ft_ui]
elif solve_term == 5:
	ft_list = [ft_core]
else:
	ft_list = [ft_av, ft_ba, ft_fn, ft_sc, ft_ui]

# Variables to penalise the lack of breaks between classes of a focus track
Z = m.addVars([(j,k,t) for j in jobs for k in jobs for t in range(1,T+1)], vtype=GRB.BINARY, name="Z")

# Objective: minimise the lack of breaks
m.setObjective(gp.quicksum(Z[j,k,t] for j in jobs for k in jobs for t in range(1,T+1)), GRB.MINIMIZE)

print("\nObjective function defined")

### Define constraints ######################################################

### B1: All classes must be assigned
m.addConstrs((X.sum(j,"*","*")==1 for j in jobs), "all jobs assigned")

### B2: At most 1 class at any venue at any time instant
for v in set(venue.values()):
	for t in range(1,T+1):
		jobs_in_venue = 0
		for j in jobs:
			for s in range(max(1,t+1-proc[j]), t+1):
				jobs_in_venue += X[j,v,s]
		m.addConstr(jobs_in_venue <= 1, "<=1 job per venue per time")

### B3: At most 1 class per instructor at any time instant
for instr in set(instructor.values()):
	jobs_under_instructor = [j for j in jobs if instructor[j]==instr]
	#print("%s %s" % (instr, jobs_under_instructor))

	for t in range(1,T+1):
		jobs_being_done = 0
		for j in jobs_under_instructor:
			for s in range(max(1,t+1-proc[j]), t+1):
				jobs_being_done += X.sum(j,"*",s)
		m.addConstr(jobs_being_done <= 1, "<=1 job per prof per time")

### B5: Venues are blocked by other pillars/HASS
for v in set(venue.values()):
	for b in blocker:
		if v == non_ESD_venue[b]:
			m.addConstrs((X.sum("*",v,t)==0 for t in range(block_start[b], block_end[b]+1)), "venues blocked by other pillars")

### B6: Time slots are blocked by EOD and the blanket block period for HASS
for t in EOD_timings:
	m.addConstrs((X.sum(j,"*",s)==0 for j in jobs for s in range(max(1,t-proc[j]+1), t+1)), "time blocked by EOD")
for t in HASS_block_timings:
	m.addConstrs((X.sum(j,"*",s)==0 for j in jobs for s in range(max(1,t-proc[j]+1), t+1)), "time blocked for HASS")

### B4 & A2: Precedence constraints within the same subject & At most one session of each subject in a day
lecture_str = re.compile("L")
for j in jobs:
	for k in jobs:
		if j != k and subject[j] == subject[k] and prec[j] < prec[k]:
			# If same cohort (e.g. both CS01) or either class is a lecture
			if class_num[j]==class_num[k] or bool(lecture_str.search(class_num[j])) or bool(lecture_str.search(class_num[k])):
				# Then constrain j to precede k
				m.addConstr(sum(X.sum(j,"*",t)*t for t in range(1,T+1)) + 1 <= sum(X.sum(k,"*",t)*t for t in range(1,T+1)), "precedence constr")
			
			# j must be before some EOD block, and k after that same block
			for eod in range(4):
				m.addConstr(sum(X.sum(j,"*",t)*t for t in range(1,T+1)) <= EOD_timings[eod] + (1-E[j,k,eod])*T)
				m.addConstr(sum(X.sum(k,"*",t)*t for t in range(1,T+1)) >= EOD_timings[eod] - (1-E[j,k,eod])*T)

			# Exactly one EOD block must be "sandwiched" at a time
			m.addConstr(sum(E[j,k,eod] for eod in range(4)) == 1)

### 2.503/2.504 blocks 2.503 and 2.504, vice versa
for key, vals in venue_overlap.items():
	for val in vals:
		m.addConstrs((X.sum("*",key,t)+(X.sum("*",val,t))<=1 for t in range(1,T+1)), "venue overlaps")

### j=>v: Sessions must occupy the venue stipulated in the original timetable
for j in jobs:
	for v in set(venue.values()):
		m.addConstrs(X[j,v,t]==0 for t in range(1,T+1) if v != venue[j])

### A3 & A4: Sessions for subjects in the same focus track do not clash
cohort_str = re.compile("C")
for j in jobs:
	for k in jobs:
		for ft in ft_list:
			# If j and k are in the same track but are different subjects
			if j != k and ft[j] == 1 and ft[k] == 1 and subject[j] != subject[k]:
				# Prevent k from clashing with j in the following two cases

				# (1) Both j and k are in the same cohort (i.e. both CS01 or both CS02)
				if class_num[j] == class_num[k] and bool(cohort_str.search(class_num[j])):
					for t in range(1,T+1):
						clash = sum(X[j,venue[j],s] for s in range(max(1,t+1-proc[j]),t+1))
						clash += sum(X[k,venue[k],s] for s in range(max(1,t+1-proc[j]),t+1))
						m.addConstr(clash <= 1, "no focus track course clashes (case 1)")
				
				# Implied: j and k are different cohorts or one is a lecture
				else:
					num_cohort_subject_j = len(set(class_num[i] for i in jobs if subject[i]==subject[j]))
					num_cohort_subject_k = len(set(class_num[i] for i in jobs if subject[i]==subject[k]))

					# (2) j from single-cohort subject, k from multi-cohort subject of the same track
					if num_cohort_subject_j == 1 and num_cohort_subject_k > 1:
						for t in range(1,T+1):
							clash = sum(X[j,venue[j],s] for s in range(max(1,t+1-proc[j]),t+1))
							clash += sum(X[k,venue[k],s] for s in range(max(1,t+1-proc[j]),t+1))
							m.addConstr(clash <= 1, "no focus track course clashes (case 2)")

'''
### Illustration for A3 & A4 constraints

Step 1:
	j,k are in the same track, but are different subjects

Step 2:
	j        | k        | can clash? | remarks
	-------------------------------------------
	1 of 1   | 1 of 1   | N          |
	1 of 1   | 1 of 1-2 | N          | a FT student shld take all CS01s
	1 of 1-2 | 1 of 1-2 | N          |
	2 of 1-2 | 2 of 1-2 | N          | or all CS02s
	-------------------------------------------
	1 of 1   | 2 of 1-2 | Y          | let other students take 
	1 of 1-2 | 2 of 1-2 | Y          | mixtures of CS01/CS02

'''

### A1: Preferably at least a 30min break between sessions for subjects in the same track
for ft in ft_list:
	for j in jobs:
		for k in jobs:
			# If j and k are in the same track
			if j != k and ft[j] == 1 and ft[k] == 1:
				# Add a break between sessions in the following two cases

				# (1) If same cohort (e.g. both CS01) or either class is a lecture
				if class_num[j] == class_num[k] or bool(lecture_str.search(class_num[j])) or bool(lecture_str.search(class_num[k])):
					for t in range(1,T+1):
						classcount_for_break = 0
						for s in range(max(1,t-proc[j]),t+1):
							classcount_for_break += X[j,venue[j],s]
							classcount_for_break += X[k,venue[k],s]
						m.addConstr((classcount_for_break <= 1 + Z[j,k,t]), "30min breaks (case 1)")
				
				# Implied: j and k are different cohorts
				# (2) j from single-cohort subject, k from multi-cohort subject of the same track
				elif subject[j] != subject[k]:
					num_cohort_subject_j = len(set(class_num[i] for i in jobs if subject[i]==subject[j]))
					num_cohort_subject_k = len(set(class_num[i] for i in jobs if subject[i]==subject[k]))

					if num_cohort_subject_j == 1 and num_cohort_subject_k > 1:
						for t in range(1,T+1):
							classcount_for_break = 0
							for s in range(max(1,t-proc[j]),t+1):
								classcount_for_break += X[j,venue[j],s]
								classcount_for_break += X[k,venue[k],s]
							m.addConstr((classcount_for_break <= 1 + Z[j,k,t]), "30min breaks (case 2)")

print("Constraints defined, starting optimization\n")
time1 = time.time() # to check solution time

### Solve and print results #################################################
m.setParam("TimeLimit", 30) # time limit to search for a solution
m.optimize()
time2 = time.time() # to check solution time

for x in X.values():
	if (x.x > 0.5): print("%s %g" % (x.varName, x.x)) # print decision variables that equal 1
print("Obj: %g" % m.objVal)
print("Problem modelled in %g seconds" % (round(time1-time0)))
print("Optimal solution found in %g seconds" % (round(time2-time1)))
