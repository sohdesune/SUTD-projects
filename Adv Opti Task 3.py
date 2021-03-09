# Specify optimisation variant (1 or 2 or 3, or 4 [all three])
solve_variant = 4

# Specify term (5 or 7, or 0 [both]) and half (1 or 2) to solve for
solve_term = 0
solve_half = 1

###################################################################################################

T = 92 # the last time-index, corresponding to the Fri 1300-1330 timeslot
EOD_timings = [24, 48, 58, 82, 92] # last time-index of each day
HASS_block_timings = [14,15,16,17,18,19,25,26,27,72,73,74,75,76,77,89,90,91] # time-indices blocked for HASS
venue_overlap = {6:[4,5], 32:[30,31]} # e.g. 1.415/1.416: [1.415, 1.416]

# weights of time-indices for Variant 3
time_ix_wt = [
	2,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,
	2,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,
	2,2,0,0,0,0,0,0,0,0,
	2,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,
	2,2,0,0,0,0,0,0,0,0
	]

# relative weights of objectives for Variant 4 (multi-objective optimisation)
obj_wt = [1, 1, 1] # [Var 1 obj, Var 2 obj, Var 3 obj]

max_solve_time = 30 # time limit for Gurobi to try to find a solution (seconds)

###################################################################################################

import gurobipy as gp
from gurobipy import GRB
import pandas as pd
import os
import warnings
import re
import time

### Read Data & Create Multidict ##################################################################

os.chdir("set your wd here")
ttdata_full = pd.read_csv('TT Data v6.csv')
ttdata_otherpillars = pd.read_csv('TT_NonESD_Data v2.csv')

assert solve_variant == 1 or solve_variant == 2 or solve_variant == 3 or solve_variant == 4, "Invalid solve_variant specified"
assert solve_term == 0 or solve_term == 5 or solve_term == 7, "Invalid solve_term specified"
assert solve_half == 1 or solve_half == 2, "Invalid solve_half specified"

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

# Multidict for ESD subjects
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

# Multidict for non-ESD subjects
blocker, non_ESD_venue, block_start, block_end = gp.multidict({
	ttdata_other_half["blocker"][i]: [
		ttdata_other_half["venue_ix"][i],
		ttdata_other_half["release_date"][i],
		ttdata_other_half["deadline"][i]
	] for i in range(len(ttdata_other_half))
})

print("\nConverted csv data into multidict")


### Common Variables ##############################################################################

print("Problem space: %d ESD jobs x %d venues x %d time indices = %d\n"
% (len(jobs), len(set(venue.values())), T, len(jobs)*len(set(venue.values()))*T))

print("Creating Gurobi model:")
m = gp.Model("Adv Opti Task 3 by Kutosotase 2021")
X = m.addVars([(j,v,t) for j in jobs for v in set(venue.values()) for t in range(1,T+1)], vtype=GRB.BINARY, name="X")

time0 = time.time() # to check solution time

# Variables to ensure classes of a subject are spread out across days
EOD_splits = [(j,k,t) for j in jobs for k in jobs for t in range(4)]
E = m.addVars(EOD_splits, vtype=GRB.BINARY, name="E")

# Variables to prevent clashes between focus track subjects
if solve_term  == 0:
	ft_list = [ft_core, ft_av, ft_ba, ft_fn, ft_sc, ft_ui]
	ft_names = ["Term 5", "Aviation", "BA/OR", "Financial Services", "Supply Chain & Logistics", "Urban Infrastructure Systems"]
elif solve_term == 5:
	ft_list = [ft_core]
	ft_names = ["Term 5"]
else:
	ft_list = [ft_av, ft_ba, ft_fn, ft_sc, ft_ui]
	ft_names = ["Aviation", "BA/OR", "Financial Services", "Supply Chain & Logistics", "Urban Infrastructure Systems"]


### Additional Variables ##########################################################################

# Variant 1 - Provide 30min breaks between classes of a focus track
if solve_variant == 1 or solve_variant == 4:

	# Include HASS blockout timings in 30min break provision
	HASS_block_timings += [13, 20, 24, 28, 71, 78, 88]

	# Z = number of absences of 30min breaks
	Z = m.addVars([(j,k,t) for j in jobs for k in jobs for t in range(1,T+1)], vtype=GRB.BINARY, name="Z")

# Variant 2 - Compress each focus track into as few days as possible
if solve_variant == 2 or solve_variant == 4:

	# Number of cohort numbers (e.g. CS01, CS02)
	cohort_str = re.compile("C")
	num_cohort_total = len(set(class_num[j] for j in jobs if bool(cohort_str.search(class_num[j]))))

	# D = number of days taken by a focus track, excluding HASS
	D = m.addVars([(f,l,d) for f in range(len(ft_list)) for l in range(1, num_cohort_total+1) for d in range(1,5+1)], vtype=GRB.BINARY, name="D")
	Dmax = m.addVar(vtype=GRB.INTEGER, name="Dmax")

# Variant 3 - Prefer certain time slots (punish classes at undesirable timings)
if solve_variant == 3 or solve_variant == 4:

	pass # no additional variables to define


### Objective Function ############################################################################

if solve_variant == 1: # Variant 1 - Minimise the lack of breaks
	m.setObjective(gp.quicksum(Z[j,k,t] for j in jobs for k in jobs for t in range(1,T+1)), GRB.MINIMIZE)

elif solve_variant == 2: # Variant 2 - Minimise largest number of days taken per focus track
	m.setObjective(Dmax, GRB.MINIMIZE)

elif solve_variant == 3: # Variant 3 - Minimise sum of weighted start times
	m.setObjective(gp.quicksum(X.sum(j,"*",t)*time_ix_wt[t-1] for j in jobs for t in range(1,T+1)), GRB.MINIMIZE)

else: # Variant 4 - Minimise all three objectives
	obj1 = gp.quicksum(Z[j,k,t] for j in jobs for k in jobs for t in range(1,T+1))
	obj2 = Dmax
	obj3 = gp.quicksum(X.sum(j,"*",t)*time_ix_wt[t-1] for j in jobs for t in range(1,T+1))

	m.setObjective(obj1*obj_wt[0] + obj2*obj_wt[1] + obj3*obj_wt[2], GRB.MINIMIZE)

print("\nObjective function defined for Variant %g" % solve_variant)


### Common Constraints ############################################################################

### j=>v: Sessions must occupy the venue stipulated in the original timetable
for j in jobs:
	for v in set(venue.values()):
		m.addConstrs(X[j,v,t]==0 for t in range(1,T+1) if v != venue[j])

### All classes must be assigned
m.addConstrs((X.sum(j,"*","*")==1 for j in jobs), "all jobs assigned")

### At most 1 class at any venue at any time instant
for v in set(venue.values()):
	for t in range(1,T+1):
		jobs_in_venue = 0
		for j in jobs:
			for s in range(max(1,t+1-proc[j]), t+1):
				jobs_in_venue += X[j,v,s]
		m.addConstr(jobs_in_venue <= 1, "<=1 job per venue per time")

### At most 1 class per instructor at any time instant
for instr in set(instructor.values()):
	jobs_under_instructor = [j for j in jobs if instructor[j]==instr]

	for t in range(1,T+1):
		jobs_being_done = 0
		for j in jobs_under_instructor:
			for s in range(max(1,t+1-proc[j]), t+1):
				jobs_being_done += X.sum(j,"*",s)
		m.addConstr(jobs_being_done <= 1, "<=1 job per prof per time")

### Venues are blocked by other pillars/HASS
for v in set(venue.values()):
	for b in blocker:
		if v == non_ESD_venue[b]:
			m.addConstrs((X.sum("*",v,t)==0 for t in range(block_start[b], block_end[b]+1)), "venues blocked by other pillars")

### Time slots are blocked by EOD and the blanket block period for HASS
for t in EOD_timings:
	m.addConstrs((X.sum(j,"*",s)==0 for j in jobs for s in range(max(1,t-proc[j]+1), t+1)), "time blocked by EOD")
for t in HASS_block_timings:
	m.addConstrs((X.sum(j,"*",s)==0 for j in jobs for s in range(max(1,t-proc[j]+1), t+1)), "time blocked for HASS")

### Precedence constraints within the same subject
### & At most one session of each subject in a day
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

### Sessions for subjects in the same focus track do not clash
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
					num_cohort_subject_j = len(set(class_num[i] for i in jobs if subject[i]==subject[j] and not bool(lecture_str.search(class_num[i]))))
					num_cohort_subject_k = len(set(class_num[i] for i in jobs if subject[i]==subject[k] and not bool(lecture_str.search(class_num[i]))))

					# (2) j from single-cohort subject, k from multi-cohort subject of the same track
					if num_cohort_subject_j == 1 and num_cohort_subject_k > 1:
						# Since j & k are asymmetric, they are only compared once (jk only, no kj)
						proc_longer = max(proc[j], proc[k])
						for t in range(1,T+1):
							clash = sum(X[j,venue[j],s] for s in range(max(1,t+1-proc_longer),t+1))
							clash += sum(X[k,venue[k],s] for s in range(max(1,t+1-proc_longer),t+1))
							m.addConstr(clash <= 1, "no focus track course clashes (case 2)")


### Variant-Specific Constraints ##################################################################

### V1: Preferably at least a 30min break between sessions for subjects in the same track
if solve_variant == 1 or solve_variant == 4:
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
						num_cohort_subject_j = len(set(class_num[i] for i in jobs if subject[i]==subject[j] and not bool(lecture_str.search(class_num[i]))))
						num_cohort_subject_k = len(set(class_num[i] for i in jobs if subject[i]==subject[k] and not bool(lecture_str.search(class_num[i]))))

						if num_cohort_subject_j == 1 and num_cohort_subject_k > 1:
							# Since j & k are asymmetric, they are only compared once (jk only, no kj)
							proc_longer = max(proc[j], proc[k])
							for t in range(1,T+1):
								classcount_for_break = 0
								for s in range(max(1,t-proc_longer),t+1):
									classcount_for_break += X[j,venue[j],s]
									classcount_for_break += X[k,venue[k],s]
								m.addConstr((classcount_for_break <= 1 + Z[j,k,t]), "30min breaks (case 2)")

### V2: Dmax implementation
if solve_variant == 2 or solve_variant == 4:
	for ft in ft_list:
		f = ft_list.index(ft) # get index of focus track for D[f,l,d]

		# Check the number of cohorts in focus track f
		num_cohort_track = len(set(class_num[j] for j in jobs if ft[j]==1 and not bool(lecture_str.search(class_num[j]))))
		for l in range(1,num_cohort_track+1):

			for d in range(len(EOD_timings)):
				day_end = EOD_timings[d]
				day_start = 1 if d==0 else EOD_timings[d-1]+1

				ft_cohort_classes_in_day = 0
				for j in jobs:

					# If job j is part of the focus track
					if ft[j] == 1:

						assert bool(re.compile("[0-9][0-9]").search(class_num[j][-2:])), "Last 2 chars of class_num are not numerics"

						# First constraint: check if day d is used by cohort l of focus track f
						# If j is a lecture, always check
						if bool(lecture_str.search(class_num[j])):
							ft_cohort_classes_in_day += sum(X[j,venue[j],t] for t in range(day_start, day_end+1))
							
						# Elif j is the correct cohort number, then check
						elif int(class_num[j][-2:]) == l:
							ft_cohort_classes_in_day += sum(X[j,venue[j],t] for t in range(day_start, day_end+1))

						# Second constraint: if j is from a single-cohort subject, count the day for all cohorts
						num_cohort_subject_j = len(set(class_num[i] for i in jobs if subject[i]==subject[j] and not bool(lecture_str.search(class_num[j]))))
						if num_cohort_subject_j == 1:
							
							for l_iter in range(1,num_cohort_track+1):
								# Use T as big M
								m.addConstr((sum(X[j,venue[j],t] for t in range(day_start, day_end+1)) <= D[f, l_iter, d+1]*T), "all cohorts of FT use that day")
					
				# First constraint implementation (using T as big M)
				m.addConstr((ft_cohort_classes_in_day <= D[f,l,d+1]*T), "classes of a focus track in a day")

			# Implementation of Dmax objective
			m.addConstr((Dmax >= D.sum(f,l,"*")), "Dmax is largest num of days per cohort per focus track")

### V3
if solve_variant == 3 or solve_variant == 4:
	pass # implemented in objective function, no additional constraints required

print("Constraints defined, starting optimisation\n")
time1 = time.time() # to check solution time


### Solve & Print Results #########################################################################

m.setParam("TimeLimit", max_solve_time) # time limit to search for a solution
m.optimize()
time2 = time.time() # to check solution time

for x in X.values():
	if (x.x > 0.5):
		print("%s %g" % (x.varName, x.x)) # print decision variables that equal 1

print("\nObjective value: %g" % m.objVal)

if solve_variant == 1 or solve_variant == 4:
	print("\nNumber of back-to-back classes: %g" % sum(z.x for z in Z.values()))

if solve_variant == 2 or solve_variant == 4:
	print("\nDays per cohort per focus track:")
	def get_Dvals(d, f, l):
		dprint = d.varName.strip("D[]").split(",")
		if int(dprint[0]) == f and int(dprint[1]) == l:
			return True
		else:
			return False

	for ft in ft_list:
		f = ft_list.index(ft)
		num_cohort_track = len(set(class_num[j] for j in jobs if ft[j]==1 and not bool(lecture_str.search(class_num[j]))))
		for l in range(1,num_cohort_track+1):
			d_x = 0
			for d in D.values():
				if get_Dvals(d, f, l):
					d_x += d.x
			print("  %s Cohort %s: %g days" % (ft_names[f], l, d_x))

print("\nProblem modelled in %g seconds" % (round(time1-time0, 2)))
print("Optimal solution found in %g seconds" % (round(time2-time1, 2)))
