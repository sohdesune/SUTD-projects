# -*- coding: utf-8 -*-
"""
Created on Sun Apr 08 10:37:18 2018

@author: sohdesune
"""
import random
import math

class PrimerDesign(object):
    
    def __init__ (self, name):
        
        '''parameters for the length criterion'''
        self.max_length = 22
        self.min_length = 18
        self.penalty_length = 10
        
        '''parameters for the temperature difference criterion'''
        self.max_tdiff = 1
        self.min_tdiff = 0
        self.penalty_tdiff = 10
        
        '''parameters for the cg content criterion'''
        self.max_cg = 0.6
        self.min_cg = 0.4
        self.penalty_cg = 10
        
        '''parameters for the annealing temperature criterion'''
        self.max_temp = 65
        self.min_temp = 55
        self.penalty_temp = 10
        
        '''parameters for the run criterion'''
        self.run_threshold = 4
        self.penalty_runs = 10
        
        '''parameters for the repeat criterion'''
        self.repeat_threshold = 0
        self.penalty_repeats = 10
        
        '''parameters for the specificity criterion'''
        self.penalty_specificity = 10 
        
        '''locations where the forward primer should be chosen from'''
        self.fp_start = 0
        self.fp_end = 980
        
        '''locations where the reverse primer should be chosen from'''
        self.rp_start = 990
        self.rp_end = 2090
        
        ''' parameters for the simulated annealing portion'''
        self.initial_temperature = 200
        self.stopping_temperature = 0.01
        self.drop_fraction = 0.999
    
    def set_dna_sequence(self, input_dna):
        bases = 'atcg'
        result = ''
        
        for char in input_dna:
            if char in bases:
                result += char
        
        self.dna_sequence = result
    
    def func_select_random(self, sqtype='forward', length=20):
        
        '''ensure length is positive'''
        if length < 1:
            length = 1
        
        '''choose limits based on primer type'''
        if(sqtype == 'forward'):
            start_limit = self.fp_start 
            end_limit = self.fp_end
        elif(sqtype == 'reverse'):
            start_limit = self.rp_start 
            end_limit = self.rp_end
        else:
            return None
        
        '''choose primer randomly'''
        first_base = random.randint(start_limit, end_limit-length+1)
        last_base = first_base+length-1
        
        result_primer = self.dna_sequence[first_base : last_base+1]
        
        return result_primer

    def func_length(self, sq):
        return len(sq)
    
    def func_cg_fraction(self, sq):
        cg_amt = 0
        
        '''count the number of Cs and Gs'''
        for base in sq:
            if base == 'c' or base == 'g':
                cg_amt += 1
        
        '''get ratio of Cs and Gs to entire sequence'''
        try:
            cg_fraction = cg_amt / len(sq)
        except ZeroDivisionError:
            cg_fraction = 0.
        
        return cg_fraction
    
    def func_temperature(self, sq):
        no_of_cg = 0
        no_of_at = 0
        
        '''in dna_sequence, only a, t, c, g exist'''
        '''any sq taken from dna_sequence will also only have these four characters'''
        for base in sq:
            if base == 'c' or base == 'g':
                no_of_cg += 1
            else:
                no_of_at += 1
        
        '''following equation as per Biology slides'''
        annealing_temperature = (4 * no_of_cg) + (2 * no_of_at)
        
        return annealing_temperature

    def func_count_runs(self,sq):
        numruns = 0
        before = ['at', 'ta', 'ac', 'ca', 'tc', 'ct', 'tg', 'gt', 'cg', 'gc', 'ga', 'ag']
        after = ['a t', 't a', 'a c', 'c a', 't c', 'c t', 't g', 'g t', 'c g', 'g c', 'g a', 'a g']
        
        '''split sq at points where characters change'''
        for i in range(len(before)):
            sq = sq.replace(before[i], after[i])
        
        '''split sq at points where a space has been inserted'''
        sq_list = sq.split()
        
        '''each element in sq_list is a run of the same character'''
        '''count the number of runs that exceed the run_threshold'''
        for i in sq_list:
            if len(i) > self.run_threshold:
                numruns += 1
        
        return numruns

    def func_count_repeats(self,sq):
        repeats = 0
        di_repeats = ['at','ac','ag','ca','ct','cg','ga','gt','gc','ta','tc','tg']
        
        '''find repeats of the two-base sequences in di_repeats'''
        for item in di_repeats:
            count = 0
            for i in range(len(sq)-4):
                '''example: in tgtgtg, the first and second tg are counted'''
                '''the third tg is not, so the number of repeats is 2'''
                if sq[i:i+2] == item and sq[i+2:i+4] == item:
                    count += 1
                    j = i+2
                    while sq[j:j+2] == item and sq[j+2:j+4] == item:
                        count += 1
                        j += 2
                    count = 0
                    repeats += 1
                else:
                    count = 0
        
        return repeats

    '''cost arising from primer length'''
    def cost_length(self, sq):
        sq_len = self.func_length(sq)
        
        '''compare primer length with limits'''
        if sq_len > self.max_length:
            return (sq_len - self.max_length) * self.penalty_length
        elif sq_len < self.min_length:
            return (self.min_length - sq_len) * self.penalty_length
        else:
            return 0
    
    '''cost arising from primer annealing temperature'''
    def cost_temperature(self, sq):
        sq_temp = self.func_temperature(sq)
        
        '''compare primer annealing temperature with limits'''
        if sq_temp > self.max_temp:
            return (sq_temp - self.max_temp) * self.penalty_temp
        elif sq_temp < self.min_temp:
            return (self.min_temp - sq_temp) * self.penalty_temp
        else:
            return 0
    
    '''cost arising from amount of C and G in primer'''
    def cost_cgcontent(self, sq):
        cg_content = self.func_cg_fraction(sq)
        
        '''compare amount of C and G in primer with limits'''
        if cg_content > self.max_cg:
            return (cg_content - self.max_cg) * self.penalty_cg
        elif cg_content < self.min_cg:
            return (self.min_cg - cg_content) * self.penalty_cg
        else:
            return 0
    
    '''cost arising from difference in annealing temperature of fp and rp'''
    def cost_temperature_difference(self, fp, rp):
        fp_temp = self.func_temperature(fp)
        rp_temp = self.func_temperature(rp)
        
        '''compute temperature difference'''
        t_diff = abs(fp_temp - rp_temp)
        
        '''compare temperature difference with limit'''
        if t_diff > self.max_tdiff:
            return (t_diff - self.max_tdiff) * self.penalty_tdiff
        else:
            return 0
    
    '''cost arising from uniqueness of primer within dna_sequence'''
    def cost_specificity(self, sq):
        '''count number of occurrences of primer'''
        n_positions = self.dna_sequence.count(sq)
        
        '''exclude failure because sq not in dna_sequence'''
        if n_positions == 0:
            n_positions = 1
        
        '''compute cost'''
        result = (n_positions - 1) * self.penalty_specificity
        
        return result
    
    '''cost arising from unacceptably long runs in primer'''
    def cost_runs(self, sq):
        numruns = self.func_count_runs(sq)
        
        return numruns * self.penalty_runs
    
    '''cost arising from repeats of two-base sequences in primer'''
    def cost_repeats(self, sq):
        reps = self.func_count_repeats(sq)
        
        return reps * self.penalty_repeats

    def calculate_values(self, fp, rp):
        values_dict = {}
        
        values_dict['Length (fp)'] = self.func_length(fp)
        values_dict['Length (rp)'] = self.func_length(rp)
        values_dict['% CG content (fp)'] = self.func_cg_fraction(fp)
        values_dict['% CG content (rp)'] = self.func_cg_fraction(rp)
        values_dict['Annealing temp. (fp)'] = self.func_temperature(fp)
        values_dict['Annealing temp. (rp)'] = self.func_temperature(rp)
        values_dict['Temp. difference'] = abs(values_dict['Annealing temp. (fp)']-values_dict['Annealing temp. (rp)'])
        values_dict['Specificity (fp)'] = self.dna_sequence.count(fp)
        values_dict['Specificity (rp)'] = self.dna_sequence.count(fp)
        values_dict['Runs (fp)'] = self.func_count_runs(fp)
        values_dict['Runs (rp)'] = self.func_count_runs(rp)
        values_dict['Repeats (fp)'] = self.func_count_repeats(fp)
        values_dict['Repeats (rp)'] = self.func_count_repeats(rp)
        
        return values_dict
    
    def calculate_costs(self, fp, rp):
        cost_dict = {}
        
        cost_dict['Length (fp)'] = self.cost_length(fp)
        cost_dict['Length (rp)'] = self.cost_length(rp)
        cost_dict['% CG content (fp)'] = self.cost_cgcontent(fp)
        cost_dict['% CG content (rp)'] = self.cost_cgcontent(rp)
        cost_dict['Annealing temp. (fp)'] = self.cost_temperature(fp)
        cost_dict['Annealing temp. (rp)'] = self.cost_temperature(rp)
        cost_dict['Temp. difference'] = self.cost_temperature_difference(fp, rp)
        cost_dict['Specificity (fp)'] = self.cost_specificity(fp)
        cost_dict['Specificity (rp)'] = self.cost_specificity(rp)
        cost_dict['Runs (fp)'] = self.cost_runs(fp)
        cost_dict['Runs (rp)'] = self.cost_runs(rp)
        cost_dict['Repeats (fp)'] = self.cost_repeats(fp)
        cost_dict['Repeats (rp)'] = self.cost_repeats(rp)
        
        return cost_dict
    
    def cost_objective_function(self, fp, rp):
        '''build dictionary of respective costs'''
        cost_dict = self.calculate_costs(fp, rp)
        
        '''objective function value is the sum of all of the above costs'''
        score = sum(cost_dict.values())
        
        return score

    def cost_objective_function_info(self, fp, rp):
        print("===============================================")
        print("{:^47}\n{:^47}".format("Bio-DW 2D", "Primer Cost Report"))
        print("===============================================\n")
        print("Forward Primer: {}\nReverse Primer: {}".format(fp, rp))
        print("\nNote: The rp sequence above is the reverse-\ncomplement, i.e. on same DNA strand as fp.")
        
        values_dict = self.calculate_values(fp, rp)
        cost_dict = self.calculate_costs(fp, rp)
        print("\n=== Cost Function Score ===")
        score = sum(cost_dict.values())
        print(round(score,2))
        
        '''0 is not used in case of floating-point error'''
        if score < 1e-3:
            print("Excellent! All criteria met!")
        
        else:
            print("\n=== Unmet Criteria ===")
            print("{:<25}{:>11}{:>11}\n-----------------------------------------------".format("Criterion", "Value", "Cost"))
            for key in cost_dict:
                if cost_dict[key] > 1e-3:
                    print('{:<25}{:>11.2f}{:>11.2f}'.format(key, values_dict[key], cost_dict[key]))

        print("\n=== Met Criteria ===")
        print("{:<25}{:>11}{:>11}\n-----------------------------------------------".format("Criterion","Value", "Cost"))
        for key in cost_dict:
            if cost_dict[key] < 1e-3:
                print('{:<25}{:>11.2f}{:>11.2f}'.format(key, values_dict[key], cost_dict[key]))

    def get_primer_neighbour(self, current_primer):
        '''obtain position of current primer'''
        current_index = self.dna_sequence.find(current_primer)
        
        '''generate new primer randomly'''
        '''choose length within desirable range, no point exceeding'''
        new_length = random.randint(self.min_length, self.max_length)
        first_base = random.randint(current_index-10, current_index+10)
        
        '''prevent choosing bad primer due to being near start or end of dna_sequence'''
        while first_base < 0 or first_base > len(self.dna_sequence)-new_length:
            first_base = random.randint(current_index-10, current_index+10)
        
        '''select neighbouring primer'''
        last_base = first_base+new_length-1
        result_primer = self.dna_sequence[first_base : last_base+1]
        
        return result_primer
    
    def func_simulated_annealing(self):
        '''step 1.1 - decide simulation parameters'''
        temperature = self.initial_temperature
        stopping_temperature = self.stopping_temperature
        drop = self.drop_fraction
        
        print("=== Simulation parameters ===")
        print("Initial temperature: {}\nStopping temperature: {}\nDrop Fraction: {}"
                  .format(temperature, stopping_temperature, drop))
        
        '''step 1.2 - start by choosing possible solution'''
        fp_current = self.func_select_random(sqtype='forward', length = 20)
        rp_current = self.func_select_random(sqtype='reverse', length = 20)
        cost = self.cost_objective_function(fp_current, rp_current)
        current_cost = cost
        
        print("Initial FP (random): {}\nInitial RP (random): {}\nCost: {}\n"
                  .format(fp_current, rp_current, current_cost))
        
        i = 0 #to count iterations completed
        
        while temperature > stopping_temperature:
            i += 1
            '''step 2 - obtain new possible solution'''
            fp_new = self.get_primer_neighbour(fp_current)
            rp_new = self.get_primer_neighbour(rp_current)
            print("Running iteration {:>4}: fp ({:<22}), rp ({:<22})"
                      .format(i, fp_new, rp_new))
            new_cost = self.cost_objective_function(fp_new, rp_new)
            
            '''step 3 - if new solution has lower cost, accept immediately'''
            if new_cost < current_cost:
                fp_current = fp_new
                rp_current = rp_new
                current_cost = new_cost
            
            else:
                '''step 4a - calculate an acceptance probability'''
                delta = new_cost - current_cost
                acceptance_probability = math.exp(-delta / temperature)
                
                '''step 4b - roll die'''
                num = random.random()
                
                '''step 4c - check if accept new solution'''
                if acceptance_probability > num:
                    fp_current = fp_new
                    rp_current = rp_new
                    current_cost = new_cost
            
            '''decrease temperature by a small factor'''
            temperature *= drop
        
        '''after arriving at solution, display results'''
        print("\n=== Final solution ===")
        print("Forward primer: {}\nReverse primer: {}\nCost: {}"
                  .format(fp_current, rp_current, current_cost))
        print("No. of iterations: {}".format(i))
        print("\nDetails:")
        self.cost_objective_function_info(fp_current, rp_current)

'''
===============================================================================
===============================================================================

                              USING THE CODE

===============================================================================
===============================================================================
'''
# key in DNA sequence here; change fp_start/end and rp_start/end accordingly
dna_sequence = ''
my_primer = PrimerDesign("my_primer")
my_primer.set_dna_sequence(dna_sequence)

answer = input("(A) calculate cost or (B) simulate annealing?\nChoose A/B: ")

if answer == "A": # display cost objective function results
    fwd = input("Type in forward primer: ")
    rev = input("Type in reverse primer: ")
    my_primer.cost_objective_function_info(fwd, rev)

elif answer == "B": # perform simulated annealing
    my_primer.func_simulated_annealing()

else:
    print("Invalid input.")