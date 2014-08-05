#!/usr/bin/env python

# DI Johannes Kasberger
# (C) 2014 Reelworx GmbH, http://reelworx.at
# https://github.com/schugabe/isotope-mass-calculator

from __future__ import print_function
import sys, getopt
import csv
import subprocess
import os
from math import pow
from decimal import *
from collections import defaultdict

def usage(name):
    print("Enumerates all sums of isotope masses that yield a certain molecule mass")
    print("The isotopes used are loaded from a file (default input.csv).")
    print("The format of each line in the input file is name;mass;abundance;")
    print("For each found solution a relevance value is computed")
    print("If two isotopes have the same mass one must be excluded since the algorithm currently doesn't support this case")
    print()
    print("USAGE: "+name+" options")
    print("OPTIONS:")
    print(" -t Specifies the targets, Multiple Targets can be specified with -t 194,135,136")
    print(" -i Specifies the input file used e.g. -i data.csv")
    print(" -o The output can be stored as csv files. For each target one csv file is generated")
    print("    e.g. -o output.csv generates output_194.csv")
    print(" -u Specify isotopes used with -u name1,name2,... isotopes not in the list are not used")
    print(" -e Exclude isotopes with -e name1,name2,...")
    print()
    print("EXAMPLES")
    print(" "+name+" -t 194")
    print("  Loads input.csv in the same directory, uses all isotopes in the file to find combinations")
    print("  that sum up to 194 and prints the result")
    print()
    print(" "+name+" -t 194 -o output.csv")
    print("  Loads input.csv in the same directory, uses all isotopes in the file to find combinations")
    print("  that sum up to 194 and stores the result in output_194.csv")
    print()
    print(" "+name+" -t 194 -u 138La,137Ba")
    print("  Loads input.csv in the same directory, uses only 138La and 137Ba from the file to find combinations")
    print("  that sum up to 194 and stores the result in output_194.csv")
    print("  to load 63Cu and 65Cu use the argument -u Cu, this works for all isotopes that share the same suffix")
    print()
    print(" "+name+" -t 194 -e 138La,137Ba")
    print("  Loads input.csv in the same directory, uses all isotopes from the file except 138La and 137Ba")
    print("  to find combinations that sum up to 194 and stores the result in output_194.csv")
    print()
    
def compute_probability(current_combinations):
    probability = Decimal('1.0')
    for (n,c) in current_combinations:
        if n > 0:
            tmp_proability = propabilities[c] ** int(n)
            probability = probability * tmp_proability
    return probability
       
def build_combs(remaining_sum, i, current_combinations, add_combinations):
    
    if add_combinations: 
        current_combinations.append(add_combinations)
        
    if remaining_sum == 0 or i == len(masses):
        if i == len(masses) and remaining_sum > 0:
            return 0
        while i < len(masses):
            current_combinations.append( (0, masses[i]) )
            i += 1
        
        probability = compute_probability(current_combinations).quantize(Decimal('1.00'))
        result = ";".join(int_format % (n) for (n,c) in current_combinations)+";"+str(probability).rjust(9)+";"
        if results.get(probability):
            old_result = results[probability]
            result = old_result+"\n"+result
        results.update({probability: result})
        return 1
        
    return sum(build_combs(remaining_sum-x*masses[i], i+1, current_combinations[:], (x,masses[i])) for x in range(0, int(remaining_sum/masses[i])+1))

def main(argv):
    ifile = "input.csv"
    ofile = "output.csv"
    
    global masses
    global propabilities
    global results
    
    verbose = 0
    targets = []
    masses = []
    names = []
    propabilities = {}
    results = {}
    use_names = []
    exclude_names = []
    getcontext().prec = 28
    use_stdout = 1
    try:
        opts, args = getopt.getopt(argv[1:],"hi:o:t:vu:e:",["ifile=","ofile=","target=","use=","exclude="])
    except getopt.GetoptError:
        usage(argv[0])
        sys.exit(2)
    for opt, arg in opts:
        if opt in ("-h", "--help"):
            usage(argv[0])
            sys.exit()
        elif opt in ("-i", "--ifile"):
            ifile = arg
        elif opt in ("-o", "--ofile"):
            ofile = arg
            use_stdout = 0
        elif opt in ("-t", "--target"):
            targets = [int(x) for x in arg.split(",")]
        elif opt == "-v":
            verbose = 1
        elif opt in ("-u", "--use"):
            use_names = arg.split(",")
        elif opt in ("-e", "--exclude"):
            exclude_names = arg.split(",")
            
    if not use_stdout:
        if not os.path.isfile(ifile):
            print("Input file " + ifile + " not found")
            sys.exit(2)
    
    with open(ifile, 'rU') as csvfile:
        inputreader = csv.reader(csvfile, delimiter=';', quotechar='|')
        for row in inputreader:
            name = row[0]
            
            use_current = 0
            if len(use_names)==0:
                use_current = 1
            else:
                for use in use_names:
                    if name.endswith(use):
                        use_current = 1
                        
            if use_current and name not in exclude_names:
                names.append(name)
                mass = int(row[1])
                masses.append(mass)
                propabilities.update({int(row[1]): Decimal(row[2])})
                if verbose:
                    print("Adding " + names[-1] + " with mass " + str(masses[-1]) + " with relevance " + str(propabilities[masses[-1]]))
            else:
                if verbose:
                    print("Ignoring "+str(name))
    
    if len(masses) == 0:
        print("No isotopes added, nothing to do")
        sys.exit(0)
        
    duplicates  = 0
    
    max_length = len(max(names, key=len))
    
    global int_format
    int_format = "%"+str(max_length)+"d"
    
    duplicates_str = "-e "
    
    D = defaultdict(list)
    for i,item in enumerate(masses):
        D[item].append(i)
    for k,v in D.items():
        if len(v)>1:
            print("Error: duplicate mass %d, please exclude all but one of" % (k), end=" ")
            print(", ".join("%s" % (names[index]) for index in v))
            duplicates_str += ",".join("%s" % (names[index]) for index in v[:-1])+","
            duplicates = 1

    if duplicates:
        print(duplicates_str[:-1])
        sys.exit(0)
    
    for target in targets:
        print("Searching "+str(target), end=" ")
        solutions = build_combs(target, 0, [], None)
        print("and found %d solution(s)" % solutions)
        
        if solutions==0:
            continue
        
        tmp_filename = ofile.split('.')
        filename = ("%s_%d.%s" % (tmp_filename[0], target, tmp_filename[1]))
        
        header = ";".join("%s" % (name.rjust(max_length)) for name in names)+";Relevance;\n"
        text = "\n".join("%s" % (results[result]) for result in sorted(results, reverse=True))
        if use_stdout:
            print(header, end=" ")
            print(text)
        else:
            if verbose:
                print("Writing Result to "+filename)
            with open(filename, "w") as result_file:
                result_file.write(header)
                result_file.write(text)

if __name__ == '__main__':
    main(sys.argv)