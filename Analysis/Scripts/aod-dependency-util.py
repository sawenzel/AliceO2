#!/usr/bin/env python3

# A small util to:
# a) determine analysis executables (path) required in order to run a given analysis x
# b) determine the set of fundamental AOD tables (from data model) needed by analysis x
#
# The tool basically takes the universe of analysis executables and produces
# an easily parseable json file that can be used for further processing (graph drawing or querying with Pandas dataframes).
# The file can easily be extended with further information.
#
# An example invocation is the following:
# aod-dependency-util.py --build-from-universe `ls ${O2_ROOT}/bin/o2-analysis*` -o out.json
#
# The file out.json will contain entries such as:
#  "o2-analysis-spectra-tof-tiny": {"path": ["o2-analysis-alice3-trackselection", "o2-analysis-pid-tof", "o2-analysis-alice3-trackextension"],
#                                   "needs-data-model-AOD": ["AOD/COLLISION", "AOD/TRACKEXTRA", "AOD/TRACK", "AOD/TRACKCOV"], ..}
#
#  listing one possible set of required analysis, as well as the fundamental AOD tables used in this analysis.

import argparse
import json
import subprocess

parser = argparse.ArgumentParser(description='ALICE O2 analysis tooling',
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('--build-from-universe', nargs="+", help='Build from a list of analysis tasks, defining the universe')
parser.add_argument('--universe-json-out', help='Global JSON dump of analysis universe')
parser.add_argument('--build-from-json', help='Build from a pregenerated global JSON (result of build-from-universe)')
parser.add_argument('-o', help='Final json output file (containing calculated properties)', default='o2-analysis-dependencies.json')
args = parser.parse_args()

# some helper structures which get filled
exec_to_workflowlist={}
exec_to_tasks={}
task_to_exec={} # can be computed from first one
task_to_inputs={}
task_to_outputs={}
input_to_id={}
output_to_id={}
task_to_id={}

# a list of fundamental AOD tables (incomplete)
fundamental=["AOD/CALO", # DM
             # "AOD/CALOTRIGGER", ?
             "AOD/V0", # DM
             "AOD/COLLISION", # DM
             "AOD/BC", # DM
             "AOD/RUN2BCINFO", # DM
             "AOD/TRACK", # DM
             "AOD/TRACKEXTRA", # DM
             "AOD/TRACKCOV", # DM
             "AOD/HMPID", # DM
             "AOD/MCTRACKLABEL", # DM
             "AOD/MCCOLLISION", # DM
             "AOD/MCPARTICLE", # DM
             "AOD/MCCOLLISLABEL", # DM
             "AOD/MCCALOLABEL", # DM
             "AOD/FT0",  # DM
             "AOD/FDD", # DM
             "AOD/FV0A", # DM
             "AOD/CASCADE", #DM
             "AOD/FV0C", #DM
             "AOD/ZDC", #DM
             "AOD/FWDTRACK",   #DM
             "AOD/FWDTRACKCOV"  #DM
]
# might be missing: AmbigTracks, AmbigMFTTracks, MUON, MFTTracks, MCMFTTrackLabel

# launches an analysis and extracts workflow json
def getjson(analysis):
    output = subprocess.check_output([analysis, '--dump-workflow'])
    spec = json.loads(output)
    return spec['workflow']

def add_namedtask(task, executable):
   if exec_to_tasks.get(executable)==None:
      exec_to_tasks[executable]=[]
   exec_to_tasks[executable].append(task)
   task_to_exec[task]=executable
   for inputs in t['inputs']:
      if task_to_inputs.get(task)==None:
         task_to_inputs[task]=[]
      task_to_inputs[task].append(inputs['origin'] + '/' + inputs['description'])
      for outputs in t['outputs']:
        if task_to_outputs.get(task)==None:
           task_to_outputs[task]=[]
        task_to_outputs[task].append(outputs['origin'] + '/' + outputs['description'])

# def getlistofp
if args.build_from_universe:  #<---- this is build from a given list of executables
  globalspec = []
  for a in args.build_from_universe:
     js = getjson(a)
     if args.universe_json_out:
       globalspec.append({'executable' : a, 'spec' : js})
     for t in js:
       task=t['name']
       add_namedtask(task, a)
     if args.universe_json_out:
       with open(args.universe_json_out, 'w') as outfile:
         json.dump(globalspec, outfile)

if args.build_from_json:  #<------ this is building from a global json (generated with universe_json_out)
    with open(args.build_from_json, 'r') as fp:
      globaljson=json.load(fp)
      for entry in globaljson:
        js = entry['spec']
        for t in js:
          add_namedtask(t['name'], entry['executable'])

all_inputs=set([ x for _,l in task_to_inputs.items() for x in l ])
all_outputs=set([ x for _,l in task_to_outputs.items() for x in l ])

# try to find out if input is fundamental
provided_by = {}  # keeps track of which input is provided by which task (is a list)
for k, outputs in task_to_outputs.items():
  for o in outputs:
    if provided_by.get(o) == None:
       provided_by[o] = []
    provided_by[o].append(k)

# calculated one possible exec dependency path (multiple may exist)
def get_one_exec_dependency_path(a):
    result = []
    inputlist = task_to_inputs[a]
    for i in inputlist:
      pr = provided_by.get(i)
      if pr != None:
        t = task_to_exec[pr[0]]
        result.append(t)
        result = result + get_one_exec_dependency_path(pr[0]) # the plus may be a times

    return result

# calculates the set of fundamental data model dependencies for a given analysis
def get_all_fundamental_dependencies(a):
    result = []
    inputlist = task_to_inputs[a]
    for i in inputlist:
       if fundamental.count(i) > 0:
          result.append(i)
       else:
          pr = provided_by.get(i)
          if pr!=None:
            result = result + get_all_fundamental_dependencies(pr[0])
    return result

# now detect if fundamental or list all dependencies of this analysis
depmap={} # maps executable to possible path and fundamental AODs needed
for a in exec_to_tasks:
   t = exec_to_tasks[a]
   execpath = list(set(get_one_exec_dependency_path(t[0])))
   deps = []
   for t in exec_to_tasks[a]:
     deps = deps + get_all_fundamental_dependencies(t)
   depmap[a]={'path' : execpath, 'needs-data-model-AOD' : list(set(deps))}

   # enrich this by a yes/false annotation for each fundamental AOD
   # which may make is easier to query in pandas
   for dm in fundamental:
     depmap[a][dm]=depmap[a]['needs-data-model-AOD'].count(dm) > 0

# export the result as a json
with open(args.o, 'w') as outfile:
  json.dump(depmap,outfile)
