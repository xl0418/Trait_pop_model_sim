# run_max.py
import subprocess

# Define command and arguments
command = 'ProcessPhylo'
path2script = 'C:/Liang/Trait_pop_model_sim/abcpp/ProcessPhyloTree.R'

# Variable number of args in a list
args = ['C:/Liang/abcpp_ms10/abcpp/treedata/bw_char.nex', 'C:/Liang/UI_test_treedata']

# Build subprocess command
cmd = [command, path2script] + args

# check_output will run the command and store to result
x = subprocess.Popen(cmd, universal_newlines=True,shell=False)

print('The maximum of the numbers is:', x)