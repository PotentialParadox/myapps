'''
The callable python script
'''
import sys
import force

FILE1 = sys.argv[1]
FILE2 = sys.argv[2]
print(force.find_forces(FILE1, FILE2))
