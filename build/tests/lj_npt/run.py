from subprocess import check_output
import os

check = 320
energy = float(check_output('../../mpmc test.input | grep "volume" | grep -o "= [-.0-9]*" | grep -o [-.0-9]* | tail -1',shell=True))

if (abs(energy-check)<=40):
	print os.getcwd().split(os.sep)[-1] + ": TEST PASSED"
else:
	print os.getcwd().split(os.sep)[-1] + ": TEST FAILED"