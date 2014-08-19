#!/usr/bin/env python
#
# vim: set ts=4 noet:
#
import sys
import numpy as np
import os
import datetime
import drdcommon

#
# tool <ref.fa> <chrm> <truth_calls> <our_calls>
# 1. load chrm
# 2. Load truth
# 3. load calls
#
# FFFFFFFFFFFFFFTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
# 000000000000000777777777777776666666666660000000000
# 000000000000066666666666777777766666666666600000000
def log(msg):
	sys.stderr.write("%s >> %s\n" % (datetime.datetime.now(), msg))

def error(msg):
	sys.stderr.write("ERROR: %s\n" % msg)
	sys.exit(1)

def usage():
	sys.stderr.write("Usage: tool <ref.fa> <chrm> <truth.bed> <our_calls_canavar.bed>\n")
	sys.exit(1)

def prepareDS():
	_max = 250000000 # chrm1 is the biggest chrm
	# ds[0] ; 0: N 1: base
	# ds[1] ; 0: No call || cnv call value
	# ds[2] ; 0: No call || cnv call value
	log("Creating DS...")
	return np.zeros((3, _max), dtype=int)

def loadChrm(ds, ref, chrm):
	if not os.path.exists(ref):
		error("Cannot find reference file: %s", ref)

	log("Reading reference genome chrm: %s" % chrm)
	i = 1
	for l in drdcommon.xopen(ref):
		l = l.strip()
		if i == 1 and l[0] == '>' and l[1:] == chrm:
			continue
		if i > 1 and l[0] == '>':
			break

		for bp in l:
			if bp.upper() != 'N':
				ds[0][i] = 1
			if i % 10000000 == 0:
				sys.stderr.write("MEM: %s nbp: %s\r" % (drdcommon.memory_usage(), i))
			i += 1
	log("\n%s bp read." % i)

def loadCalls(ds, fn, idx, chrm):
	log("Loading calls from %s; idx=%s" % (fn, idx))
	chrm_found = False
	nbp = 0
	for l in drdcommon.xopen(fn):
		c, start, end, cnv = l.strip().split()
		if c == chrm:
			chrm_found = True
			for i in range(int(start), int(end)+1):
				if nbp % 1000000 == 0:
					sys.stderr.write("MEM: %s nbp: %s\r" % (drdcommon.memory_usage(), nbp))
				ds[idx][i] = round(float(cnv))
				nbp += 1

	if not chrm_found:
		error("\nCould not find chrm in file. Bailing out.")
	log("\n%s bp loaded" % nbp)

def computeMetrics(ds):
	# Distribution of cnvs
	# Distribution of window sizes
	#
	None

def doWork(ref, chrm, truth, canavar):
	ds = prepareDS()
	loadChrm(ds, ref, chrm)
	i_truth, i_canavar = 1, 2
	loadCalls(ds, truth, i_truth, chrm)
	loadCalls(ds, canavar, i_canavar, chrm)
	computeDistCopyNumbers(ds)
	log("MEM: %s" % drdcommon.memory_usage())
	log('Done')

def main():
	if len(sys.argv) != 5:
		usage()
	else:
		ref, chrm, truth, canavar = sys.argv[1:]
		doWork(ref, chrm, truth, canavar)

if __name__ == "__main__":
	main()
