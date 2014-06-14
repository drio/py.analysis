#!/bin/bash
# vim :set ts=2 noet:


max_win=5000
max_calls=100
max_truth=10

n_cores=16
i=1

RAND=`strings /dev/urandom | grep -o '[[:alnum:]]' | head -n 50 | tr -d '\n'; echo`
rm -rf all_inter
mkdir all_inter
for _o in `seq $i $n_cores`;do
  touch all_inter/$RAND.$_o
done

for ws in `seq 500 500 $max_win`;do
	for c_min in `seq 3 2 $max_calls`;do
		for c_max in `seq $max_calls -2 3`;do
			[ `expr $c_max - $c_min` -lt 2 ] && continue
			for t_min in `seq 3 1 $max_truth`;do
				for t_max in `seq $max_truth -1 3`;do
					[ `expr $t_max - $t_min` -lt 2 ] && continue
					./src/run-intersect.sh $c_min $c_max $t_min $t_max $ws inputs/french.bed inputs/hg18_drio_canavar.bed 2>/dev/null >> all_inter/$RAND.$i &
          if [ $i == $n_cores ];then
            i=0
            wait
          fi
          i=$[$i+1]
				done
			done
		done
	done
done

