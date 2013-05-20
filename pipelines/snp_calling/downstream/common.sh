min_snp_qual=50

if [ "$seq_type" == "wes" ];then
  f_cons_to_remove="NONE|INTERGENIC|INTRON|UPSTREAM|DOWNSTREAM"
  drop_this="grep -v -P $f_cons_to_remove"
  gen_size=34000000
  export min=10
  export max=300
else
  drop_this="cat -"
  gen_size=2900000000
  export min=2
  export max=100
fi

export min_num_samples=2

if [ "$seq_type" == "wes" ];then
  gen_size=34000000
fi
if [ "$seq_type" == "wgs" ];then
  gen_size=2900000000
fi

#arr=("California" "New England" "Oregon" "Wisconsin" "Yerkes")

