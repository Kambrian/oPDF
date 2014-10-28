# /gpfs/data/jvbq85/opt/bin/mybsub job.msub A2 0 1 1 1 300
# /gpfs/data/jvbq85/opt/bin/mybsub job.msub A2 1 1 0 1 300
# /gpfs/data/jvbq85/opt/bin/mybsub job.msub A2 0 0 0 1 300
# /gpfs/data/jvbq85/opt/bin/mybsub job.msub A2 1 0 0 1 300
for halo in A2 B2 C2 D2 E2
do
#   for w in 0 1
#   do
# 	for proxy in r #L
# 	do
# 	  for percent in 0.2 0.4 0.6
# 	  do
#order: clean, template, weight, contour, proxy, percent
		  /gpfs/data/jvbq85/opt/bin/mybsub job.msub $halo 1 1 0 0 #r $percent #$proxy $percent #clean, tmpl
# 		  /gpfs/data/jvbq85/opt/bin/mybsub job.msub $halo 1 0 0 0 r $percent #$proxy $percent #clean, NFW
	  # 	/gpfs/data/jvbq85/opt/bin/mybsub job.msub $halo 0 1 $w 0 #dirty, tmpl
	  # 	/gpfs/data/jvbq85/opt/bin/mybsub job.msub $halo 0 0 $w 0 #dirty, NFW
	  # 	/gpfs/data/jvbq85/opt/bin/mybsub job.msub $halo 2 1 $w 0 #branchclean, tmpl
	  # 	/gpfs/data/jvbq85/opt/bin/mybsub job.msub $halo 2 0 $w 0 #branchclean, NFW
# 	  done
# 	done
#   done
done
# for halo in A2 B2 C2 D2 E2
# do
# /gpfs/data/jvbq85/opt/bin/mybsub job.msub $halo 0 1 0 10 100
# /gpfs/data/jvbq85/opt/bin/mybsub job.msub $halo 1 1 0 10 100
# /gpfs/data/jvbq85/opt/bin/mybsub job.msub $halo 1 0 0 10 100
# /gpfs/data/jvbq85/opt/bin/mybsub job.msub ${halo}N 1 0 0 10 100
# done