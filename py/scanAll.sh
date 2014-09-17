/gpfs/data/jvbq85/opt/bin/mybsub job.msub A2 0 1 1 1 300
/gpfs/data/jvbq85/opt/bin/mybsub job.msub A2 1 1 0 1 300
/gpfs/data/jvbq85/opt/bin/mybsub job.msub A2 0 0 0 1 300
/gpfs/data/jvbq85/opt/bin/mybsub job.msub A2 1 0 0 1 300
for halo in B2 C2 D2 E2
do
/gpfs/data/jvbq85/opt/bin/mybsub job.msub $halo 0 1 0 1 300
/gpfs/data/jvbq85/opt/bin/mybsub job.msub $halo 1 1 0 1 300
/gpfs/data/jvbq85/opt/bin/mybsub job.msub $halo 0 0 0 1 300
/gpfs/data/jvbq85/opt/bin/mybsub job.msub $halo 1 0 0 1 300
done
# for halo in A2 B2 C2 D2 E2
# do
# /gpfs/data/jvbq85/opt/bin/mybsub job.msub $halo 0 1 0 10 100
# /gpfs/data/jvbq85/opt/bin/mybsub job.msub $halo 1 1 0 10 100
# /gpfs/data/jvbq85/opt/bin/mybsub job.msub $halo 1 0 0 10 100
# /gpfs/data/jvbq85/opt/bin/mybsub job.msub ${halo}N 1 0 0 10 100
# done