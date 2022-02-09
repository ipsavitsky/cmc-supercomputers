make_omp:
	xlc -qsmp=omp omp/integral_threaded.c

make_mpi:
	mpixlc_r -qsmp=omp mpi/integral_parallel.c

run_omp:
	bsub < lsf_batch

run_mpi:
	llsubmit ll_batch

