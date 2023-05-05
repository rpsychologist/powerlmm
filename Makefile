check:
	R CMD build --no-build-vignettes . && \
	R CMD check powerlmm_0.4.0.9000.tar.gz --ignore-vignettes && \
	code powerlmm.Rcheck/00check.log