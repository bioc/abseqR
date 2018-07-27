.PHONY:
vignette:
	Rscript -e "rmarkdown::render('vignettes/abseqR.Rmd', output_options = 'all')"

.PHONY:
check:
	R CMD build
	R CMD check
	Rscript -e "library(BiocCheck);pathToPkg <- '.';BiocCheck(pathToPkg)"
