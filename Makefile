.PHONY:
all: bccheck install vignette

.PHONY:
build: clean
	R CMD build .

.PHONY:
check: clean
	Rscript -e "devtools::check()"

.PHONY:
bccheck: check
	Rscript -e "BiocCheck::BiocCheck('.')"

.PHONY:
install:
	R CMD INSTALL .

.PHONY:
vignette:
	Rscript -e "rmarkdown::render('vignettes/abseqR.Rmd', output_options = 'all')"
	rm -rf vignettes/abseqR_example/ vignettes/refined_comparison/

.PHONY:
clean:
	rm -rf vignettes/abseqR_example/ vignettes/refined_comparison/

.PHONY:
help:
	@echo -e "Available commands:"
	@echo -e "\tall      ... build, check, bccheck"
	@echo -e "\tbuild    ... builds"
	@echo -e "\tcheck    ... checks"
	@echo -e "\tbccheck  ... BiocCheck"
	@echo -e "\tinstall  ... installs"
	@echo -e "\tvignette ... builds vignette"
	@echo -e "\tclean    ... removes vignette artefacts"
