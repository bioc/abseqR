.PHONY:
all: bccheck install vignette

.PHONY:
build: clean
	R CMD build .

.PHONY:
check: build
	R CMD check .

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
	@echo "all      ... builds, checks, BiocCheck"
	@echo "build    ... builds"
	@echo "check    ... checks"
	@echo "bccheck  ... BiocCheck"
	@echo "install  ... installs"
	@echo "vignette ... builds vignette"
	@echo "clean    ... removes vignette artefacts"
