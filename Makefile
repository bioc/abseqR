.PHONY: all build check bccheck full_check install vignette clean help

all: full_check

build: clean
	R CMD build .

check: clean
	Rscript -e "devtools::check()"

bccheck:
	Rscript -e "BiocCheck::BiocCheck('.')"

full_check: check bccheck

install:
	R CMD INSTALL .

vignette:
	Rscript -e "rmarkdown::render('vignettes/abseqR.Rmd', output_options = 'all')"
	rm -rf vignettes/ex/ vignettes/refined_comparison/

clean:
	rm -rf vignettes/ex/ vignettes/refined_comparison/ vignettes/abseqR.html

help:
	@echo -e "Available commands:"
	@echo -e "\tall         ... build, check, bccheck"
	@echo -e "\tbuild       ... builds"
	@echo -e "\tfull_check  ... check and bccheck"
	@echo -e "\tcheck       ... checks"
	@echo -e "\tbccheck     ... BiocCheck"
	@echo -e "\tinstall     ... installs"
	@echo -e "\tvignette    ... builds vignette"
	@echo -e "\tclean       ... removes vignette artefacts"
