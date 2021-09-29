
main.pdf: WOR.tex refs.bib auths.tex main.png WellcomeOR_styles.sty WellcomeOR_logo_black.pdf
	pdflatex $<
	bibtex $(subst .tex,,$<)
	pdflatex $<
	pdflatex $<

R = $(strip Rscript $^ $(1) $@)

REFSPEC ?= 0.99
REFCOV  ?= 0.2

main.png: fig_main.R
	$(call R,${REFSPEC} ${REFCOV})