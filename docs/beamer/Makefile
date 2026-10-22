TEX=na64mini_meiga_results.tex
PDF=na64mini_meiga_results.pdf

.PHONY: figures geometry pdf clean

figures:
	python3 generate_beamer_figures.py

geometry:
	bash render_meiga_geometry.sh || true

pdf: figures
	@if command -v latexmk >/dev/null 2>&1; then \
		latexmk -pdf -interaction=nonstopmode -halt-on-error $(TEX); \
	else \
		pdflatex -interaction=nonstopmode -halt-on-error $(TEX); \
		pdflatex -interaction=nonstopmode -halt-on-error $(TEX); \
	fi

clean:
	rm -f *.aux *.bbl *.bcf *.blg *.fdb_latexmk *.fls *.log *.nav *.out *.run.xml *.snm *.toc *.vrb
