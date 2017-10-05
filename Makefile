# Add to latexmkrc
# Custom dependency and function for nomencl package 
# add_cus_dep( 'nlo', 'nls', 0, 'makenlo2nls' );
# sub makenlo2nls {
# system( "makeindex -s nomencl.ist -o \"$_[0].nls\" \"$_[0].nlo\"" );
#}

#add_cus_dep('glo', 'gls', 0, 'makeglo2gls');
#sub makeglo2gls {
#    system("makeindex -s '$_[0]'.ist -t '$_[0]'.glg -o '$_[0]'.gls '$_[0]'.glo");
#}



# Project-specific settings
DOCNAME = thesis

.PHONY: $(DOCNAME).pdf all clean

# The first rule in a Makefile is the one executed by default ("make"). It
# should always be the "all" rule, so that "make" and "make all" are identical.
all: $(DOCNAME).pdf

# CUSTOM BUILD RULES

# In case you didn't know, '$@' is a variable holding the name of the target,
# and '$<' is a variable holding the (first) dependency of a rule.
# "raw2tex" and "dat2tex" are just placeholders for whatever custom steps
# you might have.

%.tex: %.raw
	./raw2tex $< > $@

%.tex: %.dat
	./dat2tex $< > $@

# MAIN LATEXMK RULE

# -pdf tells latexmk to generate PDF directly (instead of DVI).
# -pdflatex="" tells latexmk to call a specific backend with specific options.
# -use-make tells latexmk to call make for generating missing files.

# -interaction=nonstopmode keeps the pdflatex backend from stopping at a
# missing file reference and interactively asking you for an alternative.

$(DOCNAME).pdf: $(DOCNAME).tex
	latexmk -pdf -pdflatex="pdflatex -interaction=scrool-mode" -use-make $(DOCNAME).tex

clean:
	latexmk -CA
