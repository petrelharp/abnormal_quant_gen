.PHONY: all, clean

all: writeup.pdf

writeup.pdf : refs.bib review-responses.tex

writeup-only.pdf : writeup.pdf
	pdfjam --outfile $@ $< 1-27

cover-letter.pdf : writeup.pdf
	pdfjam --outfile $@ $< 28

response-to-reviewers.pdf : writeup.pdf
	pdfjam --outfile $@ $< 29-

diff-to-submission.pdf : writeup-diffb8f1c941105b76b6e99932948262386888d2d172.pdf
	cp $< $@

diff-to-first-submission.pdf : writeup-diffc8561ba72d021df4d4f08e4e97b5b023853053d5.pdf
	cp $< $@

writeup-diff%.tex : writeup.tex refs.bib review-responses.tex
	latexdiff-git -r $* --force writeup.tex

clean: 
	-rm *.aux *.log *.bbl *.blg

%.pdf : %.tex %.bbl
	while ( pdflatex -shell-escape $<;  grep -q "Rerun to get" $*.log ) do true ; done
	touch $*.bbl
	touch $@

%.aux : %.tex
	-pdflatex -shell-escape $<

%.bbl : %.aux refs.bib
	bibtex $<

%.png : %.pdf
	convert -density 300 $< -flatten $@

%.pdf : %.svg
	inkscape $< --export-area-drawing --export-filename=$@
	# chromium --headless --no-pdf-header-footer --print-to-pdf=$@ $<
	# ./svg2pdf.sh $< $@

%.pdf : %.eps
	# inkscape $< --export-filename=$@
	epspdf $<

%.pdf : %.ink.svg
	inkscape $< --export-filename=$@

