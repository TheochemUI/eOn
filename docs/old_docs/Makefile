all: eon.pdf eonInstall.pdf

eon.pdf: eon.tex
	pdflatex eon.tex
eonInstall.pdf: eonInstall.tex
	pdflatex eonInstall.tex

.PHONY : clean
clean:
	rm eon.aux eon.log eon.pdf eon.out
	rm eonInstall.aux eonInstall.log eonInstall.pdf eonInstall.out
