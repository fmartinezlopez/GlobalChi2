.PHONY : main clean

main:
	latexmk -cd -pdf -bibtex -output-directory='.latexmk' paper/main.tex
	rm -f main.pdf
	cp paper/.latexmk/main.pdf .

clean:
	rm -rf paper/.latexmk main.pdf
