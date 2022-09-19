all:
	mpost --tex=latex fig.mp
	for i in *.eps ; do epstopdf "$$i" ; done
	-rm *.eps
	for i in *.pdf ; do pdftops -eps "$$i" ; done
	for i in *.pdf ; do pdftoppm -r 500 -png "$$i" > "$${i%.pdf}.png" ; done
	-mkdir png eps pdf
	-mv *.png png/
	-mv *.eps eps/
	-mv *.pdf pdf/

clean:
	-rm *.eps *.png *.mpx *.log
