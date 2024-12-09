all: log.txt

update-sources:
	-touch run.py

log.txt: run.py
	make clean
	bash run.sh

clean:
	-rm -rfv locks results* __pycache__ log.full.txt log.check.txt*
