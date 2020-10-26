refdih = 3

.PHONY: tests clean

phe-ala.stp: phe-ala.itp
	echo "$(refdih)\n0\n1\n0" | $(ITP_TO_STP) -f $< > /dev/null 2>&1 

phe-ala-gen.dat: phe-ala.stp phe-ala.gro
	$(PROFILER_GEN) -c phe-ala.gro -t phe-ala.stp -dr 0 10 360 -dk 5000 -op phe-ala-gen -min 1 -nsteps 500000 -dele 1e-04 > /dev/null 2>&1

compare.o: phe-ala-gen.dat reference.dat compare.py
	./compare.py > compare.o

tests: compare.o

clean: 
	rm -rf $(wildcard phe-ala-gen*) phe-ala.stp compare.o
