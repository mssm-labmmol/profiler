PROFILER_OPT := $(abspath ../../profilerOpt.py)

INPUTS = $(wildcard output*.inp)
OUTPUTS = $(INPUTS:.inp=.o)
COMPARES = $(subst output,compare,$(OUTPUTS))

.PHONY: clean tests prepare compare

output_%.o: output_%.inp input.stp input.g96 input.dat
	mkdir -p test_$* ; $(PROFILER_OPT) -np 1 -c input.g96 -r input.dat -t input.stp -i $< -op test_$*/output > output_$*.out 2>&1
	grep "Finished execution" output_$*.out > $@

# specifics
output_TORSTYPE.2.o: output_TORSTYPE.2.inp input-Ryck.stp input-Ryck.g96 input-Ryck.dat
	mkdir -p test_TORSTYPE.2 ; $(PROFILER_OPT) -np 1 -c input-Ryck.g96 -r input-Ryck.dat -t input-Ryck.stp -i $< -op test_TORSTYPE.2/output > output_TORSTYPE.2.out 2>&1
	grep "Finished execution" output_TORSTYPE.2.out > output_TORSTYPE.2.o

output_PAIRS.o: output_PAIRS.inp input-Pairs.stp input.g96 input.dat
	mkdir -p test_PAIRS ; $(PROFILER_OPT) -np 1 -c input.g96 -r input.dat -t input-Pairs.stp -i $< -op test_PAIRS/output > output_PAIRS.out 2>&1
	grep "Finished execution" output_PAIRS.out > output_PAIRS.o

tests: $(OUTPUTS) compare

clean:
	rm -rf $(wildcard output*) test_* prepare.o compare_*.o

prepare.o: sweeper.py
	./sweeper.py && cp output_MALG.0.inp output_PAIRS.inp > $@

prepare: prepare.o

# comparisons
compare_%.o: output_%.o
	diff -Enwbr test_$*/ reference_$*/ && echo "" > $@

compare: $(COMPARES)

reference: $(OUTPUTS)
	find . -name "test_*" -type d -exec bash -c 'foo={} ; cp -r $${foo} $${foo/test/reference}' \;

