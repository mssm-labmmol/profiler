ITP_TO_STP := $(abspath ../../utils/gromos_itp_to_stp/gromos_itp_to_stp.py
PROFILER_GEN := $(abspath ../../profilerGen.py)
SUBDIRS = gromos amber opls

export ITP_TO_STP
export PROFILER_GEN

.PHONY: tests clean $(SUBDIRS)

tests: $(SUBDIRS)

$(SUBDIRS) :
	$(MAKE) -C $@ tests

clean:
	for subdir in $(SUBDIRS); do $(MAKE) -C $$subdir clean ; done
