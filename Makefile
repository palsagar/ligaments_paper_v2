# Everything is driven from here: the solver build, the simulations,
# the analysis and the paper.

BASILISK_TARBALL := third_party/basilisk-2026-10-25.tar.gz
export BASILISK := $(CURDIR)/build/basilisk/src
QCC := $(BASILISK)/qcc
CONFIG := $(if $(filter Darwin,$(shell uname)),config.osx,config.gcc)

.PHONY: basilisk smoke linear-check paper clean

# Build the Basilisk compiler (qcc) from the archived tarball.
# The first pass fails in optional GPU subdirectories, hence `make -k`.
basilisk: $(QCC)

$(QCC): $(BASILISK_TARBALL) third_party/SHA256SUMS
	cd third_party && shasum -a 256 -c SHA256SUMS
	mkdir -p build
	tar xzf $(BASILISK_TARBALL) -C build
	cd $(BASILISK) && ln -sf $(CONFIG) config && (make -k || true) && make
	touch $(QCC)

build/ligament: sim/ligament.c $(QCC)
	$(QCC) -O2 -Wall -DFULL_INTERFACE=1 sim/ligament.c -o build/ligament -lm

# One coarse ligament (level 9, seed 1). Passes when the run reaches the
# end, the thread has broken, and the liquid volume is conserved to 0.1%.
smoke: build/ligament
	mkdir -p build/smoke
	cd build/smoke && ../ligament 9 1 2> log.txt
	awk '/^# END/ { for (i = 1; i <= NF; i++) { \
	         if ($$i == "relerr") e = $$(i+1); \
	         if ($$i == "t_pinch") p = $$(i+1) } \
	       ok = (e < 1e-3 && e > -1e-3 && p > 0) } \
	     END { exit !ok }' build/smoke/log.txt
	@grep '^# END' build/smoke/log.txt

# Linear-theory model on 20,000 random surfaces (about 3 minutes).
linear-check:
	uv run python theory/linear_check.py 20000 400

paper:
	cd paper && latexmk -pdf -outdir=build main.tex

clean:
	rm -rf build paper/build
