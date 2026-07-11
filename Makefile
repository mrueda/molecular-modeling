CXX ?= g++
CXXFLAGS ?= -std=c++11 -O2 -Wall -Wextra -pedantic
PYTHON ?= python3
PERL ?= perl

.PHONY: all build check clean run-md-cpp run-md-perl run-docking-cpp run-docking-perl

all: build

build: md_simulation/cpp/md_simulation_h20 docking/cpp/peptide_docking

md_simulation/cpp/md_simulation_h20: md_simulation/cpp/md_simulation_h20.cpp
	$(CXX) $(CXXFLAGS) -o $@ $<

docking/cpp/peptide_docking: docking/cpp/peptide_docking.cpp
	$(CXX) $(CXXFLAGS) -o $@ $<

check:
	$(PERL) -c md_simulation/perl/md_simulation_h20.pl
	$(PERL) -c docking/perl/peptide_docking.pl
	$(PERL) -c misc/discrete_fourier_transform/discrete_fourier_transform.pl
	$(PERL) -c misc/discrete_fourier_transform/fast_fourier_transform.pl
	$(PYTHON) -m py_compile md_simulation/plots/plot_3d.py docking/plots/plot_docking.py
	$(CXX) $(CXXFLAGS) -o /tmp/md_simulation_h20_check md_simulation/cpp/md_simulation_h20.cpp
	$(CXX) $(CXXFLAGS) -o /tmp/peptide_docking_check docking/cpp/peptide_docking.cpp
	$(PYTHON) -m unittest discover -s tests

run-md-cpp: md_simulation/cpp/md_simulation_h20
	./md_simulation/cpp/md_simulation_h20 > md_simulation/cpp/simulation_output.txt

run-md-perl:
	$(PERL) md_simulation/perl/md_simulation_h20.pl > md_simulation/perl/simulation_output.txt

run-docking-cpp: docking/cpp/peptide_docking
	cd docking/cpp && ./peptide_docking

run-docking-perl:
	cd docking/perl && $(PERL) peptide_docking.pl

clean:
	$(RM) md_simulation/cpp/md_simulation_h20
	$(RM) docking/cpp/peptide_docking
	$(RM) md_simulation/cpp/simulation_output.txt
	$(RM) md_simulation/perl/simulation_output.txt
	$(RM) docking/cpp/docking_trajectory.xyz
	$(RM) docking/perl/docking_trajectory.xyz
	$(RM) -r md_simulation/plots/__pycache__
	$(RM) -r docking/plots/__pycache__
	$(RM) -r tests/__pycache__
