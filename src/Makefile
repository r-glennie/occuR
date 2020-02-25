all: occu_tmb.so 

occu_tmb.so: occu_tmb.cpp
	Rscript --vanilla -e "TMB::compile('occu_tmb.cpp','-O0 -g', safebounds=FALSE)"
	
clean:
	rm -rf *o 
