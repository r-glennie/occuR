all: occu_tmb.so 

occu_tmb.so: occu_tmb.cpp
	Rscript --vanilla -e "TMB::compile('occu_tmb.cpp','-O0 -g -Wno-ignored-attributes', safebounds=FALSE)"
	
clean:
	rm -rf *o 
