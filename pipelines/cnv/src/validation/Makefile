THIS_MAKE := $(lastword $(MAKEFILE_LIST))
MAKE_DIR := $(dir $(THIS_MAKE))

.PHONY: publish

results.%:
	read a b c d <<< $$(echo "$@" | ruby  -ne 'print $$_.match(/results.(\d+)-(\d+)-(\d+)-(\d+)/)[1..4].join(" ")');\
	$(MAKE_DIR)/run-intersect.sh $$a $$b $$c $$d | tee $@


publish: Validate_CNV_pipeline.ipynb
	IP=`ls *.ipynb`; scp $$IP drio@davidr.io:public_html/notebooks; ssh drio@davidr.io "cd public_html/notebooks;ipython nbconvert $$IP"
	@echo "http://davidr.io/~drio/notebooks/Validate_CNV_pipeline.html"
