define help
TARGETS
	help                  : this.
	update                : update repo against ardmore.
	watch WATCH_FOR=files : watch for this files and update.
	forward               : create socket against remote ipython instance in ardmore
endef

export help

help:
	@echo "$$help"

update:
	@cd /Users/drio/dev/py.analysis;\
	rsync -avz --progress * is-ardmore:/stornext/snfs6/rogers/drio_scratch/dev/py.analysis/

watch:
	@echo "Watching for: $(WATCH_FOR)"
	@filewatcher '$(WATCH_FOR)' "make update"

forward:
	@(_PORT=7989;\
	_PID_TO_KILL=`ps aux | grep ssh | grep $$_PORT | grep -v grep | awk '{print $$2}' | xargs echo`;\
	[ ".$$_PID_TO_KILL" != "." ] && echo "kill -9 $$_PID_TO_KILL") ;\
	echo "ssh -N -L 7000:localhost:$$_PORT is &"

.PHONY: watch help update forward