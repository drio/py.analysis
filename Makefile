define help
TARGETS
	help                  : this.
	update                : update repo against ardmore.
	watch WATCH_FOR=files : watch for this files and update.
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

.PHONY: watch help update
