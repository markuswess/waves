DEPENDS=intro.md modelling/time-domain.md _config.yml _toc.yml
.PHONY: clean

_build: $(DEPENDS)
	jupyter-book build --all .
	./refresh_ff.sh
