DEPENDS=intro.md modelling/elastic_td.md modelling/electromagnetic_td.md modelling/acoustic_td.md _config.yml _toc.yml
.PHONY: clean

_build: $(DEPENDS)
	jupyter-book build --all .
	./refresh_ff.sh
