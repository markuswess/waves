DEPENDS=intro.md modelling/elastic_td.md modelling/electromagnetic_td.md modelling/acoustic_td.md _config.yml _toc.yml time_integration.md time_integration/second_order.md time_integration/first_order.md
.PHONY: clean

_build: $(DEPENDS)
	jupyter-book build --all .
	./refresh_ff.sh
