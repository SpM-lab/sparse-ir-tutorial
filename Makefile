build:
	rm -rf docs
	cp -r docs_template docs
	jupyter book build .