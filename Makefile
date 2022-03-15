update_header:
	find src -name "*.md" |xargs ./bin/update_metadata.py

build:	update_header
	jupyter book build --all .

upload:	build
	ghp-import -n -p -f _build/html

clean:
	rm -rf _build
