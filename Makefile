build:
	jupyter book build --all .

update_header:
	find src -name "*.md" |xargs ./bin/update_metadata.py

upload:	build
	ghp-import -n -p -f _build/html

clean:
	rm -rf _build
