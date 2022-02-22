build:
	jupyter book build .

upload:	build
	ghp-import -n -p -f _build/html
