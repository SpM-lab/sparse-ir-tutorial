build:
	jupytext --set-formats ipynb,md:myst --sync docs/*.ipynb
	jupyter book build --all .

upload:	build
	ghp-import -n -p -f _build/html
