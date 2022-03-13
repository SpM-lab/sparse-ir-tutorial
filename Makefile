build:
	rm -rf build
	cp -r src build
	jupytext --set-kernel - build/*_py.md
	jupytext --set-kernel - build/*_jl.md
	jupyter book build --all build

upload:	build
	ghp-import -n -p -f build/_build/html

clean:
	rm -rf build
