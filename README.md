# sparse-ir


## For developer
### Set up tools

```bash
pip3 install jupyter-book ghp-import jupytext
```

### Local build
The following command update all markdown files paired with notebooks and build html page. 

```bash
make build
```

### Deployment to remote

```bash
make upload
```

### References
You can add references to `references.bib`.

### Jupytext

Convert a MyST Markdown file to a notebook:

```bash
jupytext --to ipynb some_file.md
```

Convert a notebook to a MyST Markdown file:

```bash
jupytext --to md:myst some_file.ipynb
```

### Trouble shooting

* https://github.com/executablebooks/jupyter-book/issues/1541
