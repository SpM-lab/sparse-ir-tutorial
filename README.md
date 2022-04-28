# sparse-ir


## Set up tools

```bash
pip3 install jupyter-book ghp-import jupytext
```

## How to build html

1. Write a jupyter notebook. Your notebook must be placed under the  `src` directory with a name `*_py.ipynb` (python) or `*_jl.ipyb` (julia).

2. Convert the notebook to a MyST markdown file.

```bash
jupytext --to md:myst notebook_py.ipynb
```

3. Make html files

The following command builds html files. 

```bash
make build
```

## Update existing MyST Markdown file

1. Convert a MyST markdown file to a notebook:

```bash
jupytext --to ipynb notebook_py.md
```

2. Update the notebook using jupyter notebook/lab

3. Sync the markdown file and the notebook.

```bash
jupytext --sync --to md:myst notebook_py.ipynb
```


## How to commit a new/updated Markdown file

Before commiting a new MyST Markdown file, update the header of all Markdown files:

```bash
make update_header
```

Then, you will be ready to commit updated files.


## Upload html files

```bash
make upload
```

## References
You can add references to `references.bib`.

## Trouble shooting

* https://github.com/executablebooks/jupyter-book/issues/1541
