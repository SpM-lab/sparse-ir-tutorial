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

## Docker on VSCode
1. Make sure you have Docker and VSCode installed.

2. Install the ``Remote - Containers`` in VSCode Extensions.

3. Go to the repository and open it in VS Code.

```
cd /path/to/this/repository

code .
```

4. To use Docker with VS Code, execute the following command.

```
ln -s .dev/devcontainer .devcontainer
```

5. Press the green mark at the bottom left and press `` Reopen in Container`` from the command palette.
   After the build is finished, you can enter the Docker container.

## References
You can add references to `references.bib`.

## Trouble shooting

* https://github.com/executablebooks/jupyter-book/issues/1541
