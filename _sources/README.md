# sparse-ir


## For developer
### Set up tools

```bash
pip3 install jupyter-book ghp-import jupytext
```

### Local build

```bash
jupyter book build .
# open _build/html/index.html with your favourite browser
```

### Deployment to remote

```bash
ghp-import -n -p -f _build/html
```

### Trouble shooting

* https://github.com/executablebooks/jupyter-book/issues/1541