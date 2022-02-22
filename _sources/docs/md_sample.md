---
jupytext:
  text_representation:
    extension: .md
    format_name: myst
kernelspec:
  display_name: Python 3
  language: python
  name: python3
---

(md_sample)=
# Markdown sample

This pages summarizes useful information from the official document of `jupyter book`.

## Citation using bibtex
How to `` {cite:p}`Li:2020eu` `` results in {cite:p}`Li:2020eu`.
The bib file `references.bib` is located in the top directory of the repository.
Do not forget to add a list of references in the bottom of the page using the `{bibliography}` directive.


Jupyter Book uses [MathJax](http://docs.mathjax.org/) for typesetting math in your HTML book build.
This allows you to have LaTeX-style mathematics in your online content.
This page shows you a few ways to control this.

:::{seealso}
For more information about citation,
see the [official document](https://jupyterbook.org/tutorials/references.html?highlight=bibtex).
:::

:::{tip}
Some tips come here.
:::

## In-line math

To insert in-line math use the `$` symbol within a Markdown cell.
For example, the text `$this_{is}^{inline}$` will produce: $this_{is}^{inline}$.

+++

## Math blocks

You can also include math blocks for separate equations. This allows you to focus attention
on more complex or longer equations, as well as link to them in your pages. To use a block
equation, wrap the equation in either `$$` or `\begin` statements.

For example,

```latex
$$
  \int_0^\infty \frac{x^3}{e^x-1}\,dx = \frac{\pi^4}{15}
$$
```

results in:

$$
  \int_0^\infty \frac{x^3}{e^x-1}\,dx = \frac{\pi^4}{15}
$$

$$
  \int_0^\infty \frac{x^3}{e^x-1}\,dx = \frac{\pi^4}{15}
$$

+++

(md_sample:latex)=
### Latex-style math

You can use LaTeX-style math blocks.

```latex
\begin{gather*}
a_1=b_1+c_1\\
a_2=b_2+c_2-d_2+e_2
\end{gather*}

\begin{align}
a_{11}& =b_{11}&
  a_{12}& =b_{12}\\
a_{21}& =b_{21}&
  a_{22}& =b_{22}+c_{22}
\end{align}
```

which results in:

\begin{gather*}
a_1=b_1+c_1\\
a_2=b_2+c_2-d_2+e_2
\end{gather*}

\begin{align}
a_{11}& =b_{11}&
  a_{12}& =b_{12}\\
a_{21}& =b_{21}&
  a_{22}& =b_{22}+c_{22}
\end{align}

:::{seealso}
The MyST guides to [dollar math syntax](myst-parser:syntax/math), [LaTeX math syntax](myst-parser:syntax/amsmath), and [how MyST-Parser works with MathJax](myst-parser:syntax/mathjax).

For advanced use, also see how to [define MathJax TeX Macros](sphinx/tex-macros).
:::

+++

### Numbering equations

If you'd like to number equations so that you can refer to them later, use the **math directive**.
It looks like this:

````md
```{math}
:label: my_label
my_math
```
````

For example, the following code:

````md
```{math}
:label: my_label
w_{t+1} = (1 + r_{t+1}) s(w_t) + y_{t+1}
```
````

will generate

```{math}
:label: my_label
w_{t+1} = (1 + r_{t+1}) s(w_t) + y_{t+1}
```

Alternatively you can use the dollar math syntax with a prefixed label:

```md
$$
  w_{t+1} = (1 + r_{t+1}) s(w_t) + y_{t+1}
$$ (my_other_label)
```

which generates

$$
  w_{t+1} = (1 + r_{t+1}) s(w_t) + y_{t+1}
$$ (my_other_label)

:::{note}
Labels cannot start with an integer, or they won't be able to be referenced and
will throw a warning message if referenced. For example, `:label: 1` and `:label: 1eq` cannot
be referenced.
:::

### Linking to equations

If you have created an equation with a label, you can link to it from within your text
(and across pages!).

You can refer to the equation using the label that you've provided by using
the `{eq}` role. For example:

```md
- A link to an equation directive: {eq}`my_label`
- A link to a dollar math block: {eq}`my_other_label`
```

results in

- A link to an equation directive: {eq}`my_label`
- A link to a dollar math block: {eq}`my_other_label`

:::{note}
`\labels` inside LaTeX environment are not currently identified, and so cannot be referenced.
We hope to implement this in a future update (see [executablebooks/MyST-Parser#202](https://github.com/executablebooks/MyST-Parser/issues/202))!
:::

```{bibliography}
```