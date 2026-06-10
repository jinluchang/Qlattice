# `qlat_utils.qplot` — Gnuplot-Based Plotting Utilities

Source: `qlat-utils/qlat_utils/qplot.py`

> **Note:** Update this document when updating the source file.

## Outline

1. [Overview](#overview)
2. [Data Table I/O](#data-table-io)
3. [Plot Generation](#plot-generation)
4. [Display Helpers](#display-helpers)
5. [Module-Level Configuration](#module-level-configuration)
6. [Examples](#examples)

---

## Overview

`qplot` provides a lightweight interface for generating publication-quality
plots using **gnuplot** as the backend.  It handles:

- Saving and loading numeric data tables (real and complex) as text files.
- Generating gnuplot plotfiles and Makefiles in temporary directories.
- Running `make` + `gnuplot` + `pdflatex` + `pdftoppm` to produce PNG/PDF
  output.
- Displaying images inline in Jupyter notebooks.

The typical workflow is:

1. Prepare data as a dict of 2-D arrays (`dts`).
2. Write gnuplot commands (`cmds`) and plot lines (`lines`).
3. Call `plot_save()` to render to a file, or `plot_view()` to display inline.

---

## Data Table I/O

### `show_datatable`

```python
q.show_datatable(arr, *, is_return_list_of_string=False)
```

Convert a 2-D array to a whitespace-separated string (or list of lines).
Supports real and complex numbers.

### `save_datatable`

```python
q.save_datatable(arr, fn, *, is_directory_exist=False)
```

Save a 2-D numpy array to a text file.  Creates parent directories
automatically.

### `read_datatable`

```python
q.read_datatable(lines)
```

Parse a string or list of strings into a list of lists of numbers
(float or complex).  Complex numbers are represented as `real imagi`
(e.g. `1.5 2.0i`).

### `load_datatable`

```python
q.load_datatable(fn)
```

Load a data table from a text file.

### `azip`

```python
q.azip(vec, *vecs)
```

Zip multiple 1-D arrays column-wise into a 2-D array, truncating to the
shortest length.  Useful for assembling `(x, y, yerr)` tables.

---

## Plot Generation

### `plot_save`

```python
q.plot_save(
    fn=None, dts=None, cmds=None, lines=None, *,
    is_run_make=True, is_display=False, is_verbose=False,
    display_width=None,
)
```

Generate a plot and optionally render it.

| Parameter | Description |
|---|---|
| `fn` | Output filename (e.g. `"plot.png"`), or `None` for a temporary plot |
| `dts` | Dict of `{ "table.txt": [[x, y, yerr], ...], ... }` |
| `cmds` | List of gnuplot commands (e.g. `["set xlabel '$x$'"]`) |
| `lines` | List of gnuplot plot lines (e.g. `["plot", "'table.txt' w yerrorb"]`) |
| `is_run_make` | If `True`, run the Makefile to produce the plot |
| `is_display` | If `True`, display inline in Jupyter |

**Returns:** Path to the generated PNG, or to the source directory if
`is_run_make=False`.

When called with all `None` arguments, generates a demo cosine plot and
prints the full call signature for easy copy-paste customisation.

### `plot_view`

```python
q.plot_view(
    fn=None, dts=None, cmds=None, lines=None, *,
    is_verbose=False, display_width=None,
)
```

Convenience wrapper around `plot_save` with `is_display=True`.  Renders the
plot and displays it inline in a Jupyter notebook.

---

## Display Helpers

### `display_img`

```python
q.display_img(fn, *, width=None)
```

Display an image file in a Jupyter notebook using `IPython.display`.

---

## Module-Level Configuration

### `gnuplot_png_density`

```python
q.gnuplot_png_density  # default: 500
```

DPI for the `pdftoppm` PNG rasterisation step.

### `gnuplot_plotfile_header`

```python
q.gnuplot_plotfile_header  # default: ["set terminal epslatex standalone color clip lw 3"]
```

Header lines prepended to every generated gnuplot plotfile.

### `plot_save_display_width`

```python
q.plot_save_display_width  # default: None
```

Default display width (in pixels) used by `plot_view` when showing images
in Jupyter.  Can be overridden per-call via the `display_width` parameter.

---

## Examples

### Simple scatter plot with error bars

```python
import qlat_utils as q
import numpy as np

x = np.linspace(0, 2 * np.pi, 30)
y = np.sin(x)
yerr = 0.1 * np.ones_like(x)

q.plot_save(
    fn="sine_plot.png",
    dts={"sine.txt": q.azip(x, y, yerr)},
    cmds=[
        "set xlabel '$x$'",
        "set ylabel '$y$'",
        "set title 'Sine function'",
    ],
    lines=[
        "plot [:] [-1.5:1.5]",
        "sin(x) w l t '$\\sin(x)$'",
        "'sine.txt' u 1:2:3 w yerrorb t 'data'",
    ],
)
```

### Display inline in Jupyter

```python
import qlat_utils as q

q.plot_view(
    fn="my_plot.png",
    dts=dts,
    cmds=cmds,
    lines=lines,
    display_width=600,
)
```

### Load a saved data table

```python
import qlat_utils as q

data = q.load_datatable("sine.txt")
# data is a list of lists: [[x0, y0, yerr0], [x1, y1, yerr1], ...]
```
