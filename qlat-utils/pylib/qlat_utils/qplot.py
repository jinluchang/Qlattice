import numpy as np
import os
import sys
import shutil
import tempfile
import subprocess

def mk_file_dirs(fn):
    path = os.path.dirname(fn)
    if path != "":
        os.makedirs(path, exist_ok=True)

def show_number(x):
    if isinstance(x, complex):
        return f"{x.real} {x.imag}i"
    elif isinstance(x, float):
        return f"{x}"
    elif isinstance(x, float):
        return f"{x}"
    else:
        return str(x)

def show_vector(vec):
    return " ".join([ show_number(x) for x in vec ])

def show_datatable(arr, *, is_return_list_of_string = False):
    lines = [ show_vector(vec) + "\n" for vec in arr ]
    if is_return_list_of_string:
        return lines
    return "".join(lines)

def touch_file(fn, content = "", *, is_directory_exist = False):
    if not is_directory_exist:
        mk_file_dirs(fn)
    with open(fn, "w") as f:
        if isinstance(content, str):
            f.write(content)
        else:
            for s in content:
                f.write(s)

def save_datatable(arr, fn, *, is_directory_exist = False):
    """save_datatable(arr, fn), arr is (numpy) 2-D array, fn is file path name."""
    touch_file(fn,
            show_datatable(arr, is_return_list_of_string = True),
            is_directory_exist = is_directory_exist)

def read_number(s):
    if s[-1] == "i":
        return complex(0, float(s[:-1]))
    else:
        return float(s)

def read_vector(line):
    ss = line.split()
    fs = [ read_number(s) for s in ss if s != "" ]
    n = len(fs)
    xs = []
    i = 0
    while i < n:
        if i + 1 < n and isinstance(fs[i + 1], complex):
            xs.append(fs[i] + fs[i + 1])
            i += 2
        else:
            xs.append(fs[i])
            i += 1
    return xs

def read_datatable(lines):
    # return list of list of numbers (float or complex)
    if isinstance(lines, str):
        return read_datatable(lines.splitlines())
    return [ read_vector(line) for line in lines ]

def load_datatable(fn):
    with open(fn) as f:
        return read_datatable(f)

def azip(vec, *vecs):
    size_list = map(len, vecs)
    size_min = len(vec)
    for size in size_list:
        if size < size_min:
            size_min = size
    return np.array([ v[:size_min] for v in [ vec, ] + list(vecs) ]).transpose()

gnuplot_png_density = 500

def mk_tmp_dir():
    return tempfile.mkdtemp(suffix = ".dir", prefix="pyplot.")

def mk_convert_sh():
    return "\n".join([
        "for i in *.mp ; do",
        "fn=${i%.mp}",
        "rm mpost-job.* 2>&1",
        "TEX=latex mpost -jobname mpost-job $fn.mp",
        "for i in mpost-job.? ; do",
        "echo \"$i\"",
        "mv \"$i\" \"$fn\"-\"${i#mpost-job.}\".eps",
        "done",
        "done",
        "",
        ])

valid_fn_chars = ("abcdefghijklmnopqrstuvwxyz"
        + "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
        + "0123456789"
        + " ,.+-_=;:[]?{}"
        )

def check_fn(fn):
    if isinstance(fn, str):
        for c in fn:
            if c not in valid_fn_chars:
                return False
        return True
    return False

def get_plot_name(fn):
    assert fn[-4:] == ".png"
    assert check_fn(fn)
    name = fn[:-4]
    return name

def mk_makefile(fn = None):
    # fn is the target file name, e.g. plot.pdf or plot.png
    if fn is not None:
        name = get_plot_name(fn)
        target = "install-png"
    else:
        name = "plot-0"
        target = "png"
    return "\n".join([
        f"all: {target}",
        "",
        "gnuplot:",
        "\tgnuplot plotfile",
        "",
        "mpost: gnuplot",
        "\tbash ./convert.sh",
        "",
        "pdf: mpost",
        "\tepstopdf plot-0.eps",
        f"\tmv plot-0.pdf tmp ; mv tmp '{name}.pdf'",
        "",
        "png: pdf",
        f"\tpdftoppm -r {gnuplot_png_density} -png '{name}.pdf' > plot-0.png",
        f"\tmv plot-0.png tmp ; mv tmp '{name}.png'",
        "",
        "install-pdf: pdf",
        f"\tmv '{name}.pdf' ../'{name}.pdf'",
        "",
        "install-png: png",
        f"\tmv '{name}.png' ../'{name}.png'",
        "",
        "install: install-pdf install-png",
        "",
        ])

def mk_plotfile(plot_cmds, plot_lines):
    plot_prefix = [
            "set term mp color latex prologues 3 amstex",
            "set output 'plot.mp'",
            ]
    plot = plot_lines[0] + " \\\n    " + ", \\\n    ".join(plot_lines[1:])
    return "\n".join(plot_prefix + [ "", ] + plot_cmds + [ "", plot, "", ])

def populate_pyplot_folder(
        path,
        *,
        fn = None,
        dict_datatable = None,
        plot_cmds = None,
        plot_lines = None,
        ):
    if dict_datatable is None:
        dict_datatable = {}
    if plot_cmds is None:
        plot_cmds = []
    if plot_lines is not None:
        touch_file(os.path.join(path, "plotfile"),
                mk_plotfile(plot_cmds, plot_lines))
    touch_file(os.path.join(path, "convert.sh"), mk_convert_sh())
    touch_file(os.path.join(path, "Makefile"), mk_makefile(fn))
    for key, dt in dict_datatable.items():
        assert key[-4:] == ".txt"
        assert check_fn(key)
        save_datatable(dt, os.path.join(path, key))

def qremove_all(path):
    return shutil.rmtree(path, ignore_errors = True)

def mk_pyplot_folder(path = None):
    if path is None:
        path = mk_tmp_dir()
    else:
        assert isinstance(path, str)
        assert path.endswith(".pyplot.dir")
        qremove_all(path)
        os.makedirs(path)
    return path

def display_img(fn):
    from IPython.display import Image, display
    display(Image(filename = fn))

def plot_save(
        fn = None,
        dts = None,
        cmds = None,
        lines = None,
        *,
        is_run_make = True,
        is_display = False,
        is_verbose = False,
        ):
    # fn is full name of the plot or None
    # dts is dict_datatable, e.g. { "table.txt" : [ [ 0, 1, ], [ 1, 2, ], ], }
    # cmds is plot_cmds, e.g. [ "set key rm", "set size 1.0, 1.0 ", ]
    # lines is plot_lines, e.g. [ "plot", "x", ]
    target = fn
    target_fn = None
    path = None
    if target is not None:
        target_fn = os.path.basename(target)
        path = os.path.join(os.path.dirname(target),
                get_plot_name(target_fn) + ".pyplot.dir")
    path = mk_pyplot_folder(path)
    if cmds is None:
        cmds = [
                "set size 0.8, 1.0",
                "set key tm",
                "set xlabel '$x$'",
                "set ylabel '$y$'",
                ]
        print(f"cmds={cmds}")
    if dts is None:
        x = np.arange(31) * (6 / 30) - 3
        y = np.cos(x)
        yerr = 0.1 / (1 + x**2)
        dts = {
                "table.txt" : azip(x, y, yerr),
                }
        print(f"dts={dts}")
        if lines is None:
            lines = [
                    "plot [-3:3] [-1.5:1.5]",
                    "0 not",
                    "sin(x) w l t '$y = \\sin(x)$'",
                    ]
            if "table.txt" in dts:
                lines.append("'table.txt' w yerrorb t '$y = \\cos(x)$'")
            print(f"lines={lines}")
    elif lines is None:
        lines = [
                "plot [:] [:]",
                ]
        for key, val in dts.items():
            if len(val[0]) >= 3:
                lines.append(f"'{key}' u 1:2:3 w yerrorb t '{key}'")
            elif len(val[0]) == 2:
                lines.append(f"'{key}' u 1:2 w p t '{key}'")
            else:
                lines.append(f"'{key}' t '{key}'")
        print(f"lines={lines}")
    if fn is None:
        print(f"fn={fn}")
        print(f"is_run_make={is_run_make}")
        print(f"is_display={is_display}")
    populate_pyplot_folder(
            path,
            fn = target_fn,
            dict_datatable = dts,
            plot_cmds = cmds,
            plot_lines = lines,
            )
    if is_run_make:
        status = subprocess.run([ "make", "-C", path, ], capture_output = True, text = True)
        if is_verbose or status.returncode != 0:
            print("stdout:")
            print(status.stdout)
            print("stderr:")
            print(status.stderr)
        assert status.returncode == 0
        if target is None:
            path_img = os.path.join(path, "plot-0.png")
        else:
            path_img = target
        if is_display:
            display_img(path_img)
        return path_img
    else:
        assert not is_display
        # return directory that contain the sources instead of the png path
        return path

def plot_view(
        fn = None,
        dts = None,
        cmds = None,
        lines = None,
        ):
    return plot_save(
            fn = fn, dts = dts, cmds = cmds, lines = lines,
            is_run_make = True, is_display = True,
            )
