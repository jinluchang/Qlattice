# Environment variables

- `q_end_time`

  Program finish time in seconds since epoch.

  Used when `check_time_limit()`.

  Default is empty.

  `q_time_limit`

  Total running time of program in seconds.

  Used when `check_time_limit()`, but will be override if `q_end_time` is set.

  Default is `43200`.

- `q_budget`

  Default budget time in seconds.

  Used when `check_time_limit()`.

  Default is `900`.

- `q_field_init`

  Control how field objects' data are initialized.

  Choices are `fast` (default), `zero`, `random`.

- `q_mem_cache_max_size`

  Memory cache size in MB (per processes) for `qlat::vector` allocation.

  Default is `512 MB`.

- `q_num_threads`

  Number of OpenMP threads (will be override by `OMP_NUM_THREADS`).

  Default is `2`.

- `q_acc_num_threads`

  Number of `qacc` threads.

  Default is `32`.

- `q_verbose`

  Level of verbosity. Need to be more than `0` for the timing info to be shown automatically.

  Default is `-1`.

- `q_timer_mini_auto_display`

  Minimum time between auto-display of timer information summary.

  Default is `5.0 * 60.0`.

- `q_timer_mini_auto_show`

  Minimum run time for a function for its information to be shown when it start or stop.

  Default is `1.0`.

- `q_timer_max_always_show`

  Maximum number of times to always show function start or stop.

  Default is `10`.

- `q_timer_max_func_name_len`

  Maximum length for a function name before truncation.

  Default is `50`.

- `q_malloc_mmap_threshold`

  In unit of bytes.

  Default is empty. It does not alter the system setting. Suggested value is `8192`.

- `q_mk_id_node_in_shuffle_seed`

  Seed for initializing `id_node_in_shuffle`.

  Default is `4`. If start with `"seed_"`, then will be random initialization. Otherwise will be viewed as an int and used as `step_size`.

- `q_qar_multi_vol_max_size`

  Maximum size of a `qar` file in bytes. If the total size of the folder is larger, a multi-volume `qar` file will be created.

  Default is `500L * 1000L * 1000L * 1000L` (500 GB). If `q_qar_multi_vol_max_size` is negative, the size is unlimited.

  Note, `qar` never splits a single file into multiple `qar` volume. The limit may be exceeded due to the header size or a single file being too large.

## Useful options

- `OMP_STACKSIZE=8M` OpenMP option for setting per thread stack size.

- `--debug-signals` **Grid option** for intercept some errors.
