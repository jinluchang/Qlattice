# Environment variables

- `q_end_time`

  Program finish time in seconds since epoch.

  Used when `check_time_limit()`. Possible setting can be: `export q_end_time="$(($(date +%s) + 12 * 60 * 60))"` for jobs run for at most 12 hours.

  Default is empty.

  `q_time_limit`

  Total running time of program in seconds.

  Used when `check_time_limit()`, but will be overrided if `q_end_time` is set.

  Default is `43200`.

- `q_time_budget`

  Default budget time in seconds. Possible setting can be `export q_budget=$((1 * 60 * 60))`

  Used when `check_time_limit()`.

  Default is `900`.

- `q_field_init`

  Control how field objects' data are initialized.

  Choices are `fast` (default), `zero`, `random`.

- `q_mem_cache_max_size`

  Memory cache size in MB (per processes) for `qlat::vector` allocation.

  Default is `512` MB.

- `q_mem_cache_max_size_acc`

  Memory cache size in MB (per processes) for `qlat::vector` allocation with `mem_type=MemType::Acc`.

  Default is the same as `q_mem_cache_max_size`.

- `q_alloc_mem_max_size`

  Maximum allocated memory size in MB (per processes) for combined `qlat::vector` and `qlat::vector_acc` allocation. Cache size will be reduced when this limit is reached.

  Default is `256 * 1024` MB.

- `q_num_threads`

  Number of OpenMP threads (will be overrided by `OMP_NUM_THREADS`). Possible setting can be `export OMP_NUM_THREADS=16` number should be adjusted by number of cores.

  Default is `2`.

- `q_acc_num_threads`

  Number of `qacc` threads.

  Default is `32`.

- `q_verbose`

  Level of verbosity. Need to be more than `0` for the timing info to be shown automatically. Possible setting can be: `export q_verbose=2`

  Default is `-1`.

- `q_timer_mini_auto_display`

  Minimum time between auto-display of timer information summary.

  Default is `5.0 * 60.0`.

- `q_timer_mini_auto_show`

  Minimum run time for a function for its information to be shown when it start or stop.

  Default is `1.0`.

- `q_timer_max_always_show`

  Maximum number of times to always show function start or stop.

  Default is `2`.

- `q_timer_max_func_name_len`

  Maximum length for a function name before truncation.

  Default is `50`.

- `q_malloc_mmap_threshold`

  In unit of bytes.

  Default is empty. It does not alter the system setting. Possible setting can be `export q_malloc_mmap_threshold=8192`.

- `q_mk_id_node_in_shuffle_seed`

  Seed for initializing `id_node_in_shuffle`.

  Default is `4`. If start with `"seed_"`, then will be random initialization. Otherwise will be viewed as an int and used as `step_size`.

- `q_qar_multi_vol_max_size`

  Maximum size of a `qar` file in bytes. If the total size of the folder is larger, a multi-volume `qar` file will be created.

  Default is `500L * 1000L * 1000L * 1000L` (500 GB). If `q_qar_multi_vol_max_size` is negative, the size is unlimited.

  Note, `qar` never splits a single file into multiple `qar` volume. The limit may be exceeded due to the header size or a single file being too large.

- `q_write_par_limit`

  Default is `16`.

- `q_read_par_limit`

  Default is `16`.

- `q_fftw_plan_flag`

  Default is `estimate`. Possible values include `estimate`, `measure`.

- `q_mpi_alltoallv_type`

  Default is `custom`. Possible values include `native`, `custom`.

- `q_mpi_alltoallv_max_parallel_transfer`

  Default is `8`. Number of maximum parallel transfer in `mpi_alltoallv`.

- `q_default_mem_type`

  Default is `uvm`. Possible values include `cpu`, `acc`, `comm`, `uvm`.

## Useful options

- `OMP_STACKSIZE=8M` OpenMP option for setting per thread stack size.

- `--debug-signals` **Grid option** for intercept some errors.
