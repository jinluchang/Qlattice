# Shape assertion pattern

The codebase follows a consistent pattern of asserting tensor shapes at every stage of computation.

## Core rules

1. **Define expected shape upfront.** Before any computation, derive the expected shape from the lattice shape, batch dims, and other known sizes:
   ```python
   lat_shape = u.shape[:4]
   n_fields = u.shape[-1]
   per_field_shape = (*lat_shape, 3, 3, n_fields)
   full_out_shape = (*lat_shape, 4, 3, 3, m)
   ```

2. **Always assert the whole shape, not slices of it.** Compare the full `.shape` tuple against the full expected tuple. Prefer `assert u.shape == full_shape` over `assert u.shape[4:] == (4, 3, 3)` — a slice assertion can miss a silently changed lattice dimension.

3. **Assert every intermediate result immediately.** After every operation that produces a tensor, assert its shape matches the predefined expected shape. Never defer the check:
   ```python
   contribs = jnp.stack(contribs, axis=0)
   assert contribs.shape == (7, *in_per_field_shape), (
       f"func: contribs.shape={contribs.shape} != (7, {in_per_field_shape})"
   )
   ```

4. **Assert function input shapes at entry.** Validate input parameters and arguments before any computation — use `lat_shape`, `n_fields`, or hard-coded constants:
   ```python
   expected_u_shape = (*lat_shape, 4, 3, 3)
   assert u.shape == expected_u_shape, f"func: u.shape={u.shape} != {expected_u_shape}"
   assert params["conv_kernel"].shape == (n_offsets, 24, 4), ...
   assert target.shape == lat_shape, ...
   ```

5. **Assert final output shape before return.** The last thing before returning is a shape check on the result:
   ```python
   result = jnp.stack(dir_fields, axis=4)
   assert result.shape == full_out_shape, ...
   return result
   ```

6. **Use descriptive error messages** of the form `f"func_name: var.shape={var.shape} != expected"`.

## Rationale

- Catches JIT shape errors at compile time rather than silently producing wrong results.
- Makes the expected tensor contract self-documenting — the expected shape is written right next to the computation.
- JAX's JIT will dead-code eliminate these asserts in compiled code, so there is no runtime overhead in production.
