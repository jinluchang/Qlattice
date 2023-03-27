#!/usr/bin/env python3

import qlat as q

q.begin_with_mpi()

q.show_machine()

q.displayln_info(f"CHECK: finished successfully.")

q.end_with_mpi()
